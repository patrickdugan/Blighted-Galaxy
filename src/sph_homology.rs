//! BG-PHY-24: radial homology, virial and free-boundary force diagnosis.
//!
//! For fixed particle masses define X_i(s)=C+s(X_i(1)-C), where C is the
//! mass centroid. The same complete cold-EOS + Poisson-matched potential
//! is re-evaluated at each scale, with exact adaptive h/rho and MLS
//! free-surface reconstruction. The derivative is
//!
//!   dV/ds = -sum_i (X_i(1)-C) dot F_i(X(s)).
//!
//! An independent *scalar* energy evaluator, using fresh h/rho and MLS
//! re-solves, supplies central-difference and active-cap-set controls.
//! Shell budgets use only the existing conservative Hamiltonian's
//! pressure and gravitational force decomposition, never new forces.
//!
//! This module does not tune the EOS, hydrostatic reference, or gameplay.
use crate::sph3::{self,Fluid3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState,GravityModel};
use crate::sph_density_audit::{self,DensityEstimator};
use crate::sph_rk_compatible::DensityChoice;
use crate::sph_surface_gravity::{self,CoupledSnapshot};
use crate::white_dwarf::WhiteDwarfProfile;
use crate::physics::G_SI;

fn centroid(p:&[Fluid3])->V3 {
    let mass=p.iter().map(|p|p.mass_kg).sum::<f64>();
    p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/mass
}
pub fn rescale(p:&[Fluid3],scale:f64)->Vec<Fluid3> {
    assert!(scale.is_finite()&&scale>0.0);
    let c=centroid(p);
    p.iter().map(|x|{
        let mut updated=*x;
        updated.pos_m=c+(x.pos_m-c)*scale;
        updated
    }).collect()
}
#[derive(Clone,Debug)]
pub struct ScalarPotential {
    pub cold_internal_j:f64,
    pub gravity_j:f64,
    pub total_j:f64,
    pub cap_mask:Vec<bool>,
    pub minimum_omega:f64,
    pub mls_fallback_fraction:f64,
}
/// Evaluate the scalar functional independently of its analytic gradient.
/// This uses a fresh exact h/rho solve and a fresh BG-PHY-20 MLS fit, with
/// the original fixed Poisson pair potential at the resulting h_i.
pub fn scalar_potential(
    p:&[Fluid3],physics:Sph3Config,cfg:AdaptiveConfig,
    radius:f64,density:DensityChoice
)->ScalarPotential {
    assert!(p.len()>=16&&cfg.gravity_model==GravityModel::MatchedCubicSpline);
    let cfg=cfg.with_newton_density();
    let mut state=AdaptiveState::new(p,radius,cfg);
    sph_adaptive::update_adaptive_densities(p,&mut state,cfg);
    let (field,fallback) = match density {
        DensityChoice::Original => (state.densities_kg_m3.clone(),0.0),
        DensityChoice::SurfaceBlended => {
            let batch=sph_density_audit::estimate(p,&state,cfg.kernel);
            let item=batch.fields.iter().find(|x|
                x.method==DensityEstimator::SurfaceBlendedMls)
                .expect("surface-gated MLS must be available");
            (item.rho.clone(),item.fallback_fraction)
        },
    };
    let internal=p.iter().zip(&field).map(|(x,rho)|
        x.mass_kg*sph3::internal_energy_j_per_kg(*rho,physics))
        .sum::<f64>();
    let gravity=sph_surface_gravity::matched_gravity_potential_j(
        p,&state,cfg);
    let cap_mask=state.smoothing_m.iter().map(|h|
        (*h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8||
        (*h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8).collect::<Vec<_>>();
    ScalarPotential {
        cold_internal_j:internal,gravity_j:gravity,
        total_j:internal+gravity,cap_mask,
        minimum_omega:state.omega.iter().copied()
            .fold(f64::INFINITY,f64::min),
        mls_fallback_fraction:fallback,
    }
}
#[derive(Clone,Debug)]
pub struct HomologyPoint {
    pub scale:f64,
    pub rmax_km:f64,
    pub rmax_bias_pct:f64,
    pub internal_j:f64,
    pub gravity_j:f64,
    pub total_j:f64,
    pub analytic_dv_ds_j:f64,
    pub pressure_dv_ds_j:f64,
    pub gravity_dv_ds_j:f64,
    pub fd_dv_ds_j:f64,
    /// |dV_analytic-dV_FD| / (|dU/ds|+|dPhi/ds|).
    /// This avoids the singular denominator near stationary equilibrium.
    pub derivative_relative_error:f64,
    pub curvature_fd_j:f64,
    pub active_set_stable:bool,
    pub cap_fraction:f64,
    pub minimum_omega:f64,
    pub mls_fallback_fraction:f64,
    pub net_force_ratio:f64,
    pub net_torque_ratio:f64,
}
pub fn point(
    seed:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,density:DensityChoice,scale:f64,epsilon:f64,
)->HomologyPoint {
    assert!((0.30..=1.60).contains(&scale));
    assert!((1e-7..=1e-2).contains(&epsilon));
    let p=rescale(seed,scale);
    let center=centroid(seed);
    let coupled=sph_surface_gravity::snapshot(
        &p,physics,cfg,star.radius_m,density);
    let scalar=scalar_potential(&p,physics,cfg,star.radius_m,density);
    let potential_scale=scalar.cold_internal_j.abs()+
        scalar.gravity_j.abs();
    assert!((scalar.total_j-coupled.total_potential_j).abs()/
        potential_scale.max(1.0)<3e-12,
        "independent total scalar differs from analytic energy");
    let mut dp=0.0;
    let mut dg=0.0;
    for (i,x) in seed.iter().enumerate(){
        let arm=x.pos_m-center;
        dp-=arm.dot(coupled.pressure_force_n[i]);
        dg-=arm.dot(coupled.gravity_force_n[i]);
    }
    let plus=scalar_potential(&rescale(seed,scale+epsilon),
        physics,cfg,star.radius_m,density);
    let minus=scalar_potential(&rescale(seed,scale-epsilon),
        physics,cfg,star.radius_m,density);
    let fd=(plus.total_j-minus.total_j)/(2.0*epsilon);
    let curvature=(plus.total_j-2.0*scalar.total_j+minus.total_j)/
        (epsilon*epsilon);
    let stable=plus.cap_mask==scalar.cap_mask&&
        minus.cap_mask==scalar.cap_mask&&
        (plus.mls_fallback_fraction-scalar.mls_fallback_fraction).abs()<1e-12&&
        (minus.mls_fallback_fraction-scalar.mls_fallback_fraction).abs()<1e-12;
    let error=((dp+dg)-fd).abs()/(dp.abs()+dg.abs()).max(1.0);
    if stable {
        assert!(error<1e-5,
            "homology analytic/FD derivative mismatch: s={scale} error={error}");
    }
    let rmax=p.iter().map(|x|(x.pos_m-centroid(&p)).len())
        .fold(0.0_f64,f64::max);
    HomologyPoint{
        scale,rmax_km:rmax/1000.0,
        rmax_bias_pct:100.0*(rmax/star.radius_m-1.0),
        internal_j:scalar.cold_internal_j,
        gravity_j:scalar.gravity_j,total_j:scalar.total_j,
        analytic_dv_ds_j:dp+dg,pressure_dv_ds_j:dp,
        gravity_dv_ds_j:dg,fd_dv_ds_j:fd,
        derivative_relative_error:error,curvature_fd_j:curvature,
        active_set_stable:stable,
        cap_fraction:scalar.cap_mask.iter().filter(|x|**x).count() as f64/
            seed.len() as f64,
        minimum_omega:scalar.minimum_omega,
        mls_fallback_fraction:scalar.mls_fallback_fraction,
        net_force_ratio:coupled.net_force_ratio,
        net_torque_ratio:coupled.net_torque_ratio,
    }
}
#[derive(Clone,Debug)]
pub struct HomologyMinimum {
    /// A sign change in dV/ds from negative to positive was bracketed.
    pub bracketed:bool,
    pub stationary_scale:Option<f64>,
    pub stationary_rmax_bias_pct:Option<f64>,
    pub stationary_curvature_positive:bool,
    pub stationary_total_j:Option<f64>,
    pub reference_scale:f64,
    pub reference_energy_j:f64,
    pub lower_energy_j:f64,
    pub upper_energy_j:f64,
    pub samples:Vec<HomologyPoint>,
}
pub fn scan(
    seed:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,density:DensityChoice,
    lower:f64,upper:f64,step:f64,epsilon:f64,
)->HomologyMinimum {
    assert!(lower>=0.40&&upper<=1.5&&lower<1.0&&upper>1.0);
    assert!(step>0.0&&step<=0.2&&
        ((upper-lower)/step).round()>=3.0);
    let intervals=((upper-lower)/step).round() as usize;
    assert!(intervals<=60);
    let mut points=Vec::with_capacity(intervals+1);
    for index in 0..=intervals {
        let s=lower+(upper-lower)*index as f64/intervals as f64;
        points.push(point(seed,star,physics,cfg,density,s,epsilon));
    }
    let mut optimum=None;
    let mut energy=f64::INFINITY;
    for pair in points.windows(2){
        if pair[0].analytic_dv_ds_j>=0.0||
            pair[1].analytic_dv_ds_j<=0.0 {continue;}
        let mut lo=pair[0].scale;
        let mut hi=pair[1].scale;
        for _ in 0..14 {
            let mid=0.5*(lo+hi);
            // The root search has an analytic force; each midpoint
            // independently recomputes the complete density/h/MLS state.
            let p=rescale(seed,mid);
            let force=sph_surface_gravity::snapshot(
                &p,physics,cfg,star.radius_m,density);
            let c=centroid(seed);
            let derivative=-seed.iter().enumerate().map(|(i,x)|
                (x.pos_m-c).dot(force.total_force_n[i])).sum::<f64>();
            if derivative<0.0 {lo=mid;} else {hi=mid;}
        }
        let s=0.5*(lo+hi);
        let candidate=point(seed,star,physics,cfg,density,s,epsilon);
        if candidate.curvature_fd_j>0.0&&candidate.total_j<energy{
            energy=candidate.total_j;
            optimum=Some(candidate);
        }
    }
    let ref_point=point(seed,star,physics,cfg,density,1.0,epsilon);
    let r_seed=seed.iter().map(|x|(x.pos_m-centroid(seed)).len())
        .fold(0.0_f64,f64::max);
    let first=points.first().unwrap();
    let last=points.last().unwrap();
    HomologyMinimum{
        bracketed:optimum.is_some(),
        stationary_scale:optimum.as_ref().map(|x|x.scale),
        stationary_rmax_bias_pct:optimum.as_ref().map(|x|
            100.0*(x.scale*r_seed/star.radius_m-1.0)),
        stationary_curvature_positive:optimum.as_ref()
            .map(|x|x.curvature_fd_j>0.0).unwrap_or(false),
        stationary_total_j:optimum.as_ref().map(|x|x.total_j),
        reference_scale:1.0,reference_energy_j:ref_point.total_j,
        lower_energy_j:first.total_j,upper_energy_j:last.total_j,
        samples:points,
    }
}
#[derive(Clone,Copy,Debug,Default)]
pub struct ShellVirial {
    pub index:usize,
    pub count:usize,
    pub mass_fraction_lo:f64,
    pub mass_fraction_hi:f64,
    pub mean_radius_km:f64,
    /// Mass-weighted radial force divided by sum_i m_i g_ref,i.
    pub pressure_support_ratio:f64,
    pub gravity_inward_ratio:f64,
    pub net_radial_ratio:f64,
    pub pressure_tangent_rms:f64,
    pub net_tangent_rms:f64,
    pub pressure_virial_j:f64,
    pub gravity_virial_j:f64,
    pub total_virial_j:f64,
    pub partition_mean:f64,
    pub sp_density_log_bias_mean:f64,
    pub h_capped_fraction:f64,
    pub minimum_omega:f64,
    pub mls_fallback_fraction:f64,
}
/// Decompose forces at *identical frozen positions* into eight equal-mass
/// radial shell diagnostics. The original independent 1D hydrostatic
/// reference supplies expected local inward gravity at these positions.
pub fn shell_virial(
    particles:&[Fluid3],star:&WhiteDwarfProfile,
    physics:Sph3Config,cfg:AdaptiveConfig,density:DensityChoice,
    shells:usize,
)->Vec<ShellVirial> {
    let n=particles.len();
    assert!(shells>=2&&shells<=n/4&&n%shells==0);
    let center=centroid(particles);
    let cfg=cfg.with_newton_density();
    let coupled=sph_surface_gravity::snapshot(
        particles,physics,cfg,star.radius_m,density);
    let mut state=AdaptiveState::new(particles,star.radius_m,cfg);
    sph_adaptive::update_adaptive_densities(particles,&mut state,cfg);
    let batch=sph_density_audit::estimate(particles,&state,cfg.kernel);
    let fallback=batch.fields.iter().find(|x|
        x.method==DensityEstimator::SurfaceBlendedMls)
        .unwrap().fallback_fraction;
    let mut rank=(0..n).collect::<Vec<_>>();
    rank.sort_by(|&i,&j| {
        (particles[i].pos_m-center).len_sq().total_cmp(
            &(particles[j].pos_m-center).len_sq()).then(i.cmp(&j))
    });
    let per=n/shells;
    let mut output=Vec::with_capacity(shells);
    for shell in 0..shells {
        let indices=&rank[shell*per..(shell+1)*per];
        let (mut mass,mut radius,mut support,mut grav,mut total,
            mut ref_force,mut ref_squared)=(0.0,0.0,0.0,0.0,
                0.0,0.0,0.0);
        let (mut p_tan,mut net_tan,mut virial_p,mut virial_g,
            mut partition,mut log_bias)=(0.0,0.0,0.0,0.0,0.0,0.0);
        let mut caps=0;
        let mut min_omega=f64::INFINITY;
        for &i in indices {
            let x=&particles[i];
            let rvec=x.pos_m-center;
            let r=rvec.len();
            let u=rvec*(1.0/r);
            let m=x.mass_kg;
            let a_p=coupled.pressure_force_n[i]/m;
            let a_g=coupled.gravity_force_n[i]/m;
            let a_total=a_p+a_g;
            let g=G_SI*star.enclosed_mass_at(r)/(r*r);
            mass+=m;radius+=m*r;
            support+=m*a_p.dot(u);
            grav-=m*a_g.dot(u);
            total+=m*a_total.dot(u);
            ref_force+=m*g;
            ref_squared+=m*g*g;
            p_tan+=m*(a_p-u*a_p.dot(u)).len_sq();
            net_tan+=m*(a_total-u*a_total.dot(u)).len_sq();
            virial_p+=rvec.dot(coupled.pressure_force_n[i]);
            virial_g+=rvec.dot(coupled.gravity_force_n[i]);
            partition+=m*batch.partition[i];
            let reference=star.density_at(r);
            if reference>0.0 {
                log_bias+=m*(state.densities_kg_m3[i]/reference).ln();
            }
            let h=state.smoothing_m[i];
            if (h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8||
               (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8 {caps+=1;}
            min_omega=min_omega.min(state.omega[i]);
        }
        let norm=ref_force.abs().max(1e-100);
        output.push(ShellVirial{
            index:shell,count:per,mass_fraction_lo:shell as f64/shells as f64,
            mass_fraction_hi:(shell+1) as f64/shells as f64,
            mean_radius_km:radius/mass/1000.0,
            pressure_support_ratio:support/norm,
            gravity_inward_ratio:grav/norm,
            net_radial_ratio:total/norm,
            pressure_tangent_rms:(p_tan/ref_squared.max(1e-100)).sqrt(),
            net_tangent_rms:(net_tan/ref_squared.max(1e-100)).sqrt(),
            pressure_virial_j:virial_p,
            gravity_virial_j:virial_g,
            total_virial_j:virial_p+virial_g,
            partition_mean:partition/mass,
            sp_density_log_bias_mean:log_bias/mass,
            h_capped_fraction:caps as f64/per as f64,
            minimum_omega:min_omega,
            mls_fallback_fraction:fallback,
        });
    }
    output
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn setup() ->(Vec<Fluid3>,WhiteDwarfProfile,Sph3Config,AdaptiveConfig){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),32,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        (p,star,Sph3Config::cold_white_dwarf(2.0).unwrap(),cfg)
    }
    #[test]
    fn unit_scale_and_translation_preserve_full_potential(){
        let (p,star,phys,cfg)=setup();
        let a=scalar_potential(&p,phys,cfg,star.radius_m,
            DensityChoice::SurfaceBlended);
        let b=scalar_potential(&rescale(&p,1.0),phys,cfg,star.radius_m,
            DensityChoice::SurfaceBlended);
        assert!((a.total_j-b.total_j).abs()/a.total_j.abs()<1e-12);
        let mut translated=p.to_vec();
        for x in &mut translated {
            x.pos_m+=V3::new(-2.0e8,3.0e8,1.0e8);
        }
        let moved=scalar_potential(&translated,phys,cfg,star.radius_m,
            DensityChoice::SurfaceBlended);
        assert!((a.total_j-moved.total_j).abs()/a.total_j.abs()<1e-11);
    }
    #[test]
    fn virial_matches_independently_resolved_scalar_homology_derivative(){
        let (p,star,phys,cfg)=setup();
        for mode in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
            let pt=point(&p,&star,phys,cfg,mode,1.0,1e-4);
            assert!(pt.derivative_relative_error<1e-5);
            assert!((pt.analytic_dv_ds_j-
                pt.pressure_dv_ds_j-pt.gravity_dv_ds_j).abs()/
                (pt.pressure_dv_ds_j.abs()+pt.gravity_dv_ds_j.abs())<1e-12);
            let b=shell_virial(&p,&star,phys,cfg,mode,8);
            assert_eq!(b.len(),8);
            let sum=b.iter().map(|x|x.total_virial_j).sum::<f64>();
            assert!((sum+pt.analytic_dv_ds_j).abs()/
                (pt.pressure_dv_ds_j.abs()+pt.gravity_dv_ds_j.abs())<1e-10);
        }
    }
}
