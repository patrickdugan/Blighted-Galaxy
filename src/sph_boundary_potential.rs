//! BG-PHY-25: independently re-solved scalar cold-EOS plus matched-gravity
//! potential for candidate boundary enthalpy density fields.
//!
//! Scalar energy is a diagnostic only. The candidate pressure force is
//! measured via central differences of the SAME physical cold-EOS energy,
//! including the implicit h root, volumes and ghost-boundary geometry.
//! The method never inserts an unpaired pressure force into dynamics.
//! This provides a conservative Hamiltonian target if a candidate succeeds.
use crate::physics::G_SI;
use crate::sph3::{self,Fluid3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use crate::sph_kernel;
use crate::sph_free_boundary::{self,EnthalpyCandidate,BoundaryGeometry};
use crate::sph_density_audit::{self,DensityEstimator};
use crate::sph_surface_gravity;
use crate::sph_quadratic_audit;
use crate::white_dwarf::WhiteDwarfProfile;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Candidate {
    Original,
    SurfaceMls,
    EnthalpyAffineBlend,
    VacuumGhostBlend,
    MomentRatioGhostBlend,
}
impl Candidate {
    pub fn label(self)->&'static str{
        match self{
            Self::Original=>"original",
            Self::SurfaceMls=>"surface_mls",
            Self::EnthalpyAffineBlend=>"enthalpy_affine_blend",
            Self::VacuumGhostBlend=>"vacuum_ghost_blend",
            Self::MomentRatioGhostBlend=>"moment_ratio_ghost_blend",
        }
    }
    pub fn all()->[Self;5]{
        [Self::Original,Self::SurfaceMls,Self::EnthalpyAffineBlend,
            Self::VacuumGhostBlend,Self::MomentRatioGhostBlend]
    }
}
fn centroid(p:&[Fluid3])->V3{
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/m
}
fn gradient_mode(mode:Candidate)->Option<crate::sph_rk_compatible::DensityChoice>{
    match mode{
        Candidate::Original=>Some(crate::sph_rk_compatible::DensityChoice::Original),
        Candidate::SurfaceMls=>Some(crate::sph_rk_compatible::DensityChoice::SurfaceBlended),
        _=>None,
    }
}
#[derive(Clone,Debug)]
pub struct ScalarState{
    pub rho:Vec<f64>,
    pub geometry:Vec<BoundaryGeometry>,
    pub scalar_internal_j:f64,
    pub matched_gravity_j:f64,
    pub total_j:f64,
    pub fallback_fraction:f64,
    pub minimum_omega:f64,
    pub h_cap_fraction:f64,
    pub density_global_log_rms:f64,
    pub density_outer_log_rms:f64,
    pub density_last_log_rms:f64,
    pub enthalpy_gradient_rms:f64,
    pub enthalpy_outer_gradient_rms:f64,
    pub last_shell_partition:f64,
    pub last_shell_normal_alignment:f64,
    pub last_shell_moment_ratio:f64,
    pub last_shell_partition_depth_h:f64,
    pub last_shell_moment_depth_h:f64,
    pub last_shell_partition_below_half_fraction:f64,
    /// Independent 1D stellar radius used only for this score.
    pub last_shell_true_reference_depth_h:f64,
    pub last_shell_partition_depth_abs_error_h:f64,
    pub last_shell_moment_depth_abs_error_h:f64,
}
/// A full state evaluated from scratch for the supplied geometry.
/// The reference density and gravity are used *only for scoring*.
pub fn state(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,with_diagnostics:bool,
)->ScalarState {
    assert!(p.len()>=16&&physics.valid()&&cfg.valid());
    let cfg=cfg.with_newton_density();
    let mut adaptive=AdaptiveState::new(p,star.radius_m,cfg);
    sph_adaptive::update_adaptive_densities(p,&mut adaptive,cfg);
    // The original and BG-PHY-20 scalars do not require ghost geometry
    // during a scalar-energy perturbation. Reserve the expensive boundary
    // moment calculation for new candidates and diagnostic snapshots.
    let (rho,fallback,geometry)=match mode{
        Candidate::Original=>(
            adaptive.densities_kg_m3.clone(),0.0,Vec::new()),
        Candidate::SurfaceMls=>{
            let batch=sph_density_audit::estimate(p,&adaptive,cfg.kernel);
            let field=batch.fields.iter().find(|x|
                x.method==DensityEstimator::SurfaceBlendedMls).unwrap();
            (field.rho.clone(),field.fallback_fraction,Vec::new())
        }
        Candidate::EnthalpyAffineBlend|
        Candidate::VacuumGhostBlend|
        Candidate::MomentRatioGhostBlend=>{
            let batch=sph_free_boundary::reconstruct(
                p,&adaptive,cfg.kernel,star.eos);
            let wanted=match mode{
                Candidate::EnthalpyAffineBlend=>EnthalpyCandidate::BlendedH,
                Candidate::VacuumGhostBlend=>EnthalpyCandidate::GhostBlendedH,
                Candidate::MomentRatioGhostBlend=>
                    EnthalpyCandidate::MomentRatioGhostBlend,
                _=>unreachable!(),
            };
            let field=batch.fields.iter().find(|x|
                x.mode==wanted).unwrap();
            (field.rho.clone(),field.fallback_fraction,batch.geometry)
        }
    };
    let internal=p.iter().zip(&rho).map(|(x,r)|
        x.mass_kg*sph3::internal_energy_j_per_kg(*r,physics))
        .sum::<f64>();
    let gravity=sph_surface_gravity::matched_gravity_potential_j(
        p,&adaptive,cfg);
    let n=p.len();
    let cap=adaptive.smoothing_m.iter().filter(|&&h|
        (h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8||
        (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8).count();
    let mut result=ScalarState{
        rho,geometry,scalar_internal_j:internal,
        matched_gravity_j:gravity,total_j:internal+gravity,
        fallback_fraction:fallback,
        minimum_omega:adaptive.omega.iter().copied().fold(
            f64::INFINITY,f64::min),
        h_cap_fraction:cap as f64/n as f64,
        density_global_log_rms:0.0,density_outer_log_rms:0.0,
        density_last_log_rms:0.0,enthalpy_gradient_rms:0.0,
        enthalpy_outer_gradient_rms:0.0,last_shell_partition:0.0,
        last_shell_normal_alignment:0.0,
        last_shell_moment_ratio:0.0,
        last_shell_partition_depth_h:0.0,
        last_shell_moment_depth_h:0.0,
        last_shell_partition_below_half_fraction:0.0,
        last_shell_true_reference_depth_h:0.0,
        last_shell_partition_depth_abs_error_h:0.0,
        last_shell_moment_depth_abs_error_h:0.0,
    };
    if with_diagnostics {
        if result.geometry.is_empty() {
            result.geometry=sph_free_boundary::reconstruct(
                p,&adaptive,cfg.kernel,star.eos).geometry;
        }
        let c=centroid(p);
        let mut rank=(0..n).collect::<Vec<_>>();
        rank.sort_by(|&i,&j|(p[i].pos_m-c).len_sq()
            .total_cmp(&(p[j].pos_m-c).len_sq()).then(i.cmp(&j)));
        let q=n/4;
        let last=n/8;
        let mut full=0.0;
        let mut outer=0.0;
        let mut final_shell=0.0;
        let mut mass=0.0;
        let mut outer_mass=0.0;
        let mut last_mass=0.0;
        for (index,&i) in rank.iter().enumerate(){
            let r=(p[i].pos_m-c).len();
            let ref_rho=star.density_at(r);
            if ref_rho>0.0 {
                let squared=(result.rho[i]/ref_rho).ln().powi(2);
                full+=p[i].mass_kg*squared;
                mass+=p[i].mass_kg;
                if index>=n-q{
                    outer+=p[i].mass_kg*squared;
                    outer_mass+=p[i].mass_kg;
                }
                if index>=n-last{
                    final_shell+=p[i].mass_kg*squared;
                    last_mass+=p[i].mass_kg;
                }
            }
        }
        result.density_global_log_rms=(full/mass.max(1e-99)).sqrt();
        result.density_outer_log_rms=(outer/outer_mass.max(1e-99)).sqrt();
        result.density_last_log_rms=(final_shell/last_mass.max(1e-99)).sqrt();
        let grad=sph_quadratic_audit::audit_with_density_samples(
            p,&adaptive,&result.rho,cfg.kernel,star);
        result.enthalpy_gradient_rms=grad.global.vector_error[5];
        result.enthalpy_outer_gradient_rms=grad.outer.vector_error[5];
        let mut part=0.0;
        let mut align=0.0;
        let (mut ratio,mut depth_s,mut depth_ratio,mut below_half)=
            (0.0,0.0,0.0,0.0);
        let (mut reference_depth,mut s_depth_error,mut ratio_depth_error)=
            (0.0,0.0,0.0);
        for &i in &rank[n-last..] {
            let radial=(p[i].pos_m-c).normalized();
            let g=result.geometry[i];
            part+=g.support;
            align+=g.outward.dot(radial);
            ratio+=g.first_to_zeroth_ratio;
            depth_s+=g.plane_distance_h;
            depth_ratio+=g.ratio_plane_distance_h;
            if g.support<0.5 {below_half+=1.0;}
            let radial_distance=(p[i].pos_m-c).len();
            let reference=(star.radius_m-radial_distance)/
                adaptive.smoothing_m[i];
            reference_depth+=reference;
            s_depth_error+=(g.plane_distance_h-reference).abs();
            ratio_depth_error+=
                (g.ratio_plane_distance_h-reference).abs();
        }
        result.last_shell_partition=part/last as f64;
        result.last_shell_normal_alignment=align/last as f64;
        result.last_shell_moment_ratio=ratio/last as f64;
        result.last_shell_partition_depth_h=depth_s/last as f64;
        result.last_shell_moment_depth_h=depth_ratio/last as f64;
        result.last_shell_partition_below_half_fraction=
            below_half/last as f64;
        result.last_shell_true_reference_depth_h=
            reference_depth/last as f64;
        result.last_shell_partition_depth_abs_error_h=
            s_depth_error/last as f64;
        result.last_shell_moment_depth_abs_error_h=
            ratio_depth_error/last as f64;
    }
    result
}
pub fn energy(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,
)->(f64,f64){
    let result=state(p,star,physics,cfg,mode,false);
    (result.scalar_internal_j,result.matched_gravity_j)
}
pub fn homology_geometry(p:&[Fluid3],s:f64)->Vec<Fluid3>{
    assert!(s>0.0&&s.is_finite());
    let c=centroid(p);
    p.iter().map(|x|{
        let mut y=*x;
        y.pos_m=c+(x.pos_m-c)*s;
        y
    }).collect()
}
#[derive(Clone,Debug)]
pub struct HomologyPoint{
    pub s:f64,
    pub internal_j:f64,
    pub gravity_j:f64,
    pub total_j:f64,
    pub derivative_j:f64,
    pub curvature_j:f64,
}
pub fn homology_point(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,s:f64,eps:f64,
)->HomologyPoint {
    assert!(s>eps&&eps>0.0&&eps<0.01);
    let (u,g)=energy(&homology_geometry(p,s),star,physics,cfg,mode);
    let plus=energy(&homology_geometry(p,s+eps),star,physics,cfg,mode);
    let minus=energy(&homology_geometry(p,s-eps),star,physics,cfg,mode);
    let sum=u+g;
    HomologyPoint{
        s,internal_j:u,gravity_j:g,total_j:sum,
        derivative_j:((plus.0+plus.1)-(minus.0+minus.1))/(2.0*eps),
        curvature_j:((plus.0+plus.1)-2.0*sum+
            (minus.0+minus.1))/(eps*eps),
    }
}
#[derive(Clone,Debug)]
pub struct HomologyOptimum{
    pub bracketed:bool,
    pub s_star:Option<f64>,
    pub radius_bias_pct:Option<f64>,
    pub curvature_positive:bool,
    pub scan:Vec<HomologyPoint>,
}
pub fn homology_minimum(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,
    lower:f64,upper:f64,intervals:usize,eps:f64,
)->HomologyOptimum {
    assert!(lower>0.4&&upper<=1.5&&upper>lower);
    assert!(intervals>=3&&intervals<=24);
    let samples=(0..=intervals).map(|i|{
        let scale=lower+(upper-lower)*i as f64/intervals as f64;
        homology_point(p,star,physics,cfg,mode,scale,eps)
    }).collect::<Vec<_>>();
    let mut best=None;
    let mut energy_min=f64::INFINITY;
    for pair in samples.windows(2){
        if pair[0].derivative_j>=0.0||pair[1].derivative_j<=0.0{
            continue;
        }
        let (mut lo,mut hi)=(pair[0].s,pair[1].s);
        for _ in 0..10 {
            let mid=0.5*(lo+hi);
            let d=homology_point(p,star,physics,cfg,mode,mid,eps).derivative_j;
            if d<0.0 {lo=mid;} else {hi=mid;}
        }
        let s=0.5*(lo+hi);
        let q=homology_point(p,star,physics,cfg,mode,s,eps);
        if q.curvature_j>0.0&&q.total_j<energy_min{
            energy_min=q.total_j;best=Some(q);
        }
    }
    let c=centroid(p);
    let rmax=p.iter().map(|x|(x.pos_m-c).len())
        .fold(0.0_f64,f64::max);
    HomologyOptimum {
        bracketed:best.is_some(),
        s_star:best.as_ref().map(|x|x.s),
        radius_bias_pct:best.as_ref().map(|x|
            100.0*(x.s*rmax/star.radius_m-1.0)),
        curvature_positive:best.as_ref().map(|x|x.curvature_j>0.0)
            .unwrap_or(false),
        scan:samples,
    }
}
/// Variational scalar finite-difference forces along shell radial unit
/// vectors, evaluated with the complete density and h re-solve.
#[derive(Clone,Copy,Debug)]
pub struct ShellPressure{
    pub shell:usize,
    pub pressure_ratio:f64,
    pub gravity_inward_ratio:f64,
    pub net_ratio:f64,
    pub pressure_virial_j:f64,
    pub gravity_virial_j:f64,
    pub mean_radius_km:f64,
}
pub fn shell_radial_forces(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,shells:usize,eps:f64,
)->Vec<ShellPressure>{
    let n=p.len();
    assert!(n%shells==0&&shells>=2&&shells<=n/4);
    let c=centroid(p);
    let mut rank=(0..n).collect::<Vec<_>>();
    rank.sort_by(|&i,&j|(p[i].pos_m-c).len_sq()
        .total_cmp(&(p[j].pos_m-c).len_sq()).then(i.cmp(&j)));
    let scale=eps*star.radius_m;
    let mut output=Vec::new();
    for k in 0..shells{
        let q=&rank[k*n/shells..(k+1)*n/shells];
        let mut plus=p.to_vec();
        let mut minus=p.to_vec();
        let (mut ref_radial,mut average_r)=(0.0,0.0);
        for &i in q {
            let d=p[i].pos_m-c;
            let unit=d.normalized();
            plus[i].pos_m+=unit*scale;
            minus[i].pos_m+=unit*(-scale);
            let r=d.len();
            let gravity=G_SI*star.enclosed_mass_at(r)/(r*r);
            ref_radial+=p[i].mass_kg*gravity;
            average_r+=r;
        }
        let (up,gp)=energy(&plus,star,physics,cfg,mode);
        let (um,gm)=energy(&minus,star,physics,cfg,mode);
        let dp=-(up-um)/(2.0*scale);
        let dg=-(gp-gm)/(2.0*scale);
        // A second probe displaces only this shell proportionally
        // to its physical radius, measuring the exact virial.
        let mut plus_v=p.to_vec();
        let mut minus_v=p.to_vec();
        for &i in q {
            let d=p[i].pos_m-c;
            plus_v[i].pos_m+=d*eps;
            minus_v[i].pos_m-=d*eps;
        }
        let (vup,vgp)=energy(&plus_v,star,physics,cfg,mode);
        let (vum,vgm)=energy(&minus_v,star,physics,cfg,mode);
        let virial_p=-(vup-vum)/(2.0*eps);
        let virial_g=-(vgp-vgm)/(2.0*eps);
        output.push(ShellPressure{
            shell:k,pressure_ratio:dp/ref_radial,
            gravity_inward_ratio:-dg/ref_radial,
            net_ratio:(dp+dg)/ref_radial,
            pressure_virial_j:virial_p,
            gravity_virial_j:virial_g,
            mean_radius_km:average_r/q.len() as f64/1000.0,
        });
    }
    output
}
/// Explicit analytic shell-force comparison for the two already
/// validated BG-PHY-24 density choices. New fields remain FD diagnostics.
pub fn validated_control_error(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,mode:Candidate,eps:f64,
)->Option<f64>{
    let prior=gradient_mode(mode)?;
    let baseline=crate::sph_homology::shell_virial(
        p,star,physics,cfg,prior,8);
    let measured=shell_radial_forces(p,star,physics,cfg,mode,8,eps);
    let mut maximum=0.0_f64;
    for (a,b) in baseline.iter().zip(measured){
        let denominator=a.pressure_support_ratio.abs()+
            a.gravity_inward_ratio.abs();
        let error=(a.pressure_support_ratio-b.pressure_ratio).abs()/
            denominator.max(1e-10);
        maximum=maximum.max(error);
        let error_g=(a.gravity_inward_ratio-b.gravity_inward_ratio).abs()/
            denominator.max(1e-10);
        maximum=maximum.max(error_g);
    }
    Some(maximum)
}
#[cfg(test)]
mod tests{
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    #[test]
    fn scalar_original_density_matches_independent_prior_hamiltonian(){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),32,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        for mode in [Candidate::Original,Candidate::SurfaceMls]{
            let q=state(&p,&star,Sph3Config::cold_white_dwarf(2.0).unwrap(),
                cfg,mode,false);
            let reference=sph_surface_gravity::independent_potential_j(
                &p,Sph3Config::cold_white_dwarf(2.0).unwrap(),cfg,
                star.radius_m,gradient_mode(mode).unwrap());
            assert!((q.total_j-reference).abs()/reference.abs()<1e-12);
        }
    }
    #[test]
    fn scalar_shell_pressure_matches_analytic_old_force(){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),32,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        let phys=Sph3Config::cold_white_dwarf(2.0).unwrap();
        for mode in [Candidate::Original,Candidate::SurfaceMls]{
            let error=validated_control_error(&p,&star,phys,cfg,mode,1e-4)
                .unwrap();
            assert!(error<2e-4,"shell finite-difference mismatch: {error}");
        }
    }
}
