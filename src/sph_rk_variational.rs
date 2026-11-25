//! BG-PHY-21B: expensive TOTAL discrete barotropic-energy gradient oracle.
//!
//! Define U(x)=sum_i m_i u_Fermi(rho_i^*(x)), with rho^* either the
//! exact-root SPH density or the BG-PHY-20 support-gated affine MLS density.
//! Every energy evaluation completely recomputes the implicit h/rho roots,
//! baseline volumes, partition-of-unity S, MLS normal matrices, blending,
//! and cold EOS. Central differences of U with respect to every physical
//! particle coordinate therefore include their TOTAL dependence up to
//! truncation error and fallback/cap non-smoothness.
//!
//! The resulting -dU/dx is a finite-difference Hamiltonian reference.
//! It is expensive O(N coordinate probes × whole density solve), so N<=32.
//! It is a research oracle, not a scalable analytic CRKSPH algorithm.
//! The existing game/stellar integrators and gravity remain untouched.
use crate::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use crate::sph_density_audit::{self,DensityEstimator};
use crate::sph_rk_compatible::{self,DensityChoice};

fn cross(a:V3,b:V3)->V3{
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn component(v:V3,k:usize)->f64{
    match k{0=>v.x,1=>v.y,2=>v.z,_=>panic!("bad axis")}
}
fn displaced(v:V3,k:usize,d:f64)->V3{
    match k{
        0=>V3::new(v.x+d,v.y,v.z),
        1=>V3::new(v.x,v.y+d,v.z),
        2=>V3::new(v.x,v.y,v.z+d),
        _=>panic!("bad axis"),
    }
}
fn selected_density<'a>(
    batch:&'a sph_density_audit::DensityBatch,mode:DensityChoice,
)->&'a [f64]{
    let name=match mode {
        DensityChoice::Original=>DensityEstimator::Original,
        DensityChoice::SurfaceBlended=>DensityEstimator::SurfaceBlendedMls,
    };
    &batch.fields.iter().find(|x|x.method==name)
        .expect("density estimator must exist").rho
}
/// This is the full scalar discrete cold-barotropic thermal energy.
/// The original adaptive-h SPH state is recomputed from scratch for
/// every perturbed particle configuration, including MLS and h caps.
pub fn internal_energy_j(
    p:&[Fluid3],physics:Sph3Config,
    cfg:AdaptiveConfig,reference_radius_m:f64,mode:DensityChoice,
)->f64{
    assert!(p.len()>=16&&p.len()<=32);
    let mut state=AdaptiveState::new(p,reference_radius_m,
        cfg.with_newton_density());
    sph_adaptive::update_adaptive_densities(p,&mut state,
        cfg.with_newton_density());
    let batch=sph_density_audit::estimate(p,&state,cfg.kernel);
    let rho=selected_density(&batch,mode);
    let energy=p.iter().zip(rho).map(|(x,r)|x.mass_kg*
        sph3::internal_energy_j_per_kg(*r,physics)).sum::<f64>();
    assert!(energy.is_finite()&&energy>0.0);
    energy
}
#[derive(Clone,Debug)]
pub struct VariationalGradient {
    pub potential_j:f64,
    pub forces_n:Vec<V3>,
    pub finite_difference_step_m:f64,
    /// Unprojected net force/total absolute force; translation invariance.
    pub net_force_ratio:f64,
    /// Unprojected torque/total force moment; rotation invariance.
    pub net_torque_ratio:f64,
    /// Relative difference from the pairwise compatible-work CRK pressure force.
    pub rk_discrepancy_rms:f64,
    /// Existing exact-root, variable-h SPH pressure formula against -dU/dx.
    /// This should approach zero when mode=Original, and serves as the
    /// independent thermodynamic-control result for the oracle.
    pub conventional_pressure_discrepancy_rms:f64,
}
pub fn total_gradient(
    p:&[Fluid3],physics:Sph3Config,
    cfg:AdaptiveConfig,reference_radius_m:f64,
    mode:DensityChoice,step_fraction:f64,
)->VariationalGradient {
    assert!((1e-7..=1e-2).contains(&step_fraction)&&
        reference_radius_m>0.0);
    let eps=reference_radius_m*step_fraction;
    let baseline=internal_energy_j(
        p,physics,cfg,reference_radius_m,mode);
    let mut perturbed=p.to_vec();
    let mut forces=vec![V3::default();p.len()];
    for i in 0..p.len(){
        let x=p[i].pos_m;
        let mut force=[0.0_f64;3];
        for axis in 0..3{
            perturbed[i].pos_m=displaced(x,axis,eps);
            let plus=internal_energy_j(
                &perturbed,physics,cfg,reference_radius_m,mode);
            perturbed[i].pos_m=displaced(x,axis,-eps);
            let minus=internal_energy_j(
                &perturbed,physics,cfg,reference_radius_m,mode);
            perturbed[i].pos_m=x;
            force[axis]=-(plus-minus)/(2.0*eps);
        }
        forces[i]=V3::new(force[0],force[1],force[2]);
    }
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let com=p.iter().fold(V3::default(),
        |sum,x|sum+x.pos_m*x.mass_kg)/mass;
    let mut total_force=V3::default();
    let mut total_torque=V3::default();
    let mut force_scale=0.0;
    let mut torque_scale=0.0;
    for (x,f) in p.iter().zip(&forces){
        let arm=x.pos_m-com;
        total_force+=*f;
        total_torque+=cross(arm,*f);
        force_scale+=f.len();
        torque_scale+=arm.len()*f.len();
    }
    let mut state=AdaptiveState::new(p,reference_radius_m,cfg);
    let exact_cfg=cfg.with_newton_density();
    let rk=sph_rk_compatible::pressure_snapshot(
        p,&mut state,physics,exact_cfg,mode);
    let mut conventional_state=AdaptiveState::new(
        p,reference_radius_m,exact_cfg);
    let central=Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let conventional=sph_adaptive::forces(
        &central,p,&mut conventional_state,physics,exact_cfg);
    let (mut diff,mut conventional_diff,mut denom)=(0.0,0.0,0.0);
    for (i,f) in forces.iter().enumerate(){
        let other=rk.accelerations_m_s2[i]*p[i].mass_kg;
        let old=conventional.pressure_accelerations_m_s2[i]*
            p[i].mass_kg;
        diff+=(*f-other).len_sq();
        conventional_diff+=(*f-old).len_sq();
        denom+=f.len_sq();
    }
    VariationalGradient{
        potential_j:baseline,forces_n:forces,
        finite_difference_step_m:eps,
        net_force_ratio:total_force.len()/force_scale.max(1.0),
        net_torque_ratio:total_torque.len()/torque_scale.max(1.0),
        rk_discrepancy_rms:(diff/denom.max(1.0)).sqrt(),
        conventional_pressure_discrepancy_rms:
            (conventional_diff/denom.max(1.0)).sqrt(),
    }
}
/// Separate research Hamiltonian integrator. For each force evaluation
/// it recomputes all density derivatives numerically. Gravity absent.
/// This makes cold-u_EOS consistency exact by construction, with
/// step-size-dependent numerical energy error.
#[derive(Clone)]
pub struct VariationalPilot {
    pub particles:Vec<Fluid3>,
    pub physics:Sph3Config,
    pub cfg:AdaptiveConfig,
    pub reference_radius_m:f64,
    pub mode:DensityChoice,
    pub difference_fraction:f64,
}
impl VariationalPilot {
    pub fn gradient(&self)->VariationalGradient {
        total_gradient(&self.particles,self.physics,self.cfg,
            self.reference_radius_m,self.mode,self.difference_fraction)
    }
    pub fn energy(&self)->f64{
        let kinetic=self.particles.iter()
            .map(|x|0.5*x.mass_kg*x.vel_m_s.len_sq()).sum::<f64>();
        kinetic+internal_energy_j(
            &self.particles,self.physics,self.cfg,
            self.reference_radius_m,self.mode)
    }
    pub fn advance(&mut self,dt:f64){
        assert!(dt.is_finite()&&dt.abs()<=0.02);
        if dt==0.0{return;}
        let before=self.gradient();
        for (p,force) in self.particles.iter_mut().zip(before.forces_n){
            p.vel_m_s+=force*(0.5*dt/p.mass_kg);
            p.pos_m+=p.vel_m_s*dt;
        }
        let after=self.gradient();
        for (p,force) in self.particles.iter_mut().zip(after.forces_n){
            p.vel_m_s+=force*(0.5*dt/p.mass_kg);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn fixture()->VariationalPilot{
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1.0e8,0.0,0.0),V3::default(),16,17).0;
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        cfg.eta=1.7;
        VariationalPilot{
            particles:p,physics:Sph3Config::cold_white_dwarf(2.0).unwrap(),
            cfg,reference_radius_m:star.radius_m,
            mode:DensityChoice::SurfaceBlended,difference_fraction:1e-4,
        }
    }
    #[test]
    fn thermodynamic_reference_is_rigid_motion_invariant(){
        let fixture=fixture();
        let e0=fixture.energy();
        let mut shifted=fixture.clone();
        for p in &mut shifted.particles{
            p.pos_m+=V3::new(3.0e8,-4.0e8,2.0e8);
        }
        let e1=shifted.energy();
        assert!((e1-e0).abs()/e0<1e-10);
        let mut rotated=fixture.clone();
        for p in &mut rotated.particles{
            let x=p.pos_m;
            p.pos_m=V3::new(-x.y,x.x,x.z);
        }
        let e2=rotated.energy();
        assert!((e2-e0).abs()/e0<1e-10);
    }
    #[test]
    fn total_derivative_force_is_translation_rotation_neutral(){
        let fixture=fixture();
        let gradient=fixture.gradient();
        assert!(gradient.net_force_ratio<2e-3,
            "translation symmetry discrepancy {}",gradient.net_force_ratio);
        assert!(gradient.net_torque_ratio<2e-3,
            "rotation symmetry discrepancy {}",gradient.net_torque_ratio);
        assert!(gradient.rk_discrepancy_rms.is_finite());
        assert!(gradient.conventional_pressure_discrepancy_rms.is_finite());
        assert!(gradient.forces_n.iter().all(|v|v.finite()));
    }
    #[test]
    fn original_density_agrees_with_variational_sph_pressure_force(){
        let original=fixture();
        for kernel in [crate::sph_kernel::DensityKernel::CubicM4,
                       crate::sph_kernel::DensityKernel::WendlandC4] {
            let mut control=original.clone();
            control.mode=DensityChoice::Original;
            control.cfg.kernel=kernel;
            let var=control.gradient();
            assert!(var.conventional_pressure_discrepancy_rms<1e-6,
                "ordinary variable-h SPH must match -dU/dx: kernel {:?}, residual {}",
                kernel,var.conventional_pressure_discrepancy_rms);
        }
    }
    #[test]
    fn short_variational_orbit_has_bounded_energy_error(){
        let mut fixture=fixture();
        let before=fixture.energy();
        for _ in 0..2{fixture.advance(0.002);}
        let after=fixture.energy();
        let drift=(after-before).abs()/before.abs();
        assert!(drift<2e-4,"variational finite difference drift {drift:e}");
    }
}
