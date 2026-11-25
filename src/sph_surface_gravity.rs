//! BG-PHY-23: opt-in self-gravitating cold-EOS Hamiltonian.
//!
//! H(x,v) = sum_i (m_i |v_i|^2/2 + m_i u_EOS(rho_hat_i(x)))
//!          + G/2 sum_{i<j} m_i m_j [phi(r_ij,h_i)+phi(r_ij,h_j)].
//!
//! rho_hat_i is either the original exact-root SPH density or the
//! BG-PHY-20 partition-gated affine MLS density. The pressure force
//! is BG-PHY-22's fully analytic -dU/dx, with implicit h/rho and MLS.
//! The gravity force is the existing matched Poisson pair force plus
//! its variable-h zeta/Omega correction in sph_adaptive::forces.
//! The potential is evaluated separately from the same matched phi,
//! and independently checked by full coordinate finite differences.
//!
//! Research only: no central massive body, rotation/relativity, shocks,
//! or gameplay-engine integration. Preserve the legacy defaults.
use crate::physics::G_SI;
use crate::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState,GravityModel};
use crate::sph_kernel;
use crate::sph_rk_compatible::DensityChoice;
use crate::sph_surface_variation::{self,AnalyticEnergyGradient};
use crate::sph_rk_variational;

fn matched(cfg:AdaptiveConfig)->AdaptiveConfig {
    assert!(cfg.gravity_model==GravityModel::MatchedCubicSpline,
        "BG-PHY-23 requires the same Poisson-matched gravity kernel");
    cfg.with_newton_density()
}
fn zero_central()->Massive3 {
    Massive3{pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0}
}
fn cross(a:V3,b:V3)->V3 {
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn force_rms_ratio(a:&[V3],b:&[V3])->f64 {
    assert_eq!(a.len(),b.len());
    let mut error=0.0;
    let mut ref2=0.0;
    for (x,y) in a.iter().zip(b) {
        error+=(*x-*y).len_sq();
        ref2+=y.len_sq();
    }
    (error/ref2.max(1e-100)).sqrt()
}
/// Pair potential at the *already solved* smoothing lengths.
/// Avoid independent numerical smoothing or a Plummer substitution.
pub fn matched_gravity_potential_j(
    particles:&[Fluid3],state:&AdaptiveState,cfg:AdaptiveConfig,
)->f64 {
    matched(cfg);
    assert!(state.valid_for(particles.len()));
    let mut phi=0.0;
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let r=(particles[i].pos_m-particles[j].pos_m).len();
            let pair=0.5*(
                sph_kernel::potential(r,state.smoothing_m[i],cfg.kernel)+
                sph_kernel::potential(r,state.smoothing_m[j],cfg.kernel));
            phi+=G_SI*particles[i].mass_kg*particles[j].mass_kg*pair;
        }
    }
    phi
}
#[derive(Clone,Debug)]
pub struct CoupledSnapshot {
    pub total_force_n:Vec<V3>,
    pub pressure_force_n:Vec<V3>,
    pub gravity_force_n:Vec<V3>,
    pub cold_internal_energy_j:f64,
    pub gravity_potential_j:f64,
    pub total_potential_j:f64,
    pub pressure:AnalyticEnergyGradient,
    pub minimum_omega:f64,
    pub capped_fraction:f64,
    /// Net force and torque are computed before any projection or
    /// artificial recentring of particles.
    pub net_force_ratio:f64,
    pub net_torque_ratio:f64,
    /// Relative mismatch with existing fully matched SPH force when
    /// density_choice == Original; recorded for either density field.
    pub old_sph_relative_force_error:f64,
}
pub fn snapshot(
    particles:&[Fluid3],physics:Sph3Config,cfg:AdaptiveConfig,
    star_radius_m:f64,density_choice:DensityChoice,
)->CoupledSnapshot {
    let cfg=matched(cfg);
    assert!(particles.len()>=16&&physics.valid()&&star_radius_m>0.0);
    let pressure=sph_surface_variation::total_energy_gradient(
        particles,physics,cfg,star_radius_m,density_choice);
    let mut state=AdaptiveState::new(particles,star_radius_m,cfg);
    let conventional=sph_adaptive::forces(
        &zero_central(),particles,&mut state,physics,cfg);
    let gravity_potential_j=matched_gravity_potential_j(
        particles,&state,cfg);
    let n=particles.len();
    let mut total_force_n=Vec::with_capacity(n);
    let mut gravity_force_n=Vec::with_capacity(n);
    let mut old_force_n=Vec::with_capacity(n);
    let mut net=V3::default();
    let mut torque=V3::default();
    let mut fsum=0.0;
    let mut tsum=0.0;
    let mass=particles.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=particles.iter().fold(V3::default(),
        |s,x|s+x.pos_m*x.mass_kg)/mass;
    for (i,p) in particles.iter().enumerate() {
        let gravity=conventional.gravity_accelerations_m_s2[i]*p.mass_kg;
        let total=pressure.forces_n[i]+gravity;
        gravity_force_n.push(gravity);
        old_force_n.push(conventional.accelerations_m_s2[i]*p.mass_kg);
        total_force_n.push(total);
        net+=total;
        let arm=p.pos_m-center;
        torque+=cross(arm,total);
        fsum+=total.len();
        tsum+=arm.len()*total.len();
    }
    let old_sph_relative_force_error=force_rms_ratio(
        &total_force_n,&old_force_n);
    let capped=state.smoothing_m.iter().filter(|&&h|
        (h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8||
        (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8).count();
    let cold_internal_energy_j=pressure.potential_j;
    CoupledSnapshot{
        total_force_n,
        pressure_force_n:pressure.forces_n.clone(),
        gravity_force_n,
        cold_internal_energy_j,gravity_potential_j,
        total_potential_j:cold_internal_energy_j+gravity_potential_j,
        pressure,minimum_omega:state.omega.iter().copied()
            .fold(f64::INFINITY,f64::min),
        capped_fraction:capped as f64/n as f64,
        net_force_ratio:net.len()/fsum.max(1.0),
        net_torque_ratio:torque.len()/tsum.max(1.0),
        old_sph_relative_force_error,
    }
}
/// Independent scalar oracle: recompute exact h/rho, surface MLS,
/// Fermi u(rho_hat), and Poisson-matched pair potential at every
/// perturbed particle configuration. The force implementation above
/// is not called by this function.
pub fn independent_potential_j(
    particles:&[Fluid3],physics:Sph3Config,cfg:AdaptiveConfig,
    star_radius_m:f64,density_choice:DensityChoice,
)->f64 {
    let cfg=matched(cfg);
    assert!((16..=64).contains(&particles.len()));
    let internal=sph_rk_variational::internal_energy_j(
        particles,physics,cfg,star_radius_m,density_choice);
    let mut state=AdaptiveState::new(particles,star_radius_m,cfg);
    sph_adaptive::update_adaptive_densities(particles,&mut state,cfg);
    internal+matched_gravity_potential_j(particles,&state,cfg)
}
#[derive(Clone,Debug)]
pub struct CoupledOracle {
    pub force_n:Vec<V3>,
    pub potential_j:f64,
    pub step_m:f64,
    pub relative_force_difference:f64,
    pub relative_potential_difference:f64,
    pub net_force_ratio:f64,
    pub net_torque_ratio:f64,
}
pub fn independent_total_gradient(
    particles:&[Fluid3],physics:Sph3Config,cfg:AdaptiveConfig,
    star_radius_m:f64,density_choice:DensityChoice,
    step_fraction:f64,
)->CoupledOracle{
    assert!((1e-7..=1e-2).contains(&step_fraction));
    let epsilon=step_fraction*star_radius_m;
    let base=independent_potential_j(particles,physics,cfg,
        star_radius_m,density_choice);
    let analytic=snapshot(particles,physics,cfg,
        star_radius_m,density_choice);
    let mut updated=particles.to_vec();
    let mut force_n=vec![V3::default();particles.len()];
    let mut net=V3::default();
    let mut torque=V3::default();
    let mut scale=0.0;
    let mut torque_scale=0.0;
    let mass=particles.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=particles.iter().fold(V3::default(),
        |s,x|s+x.pos_m*x.mass_kg)/mass;
    for (i,p) in particles.iter().enumerate(){
        for axis in 0..3 {
            let original=p.pos_m;
            let delta=match axis {
                0=>V3::new(epsilon,0.0,0.0),
                1=>V3::new(0.0,epsilon,0.0),
                _=>V3::new(0.0,0.0,epsilon),
            };
            updated[i].pos_m=original+delta;
            let plus=independent_potential_j(
                &updated,physics,cfg,star_radius_m,density_choice);
            updated[i].pos_m=original-delta;
            let minus=independent_potential_j(
                &updated,physics,cfg,star_radius_m,density_choice);
            updated[i].pos_m=original;
            let value=-(plus-minus)/(2.0*epsilon);
            match axis{
                0=>force_n[i].x=value,
                1=>force_n[i].y=value,
                _=>force_n[i].z=value,
            };
        }
        net+=force_n[i];
        let arm=p.pos_m-center;
        torque+=cross(arm,force_n[i]);
        scale+=force_n[i].len();
        torque_scale+=arm.len()*force_n[i].len();
    }
    let relative_force_difference=force_rms_ratio(
        &analytic.total_force_n,&force_n);
    CoupledOracle{
        force_n,potential_j:base,step_m:epsilon,
        relative_force_difference,
        relative_potential_difference:
            (analytic.total_potential_j-base).abs()/base.abs().max(1.0),
        net_force_ratio:net.len()/scale.max(1.0),
        net_torque_ratio:torque.len()/torque_scale.max(1.0),
    }
}
#[derive(Clone,Debug)]
pub struct HamiltonianPilot{
    pub particles:Vec<Fluid3>,
    pub physics:Sph3Config,
    pub cfg:AdaptiveConfig,
    pub reference_radius_m:f64,
    pub density:DensityChoice,
}
impl HamiltonianPilot{
    pub fn gradient(&self)->CoupledSnapshot{
        snapshot(&self.particles,self.physics,self.cfg,
            self.reference_radius_m,self.density)
    }
    pub fn energy_j(&self)->f64{
        let potential=self.gradient().total_potential_j;
        potential+self.particles.iter()
            .map(|x|0.5*x.mass_kg*x.vel_m_s.len_sq()).sum::<f64>()
    }
    /// Fully reversible Hamiltonian kick-drift-kick integrator.
    /// No damping, fixed gravity softening, or gameplay callbacks.
    pub fn advance(&mut self,dt:f64){
        assert!(dt.is_finite()&&dt.abs()<=0.05);
        if dt==0.0{return;}
        let before=self.gradient();
        for (p,f) in self.particles.iter_mut().zip(before.total_force_n){
            p.vel_m_s+=f*(0.5*dt/p.mass_kg);
            p.pos_m+=p.vel_m_s*dt;
        }
        let after=self.gradient();
        for (p,f) in self.particles.iter_mut().zip(after.total_force_n){
            p.vel_m_s+=f*(0.5*dt/p.mass_kg);
        }
    }
    pub fn damped_advance(&mut self,dt:f64,tau_s:f64){
        assert!(tau_s.is_finite()&&tau_s>0.0);
        self.advance(dt);
        let factor=(-dt/tau_s).exp();
        let total=self.particles.iter().map(|x|x.mass_kg).sum::<f64>();
        let bulk=self.particles.iter().fold(V3::default(),
            |sum,x|sum+x.vel_m_s*x.mass_kg)/total;
        for x in &mut self.particles {
            x.vel_m_s=bulk+(x.vel_m_s-bulk)*factor;
        }
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn fixture(n:usize,kernel:sph_kernel::DensityKernel)
        ->(Vec<Fluid3>,Sph3Config,AdaptiveConfig,f64){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),n,17).0;
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==sph_kernel::DensityKernel::WendlandC4{
            cfg=cfg.with_wendland_c4();
        }
        cfg.eta=1.7;
        (p,Sph3Config::cold_white_dwarf(2.0).unwrap(),cfg,star.radius_m)
    }
    #[test]
    fn original_density_coupled_force_reproduces_existing_sph(){
        for kernel in [sph_kernel::DensityKernel::CubicM4,
                       sph_kernel::DensityKernel::WendlandC4] {
            let (p,physics,cfg,r)=fixture(16,kernel);
            let result=snapshot(&p,physics,cfg,r,DensityChoice::Original);
            assert!(result.old_sph_relative_force_error<5e-13,
                "legacy total force mismatch = {}",result.old_sph_relative_force_error);
            assert!(result.net_force_ratio<1e-11);
            assert!(result.net_torque_ratio<1e-11);
        }
    }
    #[test]
    fn complete_m4_and_c4_potentials_match_independent_fd_gradient(){
        for kernel in [sph_kernel::DensityKernel::CubicM4,
                       sph_kernel::DensityKernel::WendlandC4]{
            let (p,physics,cfg,r)=fixture(16,kernel);
            for density in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
                let result=snapshot(&p,physics,cfg,r,density);
                let oracle=independent_total_gradient(
                    &p,physics,cfg,r,density,1e-4);
                assert!(oracle.relative_force_difference<2e-5,
                    "full force mismatch {:?} {:?} {}",
                    kernel,density,oracle.relative_force_difference);
                assert!(oracle.relative_potential_difference<1e-12);
                assert!(result.net_force_ratio<1e-10);
                assert!(result.net_torque_ratio<1e-10);
            }
        }
    }
    #[test]
    fn coupled_hamiltonian_energy_and_reversal_are_bounded(){
        let (p,physics,cfg,r)=fixture(16,sph_kernel::DensityKernel::WendlandC4);
        let mut pilot=HamiltonianPilot{
            particles:p,physics,cfg,reference_radius_m:r,
            density:DensityChoice::SurfaceBlended,
        };
        let before=pilot.clone();
        let e0=pilot.energy_j();
        for _ in 0..3 {pilot.advance(0.002);}
        let drift=(pilot.energy_j()-e0).abs()/e0.abs().max(1.0);
        assert!(drift<1e-5,"coupled Hamiltonian drift={drift}");
        for _ in 0..3 {pilot.advance(-0.002);}
        let dx=pilot.particles.iter().zip(&before.particles)
            .map(|(a,b)|(a.pos_m-b.pos_m).len()).fold(0.0_f64,f64::max);
        assert!(dx/r<1e-9,"reverse position mismatch {}",dx/r);
    }
}
