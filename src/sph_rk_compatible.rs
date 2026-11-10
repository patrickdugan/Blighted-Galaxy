//! BG-PHY-21: isolated pressure-only compatible-work CRK pilot.
//!
//! The opt-in pilot uses the full fixed-volume QUERY gradients from
//! sph_rk.rs, antisymmetrized pairwise as in Frontiere et al. (2017).
//! The pressure force for an unordered i<j pair is
//!
//! F_ij = -V_i V_j (P_i+P_j) * 0.5*(grad W^R_ij-grad W^R_ji).
//!
//! V_i=m_i/rho_SPH,i are the SAME volumes used to construct RK moments.
//! The EOS pressure may optionally use the BG-PHY-20 surface-gated
//! frozen density, making this a useful pilot of the measured correction.
//!
//! For a kick of duration dt, v_mid=(v_before+v_after)/2.
//! Pair work is (v_mid,i-v_mid,j) dot (F_ij dt).
//! Assigning minus half this work to each particle's internal-energy
//! reservoir conserves total K+U to rounding precision. Pair forces
//! conserve linear momentum. Noncentral RK forces conserve angular
//! momentum only approximately; the benchmark measures the drift.
//!
//! Crucial boundary: evolving compatible-work U while evaluating P
//! using a separately estimated rho is NOT a full variational cold
//! barotrope. The missing density/position/h derivatives and the
//! evolving U vs u_EOS(rho) discrepancy are measured separately.
//! Gravity, viscosity, and gameplay engine coupling are excluded.
use crate::sph3::{self,Fluid3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use crate::sph_density_audit::{self,DensityEstimator};
use crate::sph_rk::{RkGeometry,RkIdentityAudit};

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum DensityChoice {
    Original,
    SurfaceBlended,
}
#[derive(Clone,Copy,Debug)]
pub struct RkPairForce {
    pub i:usize,
    pub j:usize,
    /// Total force applied to i. Exactly minus the force on j.
    pub on_i_newton:V3,
}
#[derive(Clone,Debug)]
pub struct RkPressureSnapshot {
    pub pairs:Vec<RkPairForce>,
    pub accelerations_m_s2:Vec<V3>,
    pub eos_u_j_per_kg:Vec<f64>,
    pub valid_rk_fraction:f64,
    pub rk: RkIdentityAudit,
    pub surface_density_log_displacement_rms:f64,
}
#[derive(Clone,Debug)]
pub struct RkPilotState {
    pub particles:Vec<Fluid3>,
    pub u_j_per_kg:Vec<f64>,
    pub sph_state:AdaptiveState,
    pub physics:Sph3Config,
    pub cfg:AdaptiveConfig,
    pub density:DensityChoice,
}
fn cross(a:V3,b:V3)->V3{
    V3::new(a.y*b.z-a.z*b.y,
        a.z*b.x-a.x*b.z,
        a.x*b.y-a.y*b.x)
}
pub fn pressure_snapshot(
    particles:&[Fluid3],state:&mut AdaptiveState,
    physics:Sph3Config,cfg:AdaptiveConfig,density:DensityChoice,
)->RkPressureSnapshot{
    assert!(particles.len()>=8&&state.valid_for(particles.len()));
    sph_adaptive::update_adaptive_densities(particles,state,cfg);
    let reconstructed=sph_density_audit::estimate(
        particles,state,cfg.kernel);
    let field=reconstructed.fields.iter().find(|x|
        x.method==match density{
            DensityChoice::Original=>DensityEstimator::Original,
            DensityChoice::SurfaceBlended=>DensityEstimator::SurfaceBlendedMls,
        }).expect("selected density candidate must exist");
    let rk=RkGeometry::from_state(particles,state,cfg.kernel);
    let identities=rk.identities(particles);
    let n=particles.len();
    let pressures=field.rho.iter().map(|rho|
        sph3::pressure_pa(*rho,physics)).collect::<Vec<_>>();
    let energies=field.rho.iter().map(|rho|
        sph3::internal_energy_j_per_kg(*rho,physics)).collect::<Vec<_>>();
    let mut log_displacement=0.0;
    for (rho,baseline) in field.rho.iter().zip(&state.densities_kg_m3){
        log_displacement+=(*rho/baseline).ln().powi(2);
    }
    let mut accelerations=vec![V3::default();n];
    let mut pairs=Vec::with_capacity(n*(n-1)/2);
    for i in 0..n {
        for j in i+1..n {
            let grad=rk.antisymmetric_pair_gradient(particles,i,j);
            // rho_SPH volumes appear in the RK moment matrix as well.
            let force=grad*
                (-rk.volumes[i]*rk.volumes[j]*(pressures[i]+pressures[j]));
            if force.len_sq()==0.0{continue;}
            accelerations[i]+=force*(1.0/particles[i].mass_kg);
            accelerations[j]-=force*(1.0/particles[j].mass_kg);
            pairs.push(RkPairForce{i,j,on_i_newton:force});
        }
    }
    RkPressureSnapshot {
        pairs,accelerations_m_s2:accelerations,eos_u_j_per_kg:energies,
        valid_rk_fraction:identities.valid_fraction,rk:identities,
        surface_density_log_displacement_rms:(log_displacement/n as f64).sqrt(),
    }
}
/// Exact compatible work for an impulse over dt with frozen pair forces.
/// Velocities all advance simultaneously, and work uses their mid-impulse
/// values, so the sum_i m_i (delta K_i + delta u_i) is zero.
pub fn compatible_kick(
    particles:&mut [Fluid3],u:&mut [f64],
    force:&RkPressureSnapshot,dt:f64,
){
    assert!(dt.is_finite()&&particles.len()==u.len()&&
        particles.len()==force.accelerations_m_s2.len());
    let mut impulse=vec![V3::default();particles.len()];
    for p in &force.pairs {
        let delta=p.on_i_newton*dt;
        impulse[p.i]+=delta;
        impulse[p.j]-=delta;
    }
    let mid=particles.iter().zip(&impulse)
        .map(|(p,dp)|p.vel_m_s+*dp*(0.5/p.mass_kg))
        .collect::<Vec<_>>();
    for pair in &force.pairs {
        let work=(mid[pair.i]-mid[pair.j]).dot(pair.on_i_newton*dt);
        u[pair.i]-=0.5*work/particles[pair.i].mass_kg;
        u[pair.j]-=0.5*work/particles[pair.j].mass_kg;
    }
    for (p,dp) in particles.iter_mut().zip(impulse){
        p.vel_m_s+=dp*(1.0/p.mass_kg);
    }
}
impl RkPilotState {
    pub fn new(
        particles:Vec<Fluid3>,physics:Sph3Config,cfg:AdaptiveConfig,
        reference_radius_m:f64,density:DensityChoice,
    )->Self {
        let mut sph_state=AdaptiveState::new(
            &particles,reference_radius_m,cfg);
        let force=pressure_snapshot(
            &particles,&mut sph_state,physics,cfg,density);
        Self {
            particles,u_j_per_kg:force.eos_u_j_per_kg,
            sph_state,physics,cfg,density,
        }
    }
    pub fn forces(&mut self)->RkPressureSnapshot {
        pressure_snapshot(&self.particles,&mut self.sph_state,
            self.physics,self.cfg,self.density)
    }
    /// Velocity Verlet; +/- dt permit time-reversal tests.
    /// Gravity is intentionally absent from this pilot.
    pub fn advance(&mut self,dt:f64){
        assert!(dt.is_finite()&&dt.abs()<=0.05);
        if dt==0.0{return;}
        let first=self.forces();
        compatible_kick(
            &mut self.particles,&mut self.u_j_per_kg,&first,0.5*dt);
        for p in &mut self.particles{
            p.pos_m+=p.vel_m_s*dt;
        }
        let last=self.forces();
        compatible_kick(
            &mut self.particles,&mut self.u_j_per_kg,&last,0.5*dt);
    }
    pub fn kinetic_j(&self)->f64{
        self.particles.iter().map(|p|0.5*p.mass_kg*p.vel_m_s.len_sq()).sum()
    }
    pub fn internal_j(&self)->f64{
        self.particles.iter().zip(&self.u_j_per_kg)
            .map(|(p,u)|p.mass_kg*u).sum()
    }
    pub fn total_energy_j(&self)->f64{
        self.kinetic_j()+self.internal_j()
    }
    pub fn momentum(&self)->V3 {
        self.particles.iter().fold(V3::default(),
            |sum,p|sum+p.vel_m_s*p.mass_kg)
    }
    pub fn angular_momentum(&self)->V3 {
        let mass=self.particles.iter().map(|p|p.mass_kg).sum::<f64>();
        let center=self.particles.iter().fold(V3::default(),
            |sum,p|sum+p.pos_m*p.mass_kg)/mass;
        let frame_v=self.momentum()*(1.0/mass);
        self.particles.iter().fold(V3::default(),|sum,p|
            sum+cross(p.pos_m-center,
                (p.vel_m_s-frame_v)*p.mass_kg))
    }
    pub fn eos_consistency_rms(&mut self)->f64{
        let snapshot=self.forces();
        let mut numerator=0.0;
        let mut denominator=0.0;
        for (i,p) in self.particles.iter().enumerate() {
            numerator+=p.mass_kg*
                (self.u_j_per_kg[i]-snapshot.eos_u_j_per_kg[i]).powi(2);
            denominator+=p.mass_kg*snapshot.eos_u_j_per_kg[i].powi(2);
        }
        (numerator/denominator.max(1.0)).sqrt()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{profile_for_mass,DegenerateElectrons};
    fn setup(n:usize,choice:DensityChoice)->RkPilotState{
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,V3::new(1e8,0.0,0.0),
            V3::default(),n,17).0;
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        cfg.eta=1.7;
        RkPilotState::new(p,Sph3Config::cold_white_dwarf(2.0).unwrap(),
            cfg,star.radius_m,choice)
    }
    #[test]
    fn reproducing_force_pairs_are_momentum_conserving(){
        for density in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
            let mut state=setup(32,density);
            let force=state.forces();
            let net=state.particles.iter().zip(&force.accelerations_m_s2)
                .fold(V3::default(),|sum,(p,a)|sum+*a*p.mass_kg);
            let scale=state.particles.iter().zip(&force.accelerations_m_s2)
                .map(|(p,a)|p.mass_kg*a.len()).sum::<f64>().max(1.0);
            assert!(net.len()/scale<1e-12);
            assert!(force.valid_rk_fraction>=0.5);
            assert!(force.rk.maximum_constant_error<1e-8);
            assert!(force.rk.maximum_linear_error<1e-8);
            assert!(force.rk.maximum_query_gradient_error<1e-8);
        }
    }
    #[test]
    fn compatible_work_conserves_total_energy_and_reverses(){
        let mut state=setup(32,DensityChoice::SurfaceBlended);
        let initial=state.clone();
        let energy=state.total_energy_j();
        let momentum=state.momentum();
        let dt=0.001;
        for _ in 0..8{state.advance(dt);}
        let drift=(state.total_energy_j()-energy).abs()/energy.abs();
        let mass=state.particles.iter().map(|p|p.mass_kg).sum::<f64>();
        let dp=(state.momentum()-momentum).len()/(mass*1e5);
        assert!(drift<1e-11,"compatible energy drift {drift:e}");
        assert!(dp<1e-12,"momentum drift {dp:e}");
        for _ in 0..8{state.advance(-dt);}
        for (now,original) in state.particles.iter().zip(&initial.particles){
            assert!((now.pos_m-original.pos_m).len()/8.4e6<1e-9);
            assert!((now.vel_m_s-original.vel_m_s).len()/1e6<1e-9);
        }
    }
    #[test]
    fn uniform_boost_does_not_change_acceleration(){
        let mut a=setup(32,DensityChoice::Original);
        let before=a.forces();
        for p in &mut a.particles{
            p.vel_m_s+=V3::new(3000.0,-7000.0,1600.0);
        }
        let after=a.forces();
        for (x,y) in before.accelerations_m_s2.iter().zip(after.accelerations_m_s2){
            assert_eq!(*x,y);
        }
    }
}
