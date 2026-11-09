//! Damped isolated 3D white-dwarf SPH initialization and relaxation.
 //!
 //! A hydrostatic 1D density law is sampled with equal-mass antipodal
 //! particles; the 3D discretization has finite particle noise. The
 //! relaxation evolves self-gravity and cold-electron pressure, with a
 //! temporary artificial velocity damping term, then reports the net
 //! SPH residual and departure from the reference density profile.
 //!
 //! Damping is intentionally dissipative: energy conservation is NOT
 //! expected during relaxation. The conservative post-relaxation solver
 //! remains in sph3.rs. This is a convergence experiment, not proof that
 //! 64–192 particles represent a realistic equilibrium white dwarf.
use crate::sph3::{
    self,Fluid3,Massive3,Sph3Config,Sph3Eos,V3,
};
use crate::white_dwarf::WhiteDwarfProfile;
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};

#[derive(Debug,Clone,Copy)]
pub struct RelaxParameters {
    pub iterations:usize,
    pub dt_s:f64,
    pub damping_time_s:f64,
}
impl RelaxParameters {
    pub fn baseline()->Self {
        Self{iterations:750,dt_s:0.001,
            damping_time_s:0.15}
    }
    pub fn valid(&self)->bool {
        self.iterations>0&&self.iterations<=50_000&&
        self.dt_s.is_finite()&&self.dt_s>0.0&&self.dt_s<=0.1&&
        self.damping_time_s.is_finite()&&self.damping_time_s>0.0
    }
}
#[derive(Debug,Clone,Copy)]
pub struct RelaxReport {
    pub iterations:usize,
    pub initial_accel_rms_m_s2:f64,
    pub final_accel_rms_m_s2:f64,
    pub initial_density_log_rms:f64,
    pub final_density_log_rms:f64,
    pub final_max_speed_relative_m_s:f64,
    pub final_rms_speed_relative_m_s:f64,
    pub max_radius_m:f64,
    pub center_shift_m:f64,
    pub mass_relative_error:f64,
}
impl RelaxReport {
    pub fn force_reduction(&self)->f64 {
        self.final_accel_rms_m_s2/self.initial_accel_rms_m_s2
    }
}
fn barycenter(p:&[Fluid3])->(f64,V3,V3) {
    let mut m=0.0;let mut pos=V3::default();let mut momentum=V3::default();
    for x in p {
        m+=x.mass_kg;
        pos+=x.pos_m*x.mass_kg;
        momentum+=x.vel_m_s*x.mass_kg;
    }
    (m,pos/m,momentum/m)
}
fn rms_acceleration(p:&[Fluid3],cfg:Sph3Config)->f64 {
    let center=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let f=sph3::accelerations(&center,p,cfg);
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    (f.parcel_accel_m_s2.iter().zip(p).map(|(a,x)|
        x.mass_kg*a.len_sq()).sum::<f64>()/m).sqrt()
}
fn log_density_error(
    p:&[Fluid3],cfg:Sph3Config,profile:&WhiteDwarfProfile,cm:V3,
)->f64 {
    let rho=sph3::densities_kg_m3(p,cfg.smoothing_m);
    let e=rho.iter().zip(p).map(|(density,x)|{
        let r=(x.pos_m-cm).len();
        // Avoid logarithmic divergence at the density-zero surface.
        let target=profile.density_at(r)
            .max(profile.rho_c_kg_m3*1e-5);
        (density/target).ln().powi(2)
    }).sum::<f64>();
    (e/p.len() as f64).sqrt()
}
pub fn relax_isolated_white_dwarf(
    p:&mut [Fluid3],cfg:Sph3Config,profile:&WhiteDwarfProfile,
    options:RelaxParameters,
)->Result<RelaxReport,&'static str> {
    if !options.valid()||p.len()<4||!cfg.valid()||
        !matches!(cfg.eos,Sph3Eos::ColdDegenerate{..}) {
        return Err("invalid relaxation model or SPH particles");
    }
    let (initial_mass,original_center,bulk_velocity)=barycenter(p);
    if !initial_mass.is_finite()||
        (initial_mass-profile.mass_kg).abs()/profile.mass_kg>1e-5 {
        return Err("SPH mass does not match hydrostatic profile");
    }
    let initial_force=rms_acceleration(p,cfg);
    let initial_density=log_density_error(p,cfg,profile,original_center);
    let decay=(-options.dt_s/options.damping_time_s).exp();
    let mut central=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    for _ in 0..options.iterations {
        // The central body is massless in this isolated stellar-relaxation
        // phase. Its formal acceleration has no effect on the parcels.
        sph3::advance(&mut central,p,cfg,options.dt_s);
        let (_,current_center,current_bulk)=barycenter(p);
        for x in p.iter_mut() {
            x.pos_m+=original_center-current_center;
            x.vel_m_s=bulk_velocity+
                (x.vel_m_s-current_bulk)*decay;
        }
        if p.iter().any(|x|
            !x.pos_m.finite()||!x.vel_m_s.finite()) {
            return Err("relaxation produced nonfinite particles");
        }
    }
    let (final_mass,center,bulk)=barycenter(p);
    let final_force=rms_acceleration(p,cfg);
    let final_density=log_density_error(p,cfg,profile,center);
    let max_v=p.iter().map(|x|(x.vel_m_s-bulk).len())
        .fold(0.0_f64,f64::max);
    let rms_v=(p.iter().map(|x|
        x.mass_kg*(x.vel_m_s-bulk).len_sq()).sum::<f64>()/
        final_mass).sqrt();
    let max_radius=p.iter().map(|x|(x.pos_m-center).len())
        .fold(0.0_f64,f64::max);
    Ok(RelaxReport {
        iterations:options.iterations,
        initial_accel_rms_m_s2:initial_force,
        final_accel_rms_m_s2:final_force,
        initial_density_log_rms:initial_density,
        final_density_log_rms:final_density,
        final_max_speed_relative_m_s:max_v,
        final_rms_speed_relative_m_s:rms_v,
        max_radius_m:max_radius,
        center_shift_m:(center-original_center).len(),
        mass_relative_error:(final_mass-initial_mass).abs()/initial_mass,
    })
}
#[derive(Debug,Clone,Copy)]
pub struct AdaptiveRelaxReport {
    pub report:RelaxReport,
    pub minimum_omega:f64,
    pub minimum_h_m:f64,
    pub maximum_h_m:f64,
    /// Number of particles whose h is clipped at a configured bound.
    pub capped_fraction:f64,
    /// Max-radius mismatch with the reference hydrostatic star.
    pub radius_bias:f64,
}
fn adaptive_log_density(
    p:&[Fluid3],state:&AdaptiveState,profile:&WhiteDwarfProfile,cm:V3,
)->f64 {
    let e=p.iter().zip(&state.densities_kg_m3).map(|(x,rho)|{
        let radius=(x.pos_m-cm).len();
        let target=profile.density_at(radius)
            .max(profile.rho_c_kg_m3*1e-5);
        (rho/target).ln().powi(2)
    }).sum::<f64>();
    (e/p.len() as f64).sqrt()
}
/// Experimental adaptive-h relaxation. The pairwise pressure force includes
/// Ω_i grad-h terms, while central and self-gravity use fixed Plummer
/// softening. Damping removes kinetic energy and is disabled during
/// subsequent undamped hold tests.
pub fn relax_isolated_adaptive(
    p:&mut [Fluid3],physics:Sph3Config,profile:&WhiteDwarfProfile,
    config:AdaptiveConfig,state:&mut AdaptiveState,
    options:RelaxParameters,
)->Result<AdaptiveRelaxReport,&'static str> {
    if !options.valid()||p.len()<4||!physics.valid()||
        !config.valid()||!state.valid_for(p.len())||
        !matches!(physics.eos,Sph3Eos::ColdDegenerate{..}) {
        return Err("invalid adaptive relaxation inputs");
    }
    let (initial_mass,initial_center,initial_bulk)=barycenter(p);
    if (initial_mass-profile.mass_kg).abs()/profile.mass_kg>1e-5 {
        return Err("incorrect stellar particle mass");
    }
    let first_force=sph_adaptive::isolated_force_rms_m_s2(
        p,state,physics,config);
    let initial_density=adaptive_log_density(p,state,profile,initial_center);
    let decay=(-options.dt_s/options.damping_time_s).exp();
    let mut central=Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    for _ in 0..options.iterations {
        sph_adaptive::advance(&mut central,p,state,physics,config,
            options.dt_s);
        let (_,center,bulk)=barycenter(p);
        for x in p.iter_mut() {
            x.pos_m+=initial_center-center;
            x.vel_m_s=initial_bulk+(x.vel_m_s-bulk)*decay;
        }
        if p.iter().any(|x|
            !x.pos_m.finite()||!x.vel_m_s.finite()) {
            return Err("adaptive relaxation diverged");
        }
    }
    let (final_mass,center,bulk)=barycenter(p);
    let final_force=sph_adaptive::isolated_force_rms_m_s2(
        p,state,physics,config);
    let final_density=adaptive_log_density(p,state,profile,center);
    let max_speed=p.iter().map(|x|(x.vel_m_s-bulk).len())
        .fold(0.0_f64,f64::max);
    let rms_speed=(p.iter().map(|x|
        x.mass_kg*(x.vel_m_s-bulk).len_sq()).sum::<f64>()/
        final_mass).sqrt();
    let max_radius=p.iter().map(|x|(x.pos_m-center).len())
        .fold(0.0_f64,f64::max);
    let capped=state.smoothing_m.iter().filter(|h|
        (**h-config.min_h_m).abs()<config.min_h_m*1e-8||
        (**h-config.max_h_m).abs()<config.max_h_m*1e-8).count();
    Ok(AdaptiveRelaxReport {
        report:RelaxReport{
            iterations:options.iterations,
            initial_accel_rms_m_s2:first_force,
            final_accel_rms_m_s2:final_force,
            initial_density_log_rms:initial_density,
            final_density_log_rms:final_density,
            final_max_speed_relative_m_s:max_speed,
            final_rms_speed_relative_m_s:rms_speed,
            max_radius_m:max_radius,
            center_shift_m:(center-initial_center).len(),
            mass_relative_error:(final_mass-initial_mass).abs()/initial_mass,
        },
        minimum_omega:state.omega.iter().copied()
            .fold(f64::INFINITY,f64::min),
        minimum_h_m:state.smoothing_m.iter().copied()
            .fold(f64::INFINITY,f64::min),
        maximum_h_m:state.smoothing_m.iter().copied()
            .fold(0.0_f64,f64::max),
        capped_fraction:capped as f64/state.smoothing_m.len() as f64,
        radius_bias:max_radius/profile.radius_m-1.0,
    })
}

/// Remove damping and track the radius, force, velocity and energy drift
/// under the same adaptive force law. This does *not* make a coarse
/// particle star automatically equivalent to the continuum hydrostatic
/// reference: stability gates must be checked independently.
#[derive(Debug,Clone,Copy)]
pub struct UndampedHoldReport {
    pub elapsed_s:f64,
    pub mass_relative_drift:f64,
    pub momentum_relative_drift:f64,
    pub energy_relative_drift:f64,
    pub final_radius_m:f64,
    pub radius_relative_change:f64,
}
pub fn adaptive_undamped_hold(
    p:&mut [Fluid3],physics:Sph3Config,
    cfg:AdaptiveConfig,state:&mut AdaptiveState,
    dt_s:f64,steps:usize,
)->Result<UndampedHoldReport,&'static str> {
    if steps==0||steps>100_000||!dt_s.is_finite()||
        dt_s<=0.0||dt_s>0.1||!cfg.valid()||
        !state.valid_for(p.len()) {return Err("invalid hold parameters");}
    let mut central=Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let (mass0,center0,_)=barycenter(p);
    let momentum0=sph3::total_momentum(&central,p);
    let energy0=sph_adaptive::adaptive_energy_j(
        &central,p,state,physics,cfg);
    let radius0=p.iter().map(|x|(x.pos_m-center0).len())
        .fold(0.0_f64,f64::max);
    for _ in 0..steps {
        sph_adaptive::advance(&mut central,p,state,physics,cfg,dt_s);
    }
    let (mass,center,_)=barycenter(p);
    let momentum=sph3::total_momentum(&central,p);
    let energy=sph_adaptive::adaptive_energy_j(
        &central,p,state,physics,cfg);
    let radius=p.iter().map(|x|(x.pos_m-center).len())
        .fold(0.0_f64,f64::max);
    Ok(UndampedHoldReport{
        elapsed_s:dt_s*steps as f64,
        mass_relative_drift:(mass-mass0).abs()/mass0,
        momentum_relative_drift:(momentum-momentum0).len()/
            (mass0*3e5),
        energy_relative_drift:(energy-energy0).abs()/
            energy0.abs().max(1.0),
        final_radius_m:radius,
        radius_relative_change:radius/radius0-1.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{profile_for_mass,DegenerateElectrons};
    fn setup(n:usize)->(WhiteDwarfProfile,Vec<Fluid3>,Sph3Config) {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let p=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let samples=sph3::seed_from_hydrostatic(
            &p,V3::new(1e8,-4e7,0.0),
            V3::new(2e4,-4e4,0.0),n);
        (p,samples,Sph3Config::cold_white_dwarf(2.0).unwrap())
    }
    #[test]
    fn damped_relaxation_preserves_mass_center_and_bulk_velocity() {
        let (profile,mut particles,cfg)=setup(64);
        let (m,center,bulk)=barycenter(&particles);
        let report=relax_isolated_white_dwarf(
            &mut particles,cfg,&profile,
            RelaxParameters{iterations:240,dt_s:0.001,
                damping_time_s:0.13}).unwrap();
        let (after_m,after_center,after_bulk)=barycenter(&particles);
        assert!((after_m-m).abs()/m<1e-13);
        assert!((after_center-center).len()<1e-5);
        assert!((after_bulk-bulk).len()<1e-6);
        assert!(report.final_accel_rms_m_s2.is_finite());
        assert!(report.final_density_log_rms.is_finite());
        assert!(report.max_radius_m<profile.radius_m*4.0);
    }
    #[test]
    fn adaptive_relaxation_preserves_mass_and_center() {
        let (profile,mut p,physics)=setup(48);
        let cfg=AdaptiveConfig::from_star_radius(profile.radius_m);
        let mut state=AdaptiveState::new(&p,profile.radius_m,cfg);
        let report=relax_isolated_adaptive(&mut p,physics,&profile,
            cfg,&mut state,RelaxParameters{
                iterations:80,dt_s:0.002,damping_time_s:0.14,
            }).unwrap();
        assert!(report.report.mass_relative_error<1e-12);
        assert!(report.report.center_shift_m<1e-5);
        assert!(report.report.final_density_log_rms.is_finite());
        assert!(report.minimum_omega>=0.02);
    }
    #[test]
    fn adaptive_hold_after_damping_preserves_bulk_momentum() {
        let (profile,mut p,physics)=setup(48);
        let cfg=AdaptiveConfig::from_star_radius(profile.radius_m);
        let mut state=AdaptiveState::new(&p,profile.radius_m,cfg);
        relax_isolated_adaptive(&mut p,physics,&profile,cfg,&mut state,
            RelaxParameters{
                iterations:30,dt_s:0.002,damping_time_s:0.14,
            }).unwrap();
        let report=adaptive_undamped_hold(&mut p,physics,cfg,&mut state,
            0.002,30).unwrap();
        assert!(report.mass_relative_drift<1e-13);
        assert!(report.momentum_relative_drift<1e-10);
        assert!(report.final_radius_m.is_finite());
        assert!(report.energy_relative_drift.is_finite());
    }
    #[test]
    fn zero_damping_time_and_wrong_eos_are_rejected() {
        let (profile,mut particles,cfg)=setup(32);
        let bad=RelaxParameters{iterations:10,dt_s:0.01,
            damping_time_s:0.0};
        assert!(relax_isolated_white_dwarf(
            &mut particles,cfg,&profile,bad).is_err());
        assert!(relax_isolated_white_dwarf(
            &mut particles,Sph3Config::white_dwarf_proxy(),
            &profile,RelaxParameters::baseline()).is_err());
    }
    #[test]
    fn isolated_pair_accelerations_have_zero_bulk_force() {
        let (profile,particles,cfg)=setup(48);
        let r0=rms_acceleration(&particles,cfg);
        assert!(r0>0.0);
        assert!(profile.mass_kg>0.0);
        let central=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let acc=sph3::accelerations(&central,&particles,cfg);
        let net=particles.iter().zip(acc.parcel_accel_m_s2.iter())
            .fold(V3::default(),|total,(x,a)|
                total+*a*x.mass_kg);
        let typical=r0*profile.mass_kg;
        assert!(net.len()/typical<1e-11);
    }
}
