//! Radial pressure/gravity budget of an isolated cold white-dwarf SPH star.
 //!
 //! Uses the same pairwise forces as the live adaptive-h solver, including
 //! the gravity and pressure grad-h terms, and groups their *radial*
 //! components by enclosed particle mass. A hydrostatic 1D comparison
 //! supplies g_ref=G*M_enclosed(r)/r² independently from the SPH solver.
 //!
 //! A shell's mean net force can cancel even if individuals have large
 //! residuals, so a mass-weighted RMS residual is reported too.
use crate::physics::G_SI;
use crate::sph3::{Fluid3,Massive3,Sph3Config,V3};
use crate::sph_adaptive::{
    self,AdaptiveConfig,AdaptiveState,
};
use crate::white_dwarf::WhiteDwarfProfile;

#[derive(Debug,Clone,Copy)]
pub struct RadialShellBalance {
    pub mass_fraction_lo:f64,
    pub mass_fraction_hi:f64,
    pub particle_count:usize,
    pub mean_radius_m:f64,
    /// Mean radially outward pressure acceleration, m/s².
    pub pressure_out_m_s2:f64,
    /// Mean radially inward gravitational acceleration, m/s².
    pub gravity_in_m_s2:f64,
    pub hydrostatic_gravity_in_m_s2:f64,
    /// Signed radial acceleration; zero at radial balance.
    pub net_radial_m_s2:f64,
    /// RMS of radial net acceleration within the shell.
    pub net_radial_rms_m_s2:f64,
    pub tangential_net_rms_m_s2:f64,
}
#[derive(Debug,Clone)]
pub struct StellarBalance {
    pub shells:Vec<RadialShellBalance>,
    pub center_of_mass_m:V3,
    pub pressure_virial_j:f64,
    pub gravitational_virial_j:f64,
    /// Pressure virial divided by the magnitude of gravity virial.
    pub pressure_to_gravity_virial:f64,
    pub net_force_rms_m_s2:f64,
    pub gravity_force_rms_m_s2:f64,
    pub residual_to_gravity_rms:f64,
    pub maximum_radial_shell_residual:f64,
}
pub fn radial_balance(
    particles:&[Fluid3],state:&mut AdaptiveState,physics:Sph3Config,
    cfg:AdaptiveConfig,reference:&WhiteDwarfProfile,n_shells:usize,
)->StellarBalance {
    assert!(n_shells>0&&n_shells<=particles.len()&&
        state.valid_for(particles.len()));
    let mass=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    assert!(mass>0.0&&mass.is_finite());
    let center=particles.iter().fold(V3::default(),|c,p|
        c+p.pos_m*p.mass_kg)/mass;
    let empty=Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let f=sph_adaptive::forces(&empty,particles,state,physics,cfg);
    let mut order=particles.iter().enumerate()
        .map(|(i,p)|((p.pos_m-center).len(),i))
        .collect::<Vec<_>>();
    order.sort_by(|a,b|a.0.total_cmp(&b.0));
    let mut bins=vec![Vec::<usize>::new();n_shells];
    // Divide by approximately equal enclosed *particle mass*, rather than
    // by radial distance (the source 1D density varies strongly with r).
    let mut cumulative=0.0;
    for (_,i) in &order {
        let center_fraction=(cumulative+0.5*particles[*i].mass_kg)/mass;
        let bin=((center_fraction*n_shells as f64).floor() as usize)
            .min(n_shells-1);
        bins[bin].push(*i);
        cumulative+=particles[*i].mass_kg;
    }
    let mut pressure_virial=0.0;
    let mut gravity_virial=0.0;
    let mut net_square=0.0;
    let mut gravity_square=0.0;
    let mut shells=Vec::with_capacity(n_shells);
    let mut max_shell_residual=0.0_f64;
    let mut running_fraction=0.0;
    for bin in bins {
        let shell_mass=bin.iter().map(|&i|particles[i].mass_kg)
            .sum::<f64>();
        let mut radius=0.0;
        let mut pressure_radial=0.0;
        let mut gravity_radial=0.0;
        let mut reference_gravity=0.0;
        let mut radial_residual_sq=0.0;
        let mut tangential_sq=0.0;
        for &i in &bin {
            let particle=&particles[i];
            let displacement=particle.pos_m-center;
            let r=displacement.len();
            let unit=displacement.normalized();
            let a_pressure=f.pressure_accelerations_m_s2[i];
            let a_gravity=f.gravity_accelerations_m_s2[i];
            let a_net=f.accelerations_m_s2[i];
            let ap=a_pressure.dot(unit);
            let ag=-a_gravity.dot(unit);
            let an=a_net.dot(unit);
            let tangent=a_net-unit*an;
            let reference_g=if r>0.0 {
                G_SI*reference.enclosed_mass_at(r)/(r*r)
            } else {0.0};
            let m=particle.mass_kg;
            radius+=m*r;
            pressure_radial+=m*ap;
            gravity_radial+=m*ag;
            reference_gravity+=m*reference_g;
            radial_residual_sq+=m*an*an;
            tangential_sq+=m*tangent.len_sq();
            pressure_virial+=m*displacement.dot(a_pressure);
            gravity_virial+=m*displacement.dot(a_gravity);
            net_square+=m*a_net.len_sq();
            gravity_square+=m*a_gravity.len_sq();
        }
        let inv=1.0/shell_mass;
        let ap=pressure_radial*inv;
        let ag=gravity_radial*inv;
        let net=ap-ag;
        let normal=ap.abs()+ag.abs();
        let relative=if normal>0.0 {net.abs()/normal}else{0.0};
        max_shell_residual=max_shell_residual.max(relative);
        let upper=running_fraction+shell_mass/mass;
        shells.push(RadialShellBalance{
            mass_fraction_lo:running_fraction,
            mass_fraction_hi:upper,
            particle_count:bin.len(),
            mean_radius_m:radius*inv,
            pressure_out_m_s2:ap,
            gravity_in_m_s2:ag,
            hydrostatic_gravity_in_m_s2:reference_gravity*inv,
            net_radial_m_s2:net,
            net_radial_rms_m_s2:(radial_residual_sq*inv).sqrt(),
            tangential_net_rms_m_s2:(tangential_sq*inv).sqrt(),
        });
        running_fraction=upper;
    }
    let net_rms=(net_square/mass).sqrt();
    let grav_rms=(gravity_square/mass).sqrt();
    let virial_ratio=if gravity_virial<0.0 {
        pressure_virial/(-gravity_virial)
    } else {f64::NAN};
    StellarBalance{
        shells,center_of_mass_m:center,
        pressure_virial_j:pressure_virial,
        gravitational_virial_j:gravity_virial,
        pressure_to_gravity_virial:virial_ratio,
        net_force_rms_m_s2:net_rms,
        gravity_force_rms_m_s2:grav_rms,
        residual_to_gravity_rms:net_rms/grav_rms.max(1e-30),
        maximum_radial_shell_residual:max_shell_residual,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::{self,StarSampling};
    use crate::white_dwarf::{profile_for_mass,DegenerateElectrons};
    fn fixture(n:usize)->(WhiteDwarfProfile,Vec<Fluid3>,Sph3Config,AdaptiveConfig) {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let particles=sph_sampling::seeded_star(&star,
            V3::new(1e8,0.0,0.0),V3::default(),n,
            StarSampling::ApproximatelyCubicShells);
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity();
        (star,particles,Sph3Config::cold_white_dwarf(2.0).unwrap(),cfg)
    }
    #[test]
    fn force_budget_has_finite_radial_shells_and_mass_fractions() {
        let (star,p,physics,cfg)=fixture(96);
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let report=radial_balance(&p,&mut state,physics,cfg,&star,8);
        assert_eq!(report.shells.len(),8);
        assert!(report.net_force_rms_m_s2.is_finite());
        assert!(report.residual_to_gravity_rms.is_finite());
        assert!(report.pressure_virial_j>0.0);
        assert!(report.gravitational_virial_j<0.0);
        assert!(report.pressure_to_gravity_virial.is_finite());
        assert!(report.shells.iter().all(|s|
            s.particle_count>0&&s.gravity_in_m_s2>0.0&&
            s.hydrostatic_gravity_in_m_s2>0.0&&
            s.net_radial_rms_m_s2.is_finite() &&
            s.tangential_net_rms_m_s2.is_finite()));
        for w in report.shells.windows(2) {
            assert!(w[0].mean_radius_m < w[1].mean_radius_m);
            assert!((w[0].mass_fraction_hi-w[1].mass_fraction_lo).abs()<1e-12);
        }
        assert!((report.shells.last().unwrap().mass_fraction_hi-1.0).abs()<1e-12);
    }
    #[test]
    fn radial_balance_is_invariant_under_uniform_translation() {
        let (star,mut p,physics,cfg)=fixture(64);
        let mut s0=AdaptiveState::new(&p,star.radius_m,cfg);
        let a=radial_balance(&p,&mut s0,physics,cfg,&star,8);
        for particle in &mut p {
            particle.pos_m+=V3::new(-3e8,4e8,6e8);
        }
        let mut s1=AdaptiveState::new(&p,star.radius_m,cfg);
        let b=radial_balance(&p,&mut s1,physics,cfg,&star,8);
        assert!((a.net_force_rms_m_s2-b.net_force_rms_m_s2).abs()/
            a.net_force_rms_m_s2<1e-10);
        assert!((a.pressure_to_gravity_virial-
            b.pressure_to_gravity_virial).abs()<1e-10);
    }
    #[test]
    fn pressure_plus_gravity_reconstructs_net_shell_acceleration() {
        let (star,p,physics,cfg)=fixture(128);
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let report=radial_balance(&p,&mut state,physics,cfg,&star,10);
        for shell in report.shells {
            assert!((shell.pressure_out_m_s2-shell.gravity_in_m_s2-
                shell.net_radial_m_s2).abs()<1e-5);
        }
    }
}
