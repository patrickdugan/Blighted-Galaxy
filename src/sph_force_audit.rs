//! BG-PHY-15: frozen hydrostatic pressure/gravity consistency audit.
//!
//! At unchanged particle positions, compare the adaptive SPH accelerations
//! against the independent spherical 1D reference g(r)=G M(<r)/r^2.
//! A hydrostatic cold barotrope has pressure acceleration +g(r) rhat
//! and gravitational acceleration -g(r) rhat. Every diagnostic is
//! normalized by the RMS reference gravity on the SAME particle positions.
//! This distinguishes radial mean bias from angular (E0/disorder) noise.
use crate::physics::G_SI;
use crate::sph3::{Fluid3,Massive3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use crate::white_dwarf::WhiteDwarfProfile;

#[derive(Debug,Clone,Copy)]
pub struct FrozenForceMetrics {
    pub particle_count:usize,
    pub reference_gravity_rms_m_s2:f64,
    /// RMS vector error of pressure acceleration against +g_ref rhat.
    pub pressure_ref_rms:f64,
    /// RMS vector error of gravity acceleration against -g_ref rhat.
    pub gravity_ref_rms:f64,
    /// RMS vector error of net acceleration against zero.
    pub net_ref_rms:f64,
    /// Mass-weighted mean radial force / mean g_ref (outward pressure).
    pub pressure_radial_ratio:f64,
    /// Mass-weighted mean inward radial gravity / mean g_ref.
    pub gravity_radial_ratio:f64,
    pub pressure_tangent_rms:f64,
    pub gravity_tangent_rms:f64,
    pub net_tangent_rms:f64,
    pub mean_radius_m:f64,
}
#[derive(Debug,Clone,Copy)]
pub struct FrozenForceShell {
    pub mass_fraction_lo:f64,
    pub mass_fraction_hi:f64,
    pub metrics:FrozenForceMetrics,
}
#[derive(Debug,Clone)]
pub struct FrozenForceAudit {
    pub center_m:V3,
    pub global:FrozenForceMetrics,
    pub shells:Vec<FrozenForceShell>,
}
#[derive(Debug,Default,Clone,Copy)]
struct Moments {
    count:usize,
    mass:f64,
    mr:f64,
    m_g:f64,
    m_g2:f64,
    m_p_rad:f64,
    m_g_rad:f64,
    m_dp2:f64,
    m_dg2:f64,
    m_net2:f64,
    m_pt2:f64,
    m_gt2:f64,
    m_nt2:f64,
}
impl Moments {
    fn add(&mut self,m:f64,r:f64,g:f64,radial:V3,
           pressure:V3,gravity:V3,net:V3) {
        let target=radial*g;
        let p_rad=pressure.dot(radial);
        let g_rad=-gravity.dot(radial);
        let n_rad=net.dot(radial);
        let pt=pressure-radial*p_rad;
        let gt=gravity+radial*g_rad;
        let nt=net-radial*n_rad;
        self.count+=1;
        self.mass+=m;
        self.mr+=m*r;
        self.m_g+=m*g;
        self.m_g2+=m*g*g;
        self.m_p_rad+=m*p_rad;
        self.m_g_rad+=m*g_rad;
        self.m_dp2+=m*(pressure-target).len_sq();
        self.m_dg2+=m*(gravity+target).len_sq();
        self.m_net2+=m*net.len_sq();
        self.m_pt2+=m*pt.len_sq();
        self.m_gt2+=m*gt.len_sq();
        self.m_nt2+=m*nt.len_sq();
    }
    fn finish(&self)->FrozenForceMetrics {
        assert!(self.count>0&&self.mass>0.0&&self.m_g2>0.0);
        let inv_ref=1.0/self.m_g2;
        FrozenForceMetrics {
            particle_count:self.count,
            reference_gravity_rms_m_s2:(self.m_g2/self.mass).sqrt(),
            pressure_ref_rms:(self.m_dp2*inv_ref).sqrt(),
            gravity_ref_rms:(self.m_dg2*inv_ref).sqrt(),
            net_ref_rms:(self.m_net2*inv_ref).sqrt(),
            pressure_radial_ratio:self.m_p_rad/self.m_g,
            gravity_radial_ratio:self.m_g_rad/self.m_g,
            pressure_tangent_rms:(self.m_pt2*inv_ref).sqrt(),
            gravity_tangent_rms:(self.m_gt2*inv_ref).sqrt(),
            net_tangent_rms:(self.m_nt2*inv_ref).sqrt(),
            mean_radius_m:self.mr/self.mass,
        }
    }
}
/// Produce a force snapshot without advancing positions or velocities.
/// Force closure refreshes densities and smoothing lengths in state.
pub fn frozen_force_audit(
    particles:&[Fluid3],state:&mut AdaptiveState,
    physics:Sph3Config,cfg:AdaptiveConfig,
    reference:&WhiteDwarfProfile,n_shells:usize,
)->FrozenForceAudit {
    assert!(!particles.is_empty()&&n_shells>0&&n_shells<=particles.len());
    assert!(state.valid_for(particles.len()));
    let mass=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    assert!(mass>0.0);
    let center=particles.iter().fold(V3::default(),
        |sum,p|sum+p.pos_m*p.mass_kg)/mass;
    let empty=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let forces=sph_adaptive::forces(&empty,particles,state,physics,cfg);
    let mut order=particles.iter().enumerate()
        .map(|(i,p)|((p.pos_m-center).len(),i))
        .collect::<Vec<_>>();
    order.sort_by(|a,b|a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut bins=vec![Moments::default();n_shells];
    let mut all=Moments::default();
    let mut cumulative=0.0;
    for (r,i) in order {
        let p=&particles[i];
        let unit=if r>0.0 {(p.pos_m-center)*(1.0/r)}
            else {V3::default()};
        let g=if r>0.0 {
            G_SI*reference.enclosed_mass_at(r)/(r*r)
        } else {0.0};
        let p_force=forces.pressure_accelerations_m_s2[i];
        let g_force=forces.gravity_accelerations_m_s2[i];
        let total=forces.accelerations_m_s2[i];
        let fraction=(cumulative+0.5*p.mass_kg)/mass;
        let bin=((fraction*n_shells as f64).floor() as usize)
            .min(n_shells-1);
        bins[bin].add(p.mass_kg,r,g,unit,p_force,g_force,total);
        all.add(p.mass_kg,r,g,unit,p_force,g_force,total);
        cumulative+=p.mass_kg;
    }
    assert!(bins.iter().all(|b|b.count>0),
        "choose fewer mass-quantile shells for this particle sample");
    FrozenForceAudit {
        center_m:center,
        global:all.finish(),
        shells:bins.iter().enumerate().map(|(i,b)|FrozenForceShell {
            mass_fraction_lo:i as f64/n_shells as f64,
            mass_fraction_hi:(i+1) as f64/n_shells as f64,
            metrics:b.finish(),
        }).collect(),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_massprofile::profile_diagnostics;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn fixture(n:usize,seed:u64)->(WhiteDwarfProfile,Vec<Fluid3>,Sph3Config,AdaptiveConfig) {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let (p,_)=seed_cubic_shells_oriented(
            &star,V3::new(1e8,0.0,0.0),V3::default(),n,seed);
        let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m).with_cubic_gravity();
        (star,p,physics,cfg)
    }
    #[test]
    fn exact_initial_outer_mass_quantile_is_inside_surface() {
        for n in [64,128,256,512,1024] {
            let (star,p,_,_)=fixture(n,17);
            let observed=profile_diagnostics(&p,&star).radius_max_m;
            let expected=star.radius_at_mass_fraction(1.0-1.0/n as f64);
            assert!((observed-expected).abs()/star.radius_m<1e-12);
            assert!(observed<star.radius_m);
        }
    }
    #[test]
    fn frozen_error_components_are_finite_and_shells_cover_mass() {
        let (star,p,phys,cfg)=fixture(96,17);
        let positions=p.iter().map(|x|x.pos_m).collect::<Vec<_>>();
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let a=frozen_force_audit(&p,&mut state,phys,cfg,&star,8);
        assert!(p.iter().zip(&positions).all(|(x,y)|(x.pos_m-*y).len()<1e-12));
        assert_eq!(a.global.particle_count,96);
        assert_eq!(a.shells.iter().map(|s|s.metrics.particle_count).sum::<usize>(),96);
        assert!(a.global.reference_gravity_rms_m_s2>0.0);
        for m in std::iter::once(a.global).chain(a.shells.iter().map(|s|s.metrics)) {
            assert!(m.pressure_ref_rms.is_finite());
            assert!(m.gravity_ref_rms.is_finite());
            assert!(m.net_ref_rms.is_finite());
            assert!(m.pressure_radial_ratio.is_finite());
            assert!(m.gravity_radial_ratio.is_finite());
            assert!(m.pressure_tangent_rms.is_finite());
            assert!(m.gravity_tangent_rms.is_finite());
            assert!(m.net_tangent_rms.is_finite());
        }
        assert!((a.shells[0].mass_fraction_lo).abs()<1e-12);
        assert!((a.shells.last().unwrap().mass_fraction_hi-1.0).abs()<1e-12);
    }
    #[test]
    fn translated_frozen_force_errors_are_invariant() {
        let (star,mut p,phys,cfg)=fixture(64,0);
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let a=frozen_force_audit(&p,&mut state,phys,cfg,&star,8).global;
        for x in &mut p {x.pos_m+=V3::new(3e8,-4e8,2e8);}
        let mut shifted=AdaptiveState::new(&p,star.radius_m,cfg);
        let b=frozen_force_audit(&p,&mut shifted,phys,cfg,&star,8).global;
        assert!((a.pressure_ref_rms-b.pressure_ref_rms).abs()<1e-9);
        assert!((a.gravity_ref_rms-b.gravity_ref_rms).abs()<1e-9);
        assert!((a.net_ref_rms-b.net_ref_rms).abs()<1e-9);
    }
}
