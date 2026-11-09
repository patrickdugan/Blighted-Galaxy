//! Two-dimensional self-gravitating tidal debris, with symmetric pair forces.
//!
//! This is a 2-D SPH-like model: surface density (kg/m²), an isothermal
//! P=c_s²*surface_density closure, cubic-spline kernel, softened pair gravity,
//! and velocity-Verlet time stepping. It is NOT a 3-D stellar EOS or full SPH
//! for degenerate matter, so "white dwarf shearing" remains a gameplay proxy.
//!
//! The compact central body participates in gravity with equal/opposite force.
//! Total linear momentum is conserved to numerical roundoff. Thermal energy
//! is not tracked because the simplified closure is isothermal.
use crate::physics::{G_SI,SiBody};
use crate::V2;
use std::f64::consts::PI;

#[derive(Clone,Copy,Debug)]
pub struct FluidParcel {
    pub pos_m:V2,
    pub vel_m_s:V2,
    pub mass_kg:f64,
}
#[derive(Clone,Copy,Debug)]
pub struct DebrisConfig {
    pub h_m:f64,
    pub gravity_softening_m:f64,
    pub sound_speed_m_s:f64,
}
impl DebrisConfig {
    pub fn white_dwarf_proxy()->Self {
        Self {h_m:3.0e6,gravity_softening_m:1.0e6,sound_speed_m_s:5.0e5}
    }
    pub fn valid(&self)->bool {
        self.h_m.is_finite()&&self.h_m>0.0&&
        self.gravity_softening_m.is_finite()&&self.gravity_softening_m>0.0&&
        self.sound_speed_m_s.is_finite()&&self.sound_speed_m_s>=0.0
    }
}
/// The 2D cubic-spline normalized kernel is compactly supported at 2h.
/// W is 1/m² and its radial derivative is 1/m³.
pub fn cubic_spline_2d(distance_m:f64,h_m:f64)->f64 {
    let q=distance_m/h_m;
    let a=10.0/(7.0*PI*h_m*h_m);
    if q<0.0||!q.is_finite() {0.0}
    else if q<1.0 {a*(1.0-1.5*q*q+0.75*q*q*q)}
    else if q<2.0 {a*0.25*(2.0-q).powi(3)}
    else {0.0}
}
pub fn spline_derivative_2d(distance_m:f64,h_m:f64)->f64 {
    let q=distance_m/h_m;
    let a=10.0/(7.0*PI*h_m*h_m*h_m);
    if q<0.0||!q.is_finite() {0.0}
    else if q<1.0 {a*(-3.0*q+2.25*q*q)}
    else if q<2.0 {-a*0.75*(2.0-q).powi(2)}
    else {0.0}
}
pub fn surface_densities_kg_m2(particles:&[FluidParcel],h_m:f64)->Vec<f64> {
    assert!(h_m.is_finite()&&h_m>0.0);
    let mut rho:Vec<f64>=particles.iter()
        .map(|p|p.mass_kg*cubic_spline_2d(0.0,h_m)).collect();
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let d=(particles[i].pos_m-particles[j].pos_m).len();
            let w=cubic_spline_2d(d,h_m);
            rho[i]+=particles[j].mass_kg*w;
            rho[j]+=particles[i].mass_kg*w;
        }
    }
    rho
}
#[derive(Debug,Clone)]
pub struct ParticleForces {
    pub body_m_s2:V2,
    pub parcel_m_s2:Vec<V2>,
    pub densities_kg_m2:Vec<f64>,
}
/// Symmetric pairwise pressure and softened gravity, including central recoil.
pub fn accelerations(
    central:&SiBody,particles:&[FluidParcel],cfg:DebrisConfig,
)->ParticleForces {
    assert!(cfg.valid());
    let n=particles.len();
    let rho=surface_densities_kg_m2(particles,cfg.h_m);
    let mut a=vec![V2::default();n];
    let mut central_a=V2::default();
    let eps2=cfg.gravity_softening_m.powi(2);

    for i in 0..n {
        let delta=particles[i].pos_m-central.pos_m;
        let r2=delta.len_sq()+eps2;
        let common=delta*(G_SI/(r2*r2.sqrt()));
        a[i]-=common*central.mass_kg;
        central_a+=common*particles[i].mass_kg;
    }
    for i in 0..n {
        for j in i+1..n {
            let delta=particles[j].pos_m-particles[i].pos_m;
            let r2=delta.len_sq()+eps2;
            let grav=delta*(G_SI/(r2*r2.sqrt()));
            a[i]+=grav*particles[j].mass_kg;
            a[j]-=grav*particles[i].mass_kg;

            let r=delta.len();
            if r>0.0 && r<2.0*cfg.h_m {
                // Position derivative ∇_i W_ij points from j to i.
                // SPH pressure acceleration has a symmetric force pair.
                let d_w_dr=spline_derivative_2d(r,cfg.h_m);
                let pi_over_rho2=cfg.sound_speed_m_s.powi(2)/rho[i];
                let pj_over_rho2=cfg.sound_speed_m_s.powi(2)/rho[j];
                let pressure_factor=-(pi_over_rho2+pj_over_rho2)*d_w_dr;
                let outward=(particles[i].pos_m-particles[j].pos_m)/r;
                let pair=outward*pressure_factor;
                a[i]+=pair*particles[j].mass_kg;
                a[j]-=pair*particles[i].mass_kg;
            }
        }
    }
    ParticleForces{body_m_s2:central_a,parcel_m_s2:a,densities_kg_m2:rho}
}
pub fn system_momentum_kg_m_s(central:&SiBody,particles:&[FluidParcel])->V2 {
    particles.iter().fold(central.vel_m_s*central.mass_kg,
        |sum,p|sum+p.vel_m_s*p.mass_kg)
}
pub fn system_mass_kg(central:&SiBody,particles:&[FluidParcel])->f64 {
    central.mass_kg+particles.iter().map(|p|p.mass_kg).sum::<f64>()
}
/// Paired velocity-Verlet: central gravitational recoil is integrated.
/// Requires a sufficiently small physical timestep for close passages.
pub fn advance(
    central:&mut SiBody,particles:&mut [FluidParcel],cfg:DebrisConfig,dt_s:f64,
) {
    assert!(cfg.valid()&&dt_s.is_finite()&&dt_s>=0.0);
    if dt_s==0.0 {return;}
    let before=accelerations(central,particles,cfg);
    central.pos_m+=central.vel_m_s*dt_s+before.body_m_s2*(0.5*dt_s*dt_s);
    central.vel_m_s+=before.body_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(before.parcel_m_s2) {
        p.pos_m+=p.vel_m_s*dt_s+a*(0.5*dt_s*dt_s);
        p.vel_m_s+=a*(0.5*dt_s);
    }
    let after=accelerations(central,particles,cfg);
    central.vel_m_s+=after.body_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(after.parcel_m_s2) {
        p.vel_m_s+=a*(0.5*dt_s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::{C_SI,SOLAR_MASS_KG};
    fn fixture()->(SiBody,Vec<FluidParcel>,DebrisConfig) {
        let central=SiBody::new(10.0*SOLAR_MASS_KG,3e4,
            V2::new(0.0,0.0),V2::new(200.0,-100.0));
        let mut parts=Vec::new();
        for i in 0..24 {
            let theta=2.0*PI*(i as f64+0.5)/24.0;
            let dir=V2::from_angle(theta);
            parts.push(FluidParcel{
                pos_m:V2::new(1.0e8,0.0)+dir*6e6,
                vel_m_s:V2::new(-3e5,4e4)+dir*2e4,
                mass_kg:0.65*SOLAR_MASS_KG/24.0,
            });
        }
        (central,parts,DebrisConfig::white_dwarf_proxy())
    }
    #[test]
    fn kernel_positive_compact_support_normalized_in_2d() {
        let h=1e6;
        assert!(cubic_spline_2d(0.0,h)>0.0);
        assert_eq!(cubic_spline_2d(2.01*h,h),0.0);
        // Trapezoidal integration of radial 2*pi*r*W(r) ≈ 1.
        let mut integral=0.0;
        for k in 0..10_000 {
            let r=(k as f64+0.5)*2.0*h/10_000.0;
            integral+=2.0*PI*r*cubic_spline_2d(r,h)*(2.0*h/10_000.0);
        }
        assert!((integral-1.0).abs()<1e-6,"integral={integral}");
    }
    #[test]
    fn density_is_strictly_positive() {
        let (_,parts,cfg)=fixture();
        assert!(surface_densities_kg_m2(&parts,cfg.h_m)
            .iter().all(|x|x.is_finite()&&*x>0.0));
    }
    #[test]
    fn pressure_derivative_points_downhill() {
        let h=1e6;
        assert!(spline_derivative_2d(h,h)<0.0);
        assert_eq!(spline_derivative_2d(2.1*h,h),0.0);
    }
    #[test]
    fn pairwise_forces_have_zero_net_force() {
        let (central,parts,cfg)=fixture();
        let f=accelerations(&central,&parts,cfg);
        let net=f.parcel_m_s2.iter().zip(&parts).fold(
            f.body_m_s2*central.mass_kg,
            |sum,(a,p)|sum+*a*p.mass_kg);
        let scale=(f.parcel_m_s2[0].len()*parts[0].mass_kg).max(1.0);
        assert!(net.len()/scale<1e-10,"net / scale={}",net.len()/scale);
    }
    #[test]
    fn isolated_softened_gravity_is_attractive() {
        let central=SiBody::new(10.0*SOLAR_MASS_KG,3e4,
            V2::default(),V2::default());
        let p=FluidParcel{pos_m:V2::new(1e8,0.0),
            vel_m_s:V2::default(),mass_kg:1e26};
        let force=accelerations(&central,&[p],DebrisConfig::white_dwarf_proxy());
        assert!(force.parcel_m_s2[0].x<0.0);
        assert!(force.body_m_s2.x>0.0);
    }
    #[test]
    fn mass_and_total_momentum_remain_constant_integration() {
        let (mut central,mut parts,cfg)=fixture();
        let mass=system_mass_kg(&central,&parts);
        let momentum=system_momentum_kg_m_s(&central,&parts);
        for _ in 0..200 {advance(&mut central,&mut parts,cfg,0.01);}
        assert!((system_mass_kg(&central,&parts)-mass)/mass<1e-14);
        let delta=(system_momentum_kg_m_s(&central,&parts)-momentum).len();
        let scale=parts[0].mass_kg*parts[0].vel_m_s.len();
        assert!(delta/scale<1e-8,
            "relative momentum drift {}",delta/scale);
    }
    #[test]
    fn zero_step_preserves_positions_and_velocities() {
        let (mut central,mut parts,cfg)=fixture();
        let original=central.pos_m;
        let other=parts[0].vel_m_s;
        advance(&mut central,&mut parts,cfg,0.0);
        assert_eq!(central.pos_m,original);
        assert_eq!(parts[0].vel_m_s,other);
    }
    #[test]
    fn self_gravity_has_finite_tracer_speeds_over_short_run() {
        let (mut central,mut parts,cfg)=fixture();
        for _ in 0..100 {advance(&mut central,&mut parts,cfg,0.015);}
        assert!(parts.iter().all(|p|p.vel_m_s.len().is_finite()
            &&p.vel_m_s.len()<C_SI*0.1));
    }
}
