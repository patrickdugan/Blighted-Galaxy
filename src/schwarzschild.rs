//! Exact equatorial timelike test-particle geodesics outside a static,
//! uncharged Schwarzschild black hole, with numerical integration in proper time.
//!
//! Dimensionless geometric units: x = r / r_g, s = c*tau/r_g,
//! r_g = GM/c², L = angular_momentum/(c*r_g), epsilon=E/(mc²).
//! Derivations: Oregon State Geometry of GR, Schwarzschild geodesics and orbits.
//! This module is *not* a Kerr solver, self-gravitating GR, or photon renderer.
use crate::physics::{C_SI,G_SI};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitStatus {
    Active,
    /// Schwarzschild coordinate time diverges as x -> 2.
    /// Stop before the horizon and hand off to a horizon-regular chart later.
    HorizonCoordinateCutoff,
}
#[derive(Debug, Clone, Copy)]
pub struct SchwarzschildOrbit {
    pub mass_kg:f64,
    pub r_g_m:f64,
    /// Areal radius in gravitational radii, r/(GM/c²).
    pub x:f64,
    /// dr/dtau divided by c.
    pub radial_u:f64,
    pub phi_rad:f64,
    pub proper_time_s:f64,
    pub coordinate_time_s:f64,
    /// Conserved specific angular momentum, L/(c r_g).
    pub ell:f64,
    /// Conserved Killing energy, E/(mc²).
    pub epsilon:f64,
    pub status:OrbitStatus,
}

const COORDINATE_CUTOFF_X:f64 = 2.01;
const MAX_DS:f64 = 0.04;

/// Effective potential squared for timelike equatorial Schwarzschild orbits.
pub fn effective_potential_squared(x:f64,ell:f64)->f64 {
    (1.0-2.0/x)*(1.0+ell*ell/(x*x))
}
pub fn gravitational_radius_m(mass_kg:f64)->Option<f64> {
    if !mass_kg.is_finite()||mass_kg<=0.0 {return None;}
    Some(G_SI*mass_kg/C_SI.powi(2))
}
pub fn isco_radius_m(mass_kg:f64)->Option<f64> {
    Some(6.0*gravitational_radius_m(mass_kg)?)
}
pub fn photon_sphere_radius_m(mass_kg:f64)->Option<f64> {
    Some(3.0*gravitational_radius_m(mass_kg)?)
}
pub fn circular_orbit_constants(x:f64)->Option<(f64,f64)> {
    if !x.is_finite()||x<=3.0 {return None;}
    let epsilon=(1.0-2.0/x)/(1.0-3.0/x).sqrt();
    let ell=x/(x-3.0).sqrt();
    Some((epsilon,ell))
}
/// Initial coordinate radial position, proper radial derivative dr/dtau,
/// and specific angular momentum. All SI inputs except normalized constants.
impl SchwarzschildOrbit {
    pub fn from_initial(
        mass_kg:f64, r_m:f64, dr_dproper_s:f64,
        angular_momentum_m2_s:f64,
    )->Option<Self> {
        let rg=gravitational_radius_m(mass_kg)?;
        if !r_m.is_finite()||r_m<=COORDINATE_CUTOFF_X*rg
            ||!dr_dproper_s.is_finite()||!angular_momentum_m2_s.is_finite() {
            return None;
        }
        let x=r_m/rg;
        let ell=angular_momentum_m2_s/(C_SI*rg);
        let radial_u=dr_dproper_s/C_SI;
        let eps=(radial_u*radial_u+effective_potential_squared(x,ell)).sqrt();
        if !eps.is_finite()||eps<=0.0 {return None;}
        Some(Self{mass_kg,r_g_m:rg,x,radial_u,phi_rad:0.0,
            proper_time_s:0.0,coordinate_time_s:0.0,
            ell,epsilon:eps,status:OrbitStatus::Active})
    }
    pub fn circular(mass_kg:f64,r_m:f64)->Option<Self> {
        let rg=gravitational_radius_m(mass_kg)?;
        let x=r_m/rg;
        let (eps,l)=circular_orbit_constants(x)?;
        Self::from_initial(mass_kg,r_m,0.0,l*C_SI*rg)
            .map(|mut orbit|{orbit.epsilon=eps;orbit})
    }
    /// Bound nonplunging orbit from Darwin p,e parameters; p is in r_g.
    /// Stable bound family requires p>6+2e, 0<=e<1.
    pub fn bound_from_p_e(mass_kg:f64,p:f64,e:f64)->Option<Self> {
        if !p.is_finite()||!e.is_finite()||e<0.0||e>=1.0||p<=6.0+2.0*e {
            return None;
        }
        let rg=gravitational_radius_m(mass_kg)?;
        let radius=rg*p/(1.0+e);
        let ell=p/(p-3.0-e*e).sqrt();
        let eps2=((p-2.0-2.0*e)*(p-2.0+2.0*e))/
            (p*(p-3.0-e*e));
        let mut orbit=Self::from_initial(mass_kg,radius,0.0,ell*C_SI*rg)?;
        orbit.epsilon=eps2.sqrt();
        Some(orbit)
    }
    pub fn radius_m(&self)->f64 {self.x*self.r_g_m}
    pub fn xy_m(&self)->crate::V2 {
        crate::V2::new(self.radius_m()*self.phi_rad.cos(),
            self.radius_m()*self.phi_rad.sin())
    }
    pub fn redshift_factor(&self)->f64 {1.0-2.0/self.x}
    pub fn current_energy_constraint_error(&self)->f64 {
        self.radial_u*self.radial_u+
            effective_potential_squared(self.x,self.ell)-self.epsilon*self.epsilon
    }
    /// Exact circular coordinate-time angular frequency in Schwarzschild
    /// spacetime, equal to sqrt(GM/r³) for any circular timelike orbit.
    pub fn circular_coordinate_omega_rad_s(&self)->Option<f64> {
        if self.x<=3.0 {return None;}
        Some((G_SI*self.mass_kg/self.radius_m().powi(3)).sqrt())
    }
    fn derivative(&self,y:[f64;4])->[f64;4] {
        let x=y[0];
        let ell2=self.ell*self.ell;
        [
            y[1],
            -1.0/(x*x)+ell2/x.powi(3)-3.0*ell2/x.powi(4),
            self.ell/(x*x),
            self.epsilon/(1.0-2.0/x),
        ]
    }
    fn rk4(&self,y:[f64;4],h:f64)->[f64;4] {
        let a=self.derivative(y);
        let yb=std::array::from_fn(|i|y[i]+a[i]*h*0.5);
        let b=self.derivative(yb);
        let yc=std::array::from_fn(|i|y[i]+b[i]*h*0.5);
        let c=self.derivative(yc);
        let yd=std::array::from_fn(|i|y[i]+c[i]*h);
        let d=self.derivative(yd);
        std::array::from_fn(|i|y[i]+h*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0)
    }
    /// Advance an active timelike orbit by proper time. Returns false after
    /// reaching the coordinate-chart cutoff. Numerical steps are capped in
    /// geometric units to control truncation error.
    pub fn advance_proper(&mut self,dt_tau_s:f64)->Result<(),&'static str> {
        if !dt_tau_s.is_finite()||dt_tau_s<0.0 {return Err("invalid proper time");}
        if self.status!=OrbitStatus::Active||dt_tau_s==0.0 {return Ok(());}
        let mut remaining=dt_tau_s*C_SI/self.r_g_m;
        if remaining>1e6 {return Err("proper interval too large");}
        while remaining>0.0 {
            if self.x<=COORDINATE_CUTOFF_X {
                self.status=OrbitStatus::HorizonCoordinateCutoff;
                break;
            }
            let mut h=remaining.min(MAX_DS);
            // Near the Schwarzschild chart's horizon, limit radial travel
            // per step to avoid evaluating t derivatives at x <= 2.
            if self.radial_u<0.0 {
                h=h.min((self.x-COORDINATE_CUTOFF_X)/
                    (2.0*self.radial_u.abs().max(1e-8)));
            }
            if h<=1e-14 {self.status=OrbitStatus::HorizonCoordinateCutoff;break;}
            let y=[self.x,self.radial_u,self.phi_rad,0.0];
            let next=self.rk4(y,h);
            if next[0]<=COORDINATE_CUTOFF_X||!next.iter().all(|v|v.is_finite()){
                // Retry a smaller interval on the next iteration.
                if h>1e-9 {
                    let reduced=h*0.25;
                    let fine=self.rk4(y,reduced);
                    if fine[0]<=COORDINATE_CUTOFF_X ||
                        !fine.iter().all(|v|v.is_finite()){
                        self.status=OrbitStatus::HorizonCoordinateCutoff;
                        break;
                    }
                    self.x=fine[0];self.radial_u=fine[1];self.phi_rad=fine[2];
                    self.proper_time_s+=reduced*self.r_g_m/C_SI;
                    self.coordinate_time_s+=fine[3]*self.r_g_m/C_SI;
                    remaining-=reduced;
                    continue;
                } else {
                    self.status=OrbitStatus::HorizonCoordinateCutoff;
                    break;
                }
            }
            self.x=next[0];
            self.radial_u=next[1];
            self.phi_rad=next[2];
            self.coordinate_time_s+=next[3]*self.r_g_m/C_SI;
            self.proper_time_s+=h*self.r_g_m/C_SI;
            remaining-=h;
        }
        Ok(())
    }
}
pub fn periapsis_advance_weak_rad(p_in_rg:f64)->Option<f64> {
    if !p_in_rg.is_finite()||p_in_rg<=6.0 {return None;}
    Some(6.0*PI/p_in_rg)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::{SOLAR_MASS_KG,C_SI};
    fn near(x:f64,y:f64,relative:f64) {
        assert!((x-y).abs()<=relative*y.abs().max(1e-25),
            "actual={x:e} expected={y:e}");
    }
    #[test]
    fn schwarzschild_scales_and_special_radii() {
        let m=10.0*SOLAR_MASS_KG;
        let rg=gravitational_radius_m(m).unwrap();
        near(isco_radius_m(m).unwrap()/rg,6.0,1e-14);
        near(photon_sphere_radius_m(m).unwrap()/rg,3.0,1e-14);
        assert!(SchwarzschildOrbit::circular(m,2.99*rg).is_none());
    }
    #[test]
    fn circular_orbit_constants_at_isco() {
        let (e,l)=circular_orbit_constants(6.0).unwrap();
        near(e,(8.0/9.0f64).sqrt(),1e-14);
        near(l,12.0f64.sqrt(),1e-14);
    }
    #[test]
    fn circular_omega_is_newtonian_coordinate_frequency() {
        let m=10.0*SOLAR_MASS_KG;let rg=gravitational_radius_m(m).unwrap();
        let o=SchwarzschildOrbit::circular(m,10.0*rg).unwrap();
        near(o.circular_coordinate_omega_rad_s().unwrap(),
            (G_SI*m/(10.0*rg).powi(3)).sqrt(),1e-14);
    }
    #[test]
    fn stable_circular_orbit_remains_stable_for_one_revolution() {
        let m=10.0*SOLAR_MASS_KG;let rg=gravitational_radius_m(m).unwrap();
        let mut o=SchwarzschildOrbit::circular(m,10.0*rg).unwrap();
        let omega_s=o.ell/(o.x*o.x)*C_SI/rg;
        let proper_period=2.0*PI/omega_s;
        for _ in 0..10_000 {o.advance_proper(proper_period/10_000.0).unwrap();}
        near(o.x,10.0,1e-8);
        near(o.phi_rad,2.0*PI,1e-6);
        assert!(o.current_energy_constraint_error().abs()<1e-9);
        near(o.coordinate_time_s/proper_period,
            o.epsilon/(1.0-0.2),1e-6);
    }
    #[test]
    fn bound_orbit_exhibits_periapsis_advance() {
        let m=10.0*SOLAR_MASS_KG;
        let rg=gravitational_radius_m(m).unwrap();
        let mut o=SchwarzschildOrbit::bound_from_p_e(m,30.0,0.2).unwrap();
        near(o.radius_m(),25.0*rg,1e-12);
        let dt=0.045*rg/C_SI;
        let mut seen_apo=false;
        let mut advance=None;
        for _ in 0..90_000 {
            let was=o.radial_u;
            o.advance_proper(dt).unwrap();
            if was>0.0&&o.radial_u<=0.0 {seen_apo=true;}
            if seen_apo&&was<0.0&&o.radial_u>=0.0 {
                advance=Some(o.phi_rad-2.0*PI);break;
            }
        }
        let delta=advance.expect("second periapsis within integration window");
        assert!(delta>0.4 && delta<1.5,
            "expected relativistic periapsis advance, got {delta}");
        assert!(o.current_energy_constraint_error().abs()<1e-7);
    }
    #[test]
    fn radial_plunge_approaches_horizon_and_stops_before_coordinate_singularity() {
        let m=10.0*SOLAR_MASS_KG;let rg=gravitational_radius_m(m).unwrap();
        let mut o=SchwarzschildOrbit::from_initial(
            m,5.0*rg,-0.1*C_SI,0.0).unwrap();
        let dt=0.02*rg/C_SI;
        for _ in 0..20_000 {
            o.advance_proper(dt).unwrap();
            if o.status!=OrbitStatus::Active {break;}
        }
        assert_eq!(o.status,OrbitStatus::HorizonCoordinateCutoff);
        assert!(o.x>2.0 && o.coordinate_time_s>o.proper_time_s);
    }
    #[test]
    fn bound_invalid_separatrix_is_rejected() {
        assert!(SchwarzschildOrbit::bound_from_p_e(10.0*SOLAR_MASS_KG,6.2,0.2).is_none());
    }
}
