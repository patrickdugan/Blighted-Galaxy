//! Equatorial Kerr test-particle and photon reference calculations.
 //!
 //! Natural units G=M=c=1 in the equations, with radius r/r_g,
 //! spin a*=Jc/(GM²), and time ct/r_g (r_g=GM/c²).
 //! Boyer-Lindquist coordinates: valid outside r_+; this is not a
 //! horizon-penetrating chart and only equatorial trajectories are modeled.
 //!
 //! Radial photon motion uses a second-order affine equation to traverse
 //! turning points without changing a sign manually.
use crate::physics::{C_SI, G_SI};
use crate::V2;
use std::f64::consts::PI;

const MAX_AFFINE_STEP: f64 = 0.075;
const HORIZON_MARGIN: f64 = 0.02;

#[derive(Clone, Copy, Debug)]
pub struct KerrGeometry {
    pub mass_kg: f64,
    /// Signed dimensionless spin; restricted away from extremality.
    pub spin: f64,
    pub r_g_m: f64,
}
impl KerrGeometry {
    pub fn new(mass_kg: f64, spin: f64) -> Option<Self> {
        if !mass_kg.is_finite() || mass_kg <= 0.0 ||
            !spin.is_finite() || spin.abs() > 0.998 {
            return None;
        }
        Some(Self { mass_kg, spin, r_g_m: G_SI * mass_kg / C_SI.powi(2) })
    }
    pub fn outer_horizon_rg(&self) -> f64 {
        1.0 + (1.0 - self.spin * self.spin).sqrt()
    }
    pub fn ergosphere_equatorial_rg(&self) -> f64 { 2.0 }
    pub fn delta(&self, r: f64) -> f64 {
        r * r - 2.0 * r + self.spin * self.spin
    }
    /// Zero-angular-momentum observer frame dragging angular frequency
    /// in Boyer-Lindquist coordinate time (radians / SI second).
    pub fn zamo_omega_rad_s(&self, radius_rg: f64) -> Option<f64> {
        if !radius_rg.is_finite() || radius_rg <= self.outer_horizon_rg() {
            return None;
        }
        let a = self.spin;
        let big_a = (radius_rg * radius_rg + a*a).powi(2) -
            a*a * self.delta(radius_rg);
        Some((2.0*a*radius_rg / big_a) * C_SI/self.r_g_m)
    }
    pub fn isco_rg(&self, direction: i8) -> Option<f64> {
        if direction != 1 && direction != -1 { return None; }
        let a = self.spin;
        let z1 = 1.0 + (1.0 - a*a).cbrt() *
            ((1.0+a).cbrt()+(1.0-a).cbrt());
        let z2 = (3.0*a*a + z1*z1).sqrt();
        let dir = f64::from(direction);
        Some(3.0+z2-dir*( (3.0-z1)*(3.0+z1+2.0*z2) ).sqrt())
    }
    pub fn photon_circular_radius_rg(&self, direction: i8) -> Option<f64> {
        if direction != 1 && direction != -1 { return None; }
        let dir = f64::from(direction);
        Some(2.0*(1.0 + ((2.0/3.0) * (-dir*self.spin).acos()).cos()))
    }
    pub fn photon_critical_b_rg(&self, direction: i8) -> Option<f64> {
        let r = self.photon_circular_radius_rg(direction)?;
        let s = f64::from(direction);
        let root = self.delta(r).sqrt();
        let denom = self.spin + s * root;
        if denom.abs() < 1e-12 { return None; }
        Some((r*r + self.spin*self.spin + s*self.spin*root)/denom)
    }
    /// Equatorial circular timelike geodesic E, L and Omega=dφ/dt.
    /// Direction +1 is positive L, -1 negative L.
    pub fn circular_timelike(&self, r: f64, direction: i8) -> Option<CircularTimelike> {
        if direction != 1 && direction != -1 || !r.is_finite() ||
            r <= self.photon_circular_radius_rg(direction)? { return None; }
        let s = f64::from(direction);
        let sqrt_r = r.sqrt();
        let root = r.powf(1.5) - 3.0*sqrt_r + 2.0*s*self.spin;
        if root<=0.0 {return None;}
        let denom = r.powf(0.75) * root.sqrt();
        let e = (r.powf(1.5)-2.0*sqrt_r+s*self.spin) / denom;
        let l = s*(r*r-2.0*s*self.spin*sqrt_r+self.spin*self.spin)/denom;
        let omega = s/(r.powf(1.5)+s*self.spin)*C_SI/self.r_g_m;
        Some(CircularTimelike{radius_rg:r,energy:e,
            angular_momentum_rg:l,omega_rad_s:omega})
    }
}
#[derive(Debug, Clone, Copy)]
pub struct CircularTimelike {
    pub radius_rg: f64,
    pub energy: f64,
    pub angular_momentum_rg: f64,
    pub omega_rad_s: f64,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhotonStatus { Active, Escaped, HorizonChartCutoff }
#[derive(Debug, Clone, Copy)]
pub struct KerrPhoton {
    pub geometry: KerrGeometry,
    /// Photon energy is normalized to E=1. b=L/E in gravitational radii.
    pub impact_b_rg: f64,
    pub radius_rg: f64,
    /// dr/dlambda with lambda an affine parameter in gravitational units.
    pub radial_velocity: f64,
    pub phi_rad: f64,
    /// Boyer-Lindquist coordinate time in units r_g/c.
    pub t_geometric: f64,
    pub affine_elapsed: f64,
    pub status: PhotonStatus,
    pub outgoing_limit_rg: f64,
}
impl KerrPhoton {
    pub fn inbound(geometry:KerrGeometry,
                   radius_rg:f64, impact_b_rg:f64) -> Option<Self> {
        if !radius_rg.is_finite() || !impact_b_rg.is_finite() ||
            radius_rg <= geometry.outer_horizon_rg() + HORIZON_MARGIN {
            return None;
        }
        let speed_sq = radial_speed_squared(geometry.spin, impact_b_rg, radius_rg);
        if speed_sq < 0.0 || !speed_sq.is_finite() {return None;}
        Some(Self {
            geometry,impact_b_rg,radius_rg,radial_velocity:-speed_sq.sqrt(),
            phi_rad:PI,t_geometric:0.0,affine_elapsed:0.0,
            status:PhotonStatus::Active,outgoing_limit_rg:radius_rg,
        })
    }
    pub fn pos_rg(&self) -> V2 {
        V2::from_angle(self.phi_rad)*self.radius_rg
    }
    pub fn pos_m(&self) -> V2 {
        self.pos_rg() * self.geometry.r_g_m
    }
    pub fn radial_constraint_residual(&self) -> f64 {
        self.radial_velocity.powi(2) -
            radial_speed_squared(self.geometry.spin,self.impact_b_rg,self.radius_rg)
    }
    fn deriv(&self, state:[f64;4]) -> [f64;4] {
        let r = state[0];
        let a = self.geometry.spin;
        let b = self.impact_b_rg;
        let delta = self.geometry.delta(r);
        let p = r*r+a*a-a*b;
        let phi_dot = ((b-a)+a*p/delta)/(r*r);
        let t_dot = (a*(b-a)+(r*r+a*a)*p/delta)/(r*r);
        let rddot = (b*b-a*a)/(r*r*r) -
            3.0*(b-a).powi(2)/r.powi(4);
        [state[1],rddot,phi_dot,t_dot]
    }
    fn rk4(&self, y:[f64;4], h:f64) -> [f64;4] {
        let k1=self.deriv(y);
        let y2=std::array::from_fn(|i|y[i]+h*k1[i]*0.5);
        let k2=self.deriv(y2);
        let y3=std::array::from_fn(|i|y[i]+h*k2[i]*0.5);
        let k3=self.deriv(y3);
        let y4=std::array::from_fn(|i|y[i]+h*k3[i]);
        let k4=self.deriv(y4);
        std::array::from_fn(|i|y[i]+
            h*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0)
    }
    pub fn advance(&mut self, interval_geometric:f64) -> Result<(), &'static str> {
        if !interval_geometric.is_finite() || interval_geometric<0.0 ||
            interval_geometric>1e5 {
            return Err("invalid affine interval");
        }
        if self.status!=PhotonStatus::Active {return Ok(());}
        let mut remaining=interval_geometric;
        while remaining>1e-13 {
            let cutoff=self.geometry.outer_horizon_rg()+HORIZON_MARGIN;
            if self.radius_rg<=cutoff {
                self.status=PhotonStatus::HorizonChartCutoff;break;
            }
            let mut h=remaining.min(MAX_AFFINE_STEP);
            if self.radial_velocity<0.0 {
                h=h.min(0.4*(self.radius_rg-cutoff)/
                    self.radial_velocity.abs().max(1e-12));
            }
            if h<1e-11 {self.status=PhotonStatus::HorizonChartCutoff;break;}
            let y=[self.radius_rg,self.radial_velocity,self.phi_rad,0.0];
            let next=self.rk4(y,h);
            if next[0]<=cutoff || !next.iter().all(|v|v.is_finite()) {
                self.status=PhotonStatus::HorizonChartCutoff;
                break;
            }
            self.radius_rg=next[0];
            self.radial_velocity=next[1];
            self.phi_rad=next[2];
            self.t_geometric+=next[3];
            self.affine_elapsed+=h;
            remaining-=h;
            if self.radius_rg>=self.outgoing_limit_rg &&
                self.radial_velocity>0.0 {
                self.status=PhotonStatus::Escaped;
                break;
            }
        }
        Ok(())
    }
}
pub fn radial_speed_squared(a:f64,b:f64,r:f64)->f64 {
    1.0+(a*a-b*b)/(r*r)+
        2.0*(b-a).powi(2)/(r*r*r)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    fn near(a:f64,b:f64,tol:f64) {
        assert!((a-b).abs()<tol,
            "actual {a:.12e}, expected {b:.12e}, tolerance {tol:.4e}");
    }
    #[test]
    fn kerr_zero_spin_reduces_to_schwarzschild_special_radii() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.0).unwrap();
        near(bh.outer_horizon_rg(),2.0,1e-14);
        near(bh.isco_rg(1).unwrap(),6.0,1e-13);
        near(bh.isco_rg(-1).unwrap(),6.0,1e-13);
        near(bh.photon_circular_radius_rg(1).unwrap(),3.0,1e-14);
        near(bh.photon_critical_b_rg(1).unwrap(),3.0*3.0f64.sqrt(),1e-12);
        near(bh.photon_critical_b_rg(-1).unwrap(),-3.0*3.0f64.sqrt(),1e-12);
        near(bh.zamo_omega_rad_s(10.0).unwrap(),0.0,1e-14);
    }
    #[test]
    fn prograde_isco_moves_inward_retrograde_outward() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.9).unwrap();
        near(bh.isco_rg(1).unwrap(),2.320883,1e-5);
        near(bh.isco_rg(-1).unwrap(),8.717352,1e-5);
        assert!(bh.photon_circular_radius_rg(1).unwrap()<3.0);
        assert!(bh.photon_circular_radius_rg(-1).unwrap()>3.0);
        assert!(bh.photon_critical_b_rg(1).unwrap()<3.0*3.0_f64.sqrt());
        assert!(bh.photon_critical_b_rg(-1).unwrap().abs()>3.0*3.0_f64.sqrt());
    }
    #[test]
    fn kerr_zamo_frame_dragging_positive_and_decreases_far_away() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.8).unwrap();
        let near_omega=bh.zamo_omega_rad_s(3.0).unwrap();
        let far_omega=bh.zamo_omega_rad_s(10.0).unwrap();
        assert!(near_omega>far_omega && far_omega>0.0);
        let reversed=KerrGeometry::new(10.0*SOLAR_MASS_KG,-0.8).unwrap();
        near(reversed.zamo_omega_rad_s(3.0).unwrap(),-near_omega,1e-12);
    }
    #[test]
    fn circular_timelike_kerr_reduces_to_schwarzschild() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.0).unwrap();
        let orbit=bh.circular_timelike(10.0,1).unwrap();
        near(orbit.energy,(1.0-2.0/10.0)/(1.0-3.0/10.0).sqrt(),1e-12);
        near(orbit.angular_momentum_rg,10.0/(10.0-3.0).sqrt(),1e-12);
        near(orbit.omega_rad_s,10.0f64.powf(-1.5)*C_SI/bh.r_g_m,1e-12);
    }
    #[test]
    fn equatorial_null_radial_constraint_holds_while_falling() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.65).unwrap();
        let mut ray=KerrPhoton::inbound(bh,45.0,2.0).unwrap();
        for _ in 0..250 {
            ray.advance(0.05).unwrap();
            if ray.status!=PhotonStatus::Active {break;}
            assert!(ray.radial_constraint_residual().abs()<1e-6);
        }
        assert!(ray.phi_rad.is_finite() && ray.t_geometric>0.0);
    }
    #[test]
    fn schwarzschild_radial_photon_plunges_with_zero_phi_deflection() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.0).unwrap();
        let mut ray=KerrPhoton::inbound(bh,10.0,0.0).unwrap();
        for _ in 0..600 {
            ray.advance(0.05).unwrap();
            if ray.status!=PhotonStatus::Active {break;}
        }
        assert_eq!(ray.status,PhotonStatus::HorizonChartCutoff);
        near(ray.phi_rad,PI,1e-12);
    }
    #[test]
    fn b_above_critical_turns_and_escapes() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.0).unwrap();
        let mut ray=KerrPhoton::inbound(bh,45.0,7.0).unwrap();
        let mut turned=false;
        for _ in 0..6000 {
            ray.advance(0.05).unwrap();
            if ray.radial_velocity>0.0 {turned=true;}
            if ray.status!=PhotonStatus::Active {break;}
        }
        assert!(turned);
        assert_eq!(ray.status,PhotonStatus::Escaped);
        assert!(ray.radial_constraint_residual().abs()<1e-5);
    }
    #[test]
    fn spin_changes_prograde_and_retrograde_light_paths() {
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.75).unwrap();
        let mut pro=KerrPhoton::inbound(bh,65.0,7.0).unwrap();
        let mut retro=KerrPhoton::inbound(bh,65.0,-7.0).unwrap();
        for _ in 0..8000 {
            pro.advance(0.05).unwrap();
            retro.advance(0.05).unwrap();
            if pro.status!=PhotonStatus::Active &&
                retro.status!=PhotonStatus::Active {break;}
        }
        assert_eq!(pro.status,PhotonStatus::Escaped);
        assert_eq!(retro.status,PhotonStatus::Escaped);
        // Compare absolute bending, discounting the direction sign.
        assert!((pro.phi_rad.abs()-retro.phi_rad.abs()).abs()>0.005);
    }
    #[test]
    fn invalid_spin_and_incoming_potential_rejected() {
        assert!(KerrGeometry::new(10.0*SOLAR_MASS_KG,1.01).is_none());
        assert!(KerrGeometry::new(0.0,0.2).is_none());
        let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,0.0).unwrap();
        assert!(KerrPhoton::inbound(bh,2.01,0.0).is_none());
    }
}
