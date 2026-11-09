//! Zero-temperature electron-degeneracy equation of state and hydrostatic
//! white-dwarf structure in SI units. A Newtonian, spherical, nonrotating
//! reference star with mean mass per electron mu_e * atomic_mass_unit.
//!
//! This is appropriate for benchmarking the existing SPH proxy. Ionic
//! thermal pressure, Coulomb corrections, crystallization, composition
//! gradients, neutronization, and general-relativistic stellar structure
//! are excluded.
use crate::physics::{G_SI,SOLAR_MASS_KG};
use std::f64::consts::PI;

pub const ELECTRON_MASS_KG:f64 = 9.109_383_713_9e-31;
pub const ATOMIC_MASS_UNIT_KG:f64 = 1.660_539_068_92e-27;
pub const HBAR_SI:f64 = 1.054_571_817e-34;
pub const C_SI:f64 = 299_792_458.0;

#[derive(Clone,Copy,Debug)]
pub struct DegenerateElectrons {pub mu_e:f64}
impl DegenerateElectrons {
    pub fn new(mu_e:f64)->Option<Self> {
        if mu_e.is_finite() && (1.0..=4.0).contains(&mu_e) {
            Some(Self{mu_e})
        } else {None}
    }
    fn h0(&self)->f64 {
        ELECTRON_MASS_KG*C_SI*C_SI/(self.mu_e*ATOMIC_MASS_UNIT_KG)
    }
    fn rho0(&self)->f64 {
        self.mu_e*ATOMIC_MASS_UNIT_KG*
            (ELECTRON_MASS_KG*C_SI/HBAR_SI).powi(3)/(3.0*PI*PI)
    }
    /// Relativity parameter x=p_F/(m_e c).
    pub fn x_for_density(&self,rho_kg_m3:f64)->f64 {
        if !rho_kg_m3.is_finite()||rho_kg_m3<=0.0 {return 0.0;}
        (rho_kg_m3/self.rho0()).cbrt()
    }
    pub fn rho_for_x(&self,x:f64)->f64 {
        if x<=0.0 {return 0.0;}
        self.rho0()*x.powi(3)
    }
    pub fn enthalpy_m2_s2(&self,rho_kg_m3:f64)->f64 {
        let x=self.x_for_density(rho_kg_m3);
        // Stable evaluation of sqrt(1+x²)-1.
        self.h0()*x*x/((1.0+x*x).sqrt()+1.0)
    }
    pub fn density_from_enthalpy(&self,h_m2_s2:f64)->f64 {
        if h_m2_s2<=0.0||!h_m2_s2.is_finite() {return 0.0;}
        let y=h_m2_s2/self.h0();
        self.rho0()*(y*(2.0+y)).powf(1.5)
    }
    /// Fully relativistic zero-temperature Fermi-electron pressure (Pa).
    /// The x<<1 series prevents catastrophic cancellation in the closed form.
    pub fn pressure_pa(&self,rho_kg_m3:f64)->f64 {
        let x=self.x_for_density(rho_kg_m3);
        if x==0.0 {return 0.0;}
        let a=ELECTRON_MASS_KG.powi(4)*C_SI.powi(5)/
            (3.0*PI*PI*HBAR_SI.powi(3));
        if x<0.1 {
            let x2=x*x;
            // Integral of t^4/sqrt(1+t²) from 0 to x.
            return a*x.powi(5)*(1.0/5.0-x2/14.0+
                x2*x2/24.0-5.0*x2*x2*x2/176.0);
        }
        let f=x*(2.0*x*x-3.0)*(1.0+x*x).sqrt()+
            3.0*x.asinh();
        a*f/8.0
    }
    /// dP/d(rho), m²/s². Sound-speed squared for barotropic fluid.
    pub fn dp_drho_m2_s2(&self,rho_kg_m3:f64)->f64 {
        let x=self.x_for_density(rho_kg_m3);
        self.h0()*x*x/(3.0*(1.0+x*x).sqrt())
    }
    /// Electron kinetic internal energy per mass, J/kg; excludes rest mass.
    pub fn internal_energy_j_kg(&self,rho_kg_m3:f64)->f64 {
        if rho_kg_m3<=0.0 {return 0.0;}
        (self.enthalpy_m2_s2(rho_kg_m3)-
            self.pressure_pa(rho_kg_m3)/rho_kg_m3).max(0.0)
    }
    pub fn nonrel_k_si(&self)->f64 {
        HBAR_SI.powi(2)/(5.0*ELECTRON_MASS_KG)*
            (3.0*PI*PI).powf(2.0/3.0)/
            (self.mu_e*ATOMIC_MASS_UNIT_KG).powf(5.0/3.0)
    }
    pub fn ultrarel_k_si(&self)->f64 {
        HBAR_SI*C_SI/4.0*(3.0*PI*PI).cbrt()/
            (self.mu_e*ATOMIC_MASS_UNIT_KG).powf(4.0/3.0)
    }
}

#[derive(Clone,Copy,Debug)]
pub struct StellarShell {
    pub radius_m:f64,
    pub enclosed_mass_kg:f64,
    pub rho_kg_m3:f64,
    pub enthalpy_m2_s2:f64,
}
#[derive(Clone,Debug)]
pub struct WhiteDwarfProfile {
    pub eos:DegenerateElectrons,
    pub rho_c_kg_m3:f64,
    pub radius_m:f64,
    pub mass_kg:f64,
    /// Monotonically increasing radius and enclosed mass.
    pub shells:Vec<StellarShell>,
}
impl WhiteDwarfProfile {
    pub fn mass_solar(&self)->f64 {self.mass_kg/SOLAR_MASS_KG}
    pub fn enclosed_mass_at(&self,radius_m:f64)->f64 {
        if radius_m<=0.0 {return 0.0;}
        if radius_m>=self.radius_m {return self.mass_kg;}
        let i=self.shells.partition_point(|s|s.radius_m<radius_m);
        if i==0 {return self.shells[0].enclosed_mass_kg;}
        let a=self.shells[i-1];let b=self.shells[i];
        let f=((radius_m-a.radius_m)/(b.radius_m-a.radius_m)).clamp(0.0,1.0);
        a.enclosed_mass_kg+(b.enclosed_mass_kg-a.enclosed_mass_kg)*f
    }
    pub fn density_at(&self,radius_m:f64)->f64 {
        if radius_m>=self.radius_m {return 0.0;}
        if radius_m<=0.0 {return self.rho_c_kg_m3;}
        let i=self.shells.partition_point(|s|s.radius_m<radius_m);
        if i==0 {return self.rho_c_kg_m3;}
        let a=self.shells[i-1];let b=self.shells[i];
        let f=((radius_m-a.radius_m)/(b.radius_m-a.radius_m)).clamp(0.0,1.0);
        a.rho_kg_m3+(b.rho_kg_m3-a.rho_kg_m3)*f
    }
    /// Inverse cumulative mass profile for reproducible equal-mass particles.
    pub fn radius_at_mass_fraction(&self,fraction:f64)->f64 {
        let f=fraction.clamp(0.0,1.0);
        if f==0.0 {return 0.0;}
        if f==1.0 {return self.radius_m;}
        let target=f*self.mass_kg;
        let i=self.shells.partition_point(|s|s.enclosed_mass_kg<target);
        if i==0 {return self.shells[0].radius_m;}
        if i>=self.shells.len() {return self.radius_m;}
        let a=self.shells[i-1];let b=self.shells[i];
        let t=(target-a.enclosed_mass_kg)/(b.enclosed_mass_kg-a.enclosed_mass_kg);
        a.radius_m+(b.radius_m-a.radius_m)*t
    }
}
fn derivatives(r:f64,h:f64,m:f64,eos:DegenerateElectrons)->(f64,f64) {
    let rho=eos.density_from_enthalpy(h);
    (-G_SI*m/(r*r),4.0*PI*r*r*rho)
}
fn rk4_shell(r:f64,h:f64,m:f64,dr:f64,eos:DegenerateElectrons)->(f64,f64) {
    let (a,b)=derivatives(r,h,m,eos);
    let (c,d)=derivatives(r+dr*0.5,h+a*dr*0.5,m+b*dr*0.5,eos);
    let (e,f)=derivatives(r+dr*0.5,h+c*dr*0.5,m+d*dr*0.5,eos);
    let (g,k)=derivatives(r+dr,h+e*dr,m+f*dr,eos);
    (h+dr*(a+2.0*c+2.0*e+g)/6.0,
     m+dr*(b+2.0*d+2.0*f+k)/6.0)
}
/// Newtonian hydrostatic dH/dr=-Gm/r², dm/dr=4πr²ρ(H).
/// A small finite-r center expansion removes the r=0 singularity.
pub fn solve_hydrostatic(
    rho_c_kg_m3:f64,eos:DegenerateElectrons,step_m:f64,
)->Option<WhiteDwarfProfile> {
    if !rho_c_kg_m3.is_finite()||rho_c_kg_m3<=0.0||
        !step_m.is_finite()||!(100.0..=100_000.0).contains(&step_m) {
        return None;
    }
    let r0=(step_m*0.05).min(200.0);
    let initial_enthalpy=eos.enthalpy_m2_s2(rho_c_kg_m3);
    let mut h=initial_enthalpy-2.0*PI*G_SI*rho_c_kg_m3*r0*r0/3.0;
    let mut m=4.0*PI*rho_c_kg_m3*r0.powi(3)/3.0;
    if h<=0.0 {return None;}
    let mut r=r0;
    let mut shells=vec![
        StellarShell{radius_m:0.0,enclosed_mass_kg:0.0,
            rho_kg_m3:rho_c_kg_m3,enthalpy_m2_s2:initial_enthalpy},
        StellarShell{radius_m:r,enclosed_mass_kg:m,
            rho_kg_m3:eos.density_from_enthalpy(h),enthalpy_m2_s2:h},
    ];
    for _ in 0..150_000 {
        let (hn,mn)=rk4_shell(r,h,m,step_m,eos);
        if !hn.is_finite()||!mn.is_finite()||mn<m {return None;}
        if hn<=0.0 {
            let frac=(h/(h-hn)).clamp(0.0,1.0);
            let surface=r+step_m*frac;
            let total_mass=m+(mn-m)*frac;
            shells.push(StellarShell{
                radius_m:surface,enclosed_mass_kg:total_mass,
                rho_kg_m3:0.0,enthalpy_m2_s2:0.0,
            });
            return Some(WhiteDwarfProfile{
                eos,rho_c_kg_m3,radius_m:surface,
                mass_kg:total_mass,shells,
            });
        }
        r+=step_m;h=hn;m=mn;
        shells.push(StellarShell{
            radius_m:r,enclosed_mass_kg:m,
            rho_kg_m3:eos.density_from_enthalpy(h),
            enthalpy_m2_s2:h,
        });
        if r>100.0e6 {return None;}
    }
    None
}
/// Find the cold Newtonian equilibrium profile closest to a requested mass.
/// Uses a bracketed density search and deterministic fixed spatial steps.
pub fn profile_for_mass(
    target_mass_kg:f64,eos:DegenerateElectrons,step_m:f64,
)->Option<WhiteDwarfProfile> {
    if !target_mass_kg.is_finite()||
        target_mass_kg<=0.0 {return None;}
    let mut lo=1e7_f64;let mut hi=1e12_f64;
    if solve_hydrostatic(lo,eos,step_m)?.mass_kg>target_mass_kg||
       solve_hydrostatic(hi,eos,step_m)?.mass_kg<target_mass_kg {
        return None;
    }
    let mut best=None;
    for _ in 0..38 {
        let mid=(lo*hi).sqrt();
        let solution=solve_hydrostatic(mid,eos,step_m)?;
        if solution.mass_kg<target_mass_kg {lo=mid;} else {hi=mid;}
        best=Some(solution);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(x:f64,y:f64,tol:f64) {
        assert!((x-y).abs()<=tol*y.abs(),
            "x={x:e}, expected={y:e}, relative_err={:e}",
            (x-y).abs()/y.abs());
    }
    #[test]
    fn degeneracy_constants_for_mu_e_2() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        close(e.nonrel_k_si(),3.161_128_6e6,1e-5);
        close(e.ultrarel_k_si(),4.934_824_9e9,1e-5);
    }
    #[test]
    fn zero_density_is_zero_pressure_and_energy() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        assert_eq!(e.pressure_pa(0.0),0.0);
        assert_eq!(e.internal_energy_j_kg(0.0),0.0);
    }
    #[test]
    fn pressure_reference_at_one_billion_kg_per_cubic_metre() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        close(e.pressure_pa(1e9),2.627275499e21,1e-6);
    }
    #[test]
    fn pressure_asymptotes_and_sound_speed() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        let rho=1e6;
        close(e.pressure_pa(rho)/(e.nonrel_k_si()*rho.powf(5.0/3.0)),
            1.0,0.01);
        let high=1e14;
        close(e.pressure_pa(high)/(e.ultrarel_k_si()*high.powf(4.0/3.0)),
            1.0,0.01);
        assert!(e.dp_drho_m2_s2(rho)>0.0);
        assert!(e.internal_energy_j_kg(rho)>0.0);
    }
    #[test]
    fn enthalpy_density_roundtrip() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        for rho in [1e6,1e8,1e9,1e11] {
            close(e.density_from_enthalpy(e.enthalpy_m2_s2(rho)),rho,1e-12);
        }
    }
    #[test]
    fn hydrostatic_mass_radius_baselines() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        let a=solve_hydrostatic(1e9,e,10_000.0).unwrap();
        close(a.mass_solar(),0.401,0.015);
        close(a.radius_m/1e6,10.958,0.01);
        let b=solve_hydrostatic(1e10,e,10_000.0).unwrap();
        close(b.mass_solar(),0.815,0.015);
        close(b.radius_m/1e6,7.085,0.01);
    }
    #[test]
    fn central_density_controls_mass_and_radius() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        let a=solve_hydrostatic(1e9,e,20_000.0).unwrap();
        let b=solve_hydrostatic(1e10,e,20_000.0).unwrap();
        assert!(b.mass_kg>a.mass_kg);
        assert!(b.radius_m<a.radius_m);
        assert!(a.shells.windows(2).all(|p|
            p[1].enclosed_mass_kg>=p[0].enclosed_mass_kg &&
            p[1].radius_m>p[0].radius_m));
    }
    #[test]
    fn solve_target_065_solar_mass_and_mass_cdf() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        let p=profile_for_mass(0.65*SOLAR_MASS_KG,e,10_000.0).unwrap();
        close(p.mass_solar(),0.65,1e-6);
        close(p.radius_m/1e6,8.4022,0.01);
        let r_half=p.radius_at_mass_fraction(0.5);
        close(p.enclosed_mass_at(r_half)/p.mass_kg,0.5,1e-9);
        assert!(p.density_at(0.0)>p.density_at(r_half));
        assert_eq!(p.density_at(p.radius_m),0.0);
    }
    #[test]
    fn pressure_derivative_matches_numeric_difference() {
        let e=DegenerateElectrons::new(2.0).unwrap();
        let r=1e9;let h=1e5;
        let numerical=(e.pressure_pa(r+h)-e.pressure_pa(r-h))/(2.0*h);
        close(e.dp_drho_m2_s2(r),numerical,1e-6);
    }
}
