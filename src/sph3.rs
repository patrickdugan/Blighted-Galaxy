//! Three-dimensional self-gravitating SPH-like stellar debris prototype.
 //!
 //! SI units throughout. A 3D M4 cubic kernel, fixed smoothing length,
 //! pair-symmetric Plummer gravity, and an isentropic 5/3 polytropic closure.
 //! This represents an idealized fluid; it does not solve nuclear burning,
 //! radiation transport, degenerate-electron relativity, or strong-field GR.
use crate::physics::{G_SI, SiBody};
use crate::white_dwarf::{DegenerateElectrons,WhiteDwarfProfile};
use std::f64::consts::PI;
use std::ops::{Add,AddAssign,Sub,SubAssign,Mul,Div};

#[derive(Clone,Copy,Debug,Default,PartialEq)]
pub struct V3 {pub x:f64,pub y:f64,pub z:f64}
impl V3 {
    pub const fn new(x:f64,y:f64,z:f64)->Self {Self{x,y,z}}
    pub fn dot(self,other:Self)->f64 {
        self.x*other.x+self.y*other.y+self.z*other.z
    }
    pub fn len_sq(self)->f64 {self.dot(self)}
    pub fn len(self)->f64 {self.len_sq().sqrt()}
    pub fn normalized(self)->Self {
        let r=self.len(); if r>0.0 {self/r} else {Self::default()}
    }
    pub fn finite(self)->bool {
        self.x.is_finite()&&self.y.is_finite()&&self.z.is_finite()
    }
}
impl Add for V3 {
    type Output=Self;
    fn add(self,o:Self)->Self{Self::new(self.x+o.x,self.y+o.y,self.z+o.z)}
}
impl AddAssign for V3 {fn add_assign(&mut self,o:Self){*self=*self+o;}}
impl Sub for V3 {
    type Output=Self;
    fn sub(self,o:Self)->Self{Self::new(self.x-o.x,self.y-o.y,self.z-o.z)}
}
impl SubAssign for V3 {fn sub_assign(&mut self,o:Self){*self=*self-o;}}
impl Mul<f64> for V3 {
    type Output=Self;
    fn mul(self,s:f64)->Self{Self::new(self.x*s,self.y*s,self.z*s)}
}
impl Div<f64> for V3 {
    type Output=Self;
    fn div(self,s:f64)->Self{Self::new(self.x/s,self.y/s,self.z/s)}
}

#[derive(Clone,Copy,Debug)]
pub struct Fluid3 {
    pub pos_m:V3,
    pub vel_m_s:V3,
    pub mass_kg:f64,
}
#[derive(Clone,Copy,Debug)]
pub struct Massive3 {
    pub pos_m:V3,
    pub vel_m_s:V3,
    pub mass_kg:f64,
    pub radius_m:f64,
}
#[derive(Clone,Copy,Debug)]
pub enum Sph3Eos {
    Polytrope,
    /// Fully relativistic zero-temperature Fermi electron EOS, mu_e=A/Z.
    ColdDegenerate {mu_e:f64},
}
#[derive(Clone,Copy,Debug)]
pub struct Sph3Config {
    pub smoothing_m:f64,
    pub gravity_softening_m:f64,
    pub polytropic_k:f64,
    pub gamma:f64,
    pub eos:Sph3Eos,
}
impl Sph3Config {
    pub fn white_dwarf_proxy()->Self {
        Self {
            smoothing_m:3.0e6,
            gravity_softening_m:1.2e6,
            polytropic_k:3.0e8,
            gamma:5.0/3.0,
            eos:Sph3Eos::Polytrope,
        }
    }
    pub fn cold_white_dwarf(mu_e:f64)->Option<Self> {
        DegenerateElectrons::new(mu_e)?;
        Some(Self {smoothing_m:1.5e6,gravity_softening_m:0.8e6,
            polytropic_k:0.0,gamma:5.0/3.0,
            eos:Sph3Eos::ColdDegenerate{mu_e}})
    }
    pub fn valid(&self)->bool {
        self.smoothing_m.is_finite()&&self.smoothing_m>0.0&&
        self.gravity_softening_m.is_finite()&&self.gravity_softening_m>0.0&&
        self.polytropic_k.is_finite()&&self.polytropic_k>=0.0&&
        self.gamma.is_finite()&&self.gamma>1.0&&
        match self.eos {
            Sph3Eos::Polytrope=>true,
            Sph3Eos::ColdDegenerate{mu_e}=>DegenerateElectrons::new(mu_e).is_some(),
        }
    }
}
/// Normalized 3D cubic spline kernel W: m^-3.
pub fn cubic_spline_3d(r:f64,h:f64)->f64 {
    assert!(h.is_finite()&&h>0.0);
    let q=r/h;
    let norm=1.0/(PI*h*h*h);
    if q<0.0||!q.is_finite() {0.0}
    else if q<1.0 {norm*(1.0-1.5*q*q+0.75*q*q*q)}
    else if q<2.0 {norm*0.25*(2.0-q).powi(3)}
    else {0.0}
}
/// Radial derivative dW/dr: m^-4.
pub fn spline_radial_derivative_3d(r:f64,h:f64)->f64 {
    assert!(h.is_finite()&&h>0.0);
    let q=r/h;
    let norm=1.0/(PI*h.powi(4));
    if q<0.0||!q.is_finite() {0.0}
    else if q<1.0 {norm*(-3.0*q+2.25*q*q)}
    else if q<2.0 {-norm*0.75*(2.0-q).powi(2)}
    else {0.0}
}
pub fn densities_kg_m3(particles:&[Fluid3],h:f64)->Vec<f64> {
    let mut density:Vec<f64>=particles.iter()
        .map(|p|p.mass_kg*cubic_spline_3d(0.0,h)).collect();
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let w=cubic_spline_3d((particles[i].pos_m-
                particles[j].pos_m).len(),h);
            density[i]+=particles[j].mass_kg*w;
            density[j]+=particles[i].mass_kg*w;
        }
    }
    density
}
pub fn pressure_pa(rho:f64,cfg:Sph3Config)->f64 {
    match cfg.eos {
        Sph3Eos::Polytrope=>cfg.polytropic_k*rho.powf(cfg.gamma),
        Sph3Eos::ColdDegenerate{mu_e}=>
            DegenerateElectrons::new(mu_e).unwrap().pressure_pa(rho),
    }
}
pub fn internal_energy_j_per_kg(rho:f64,cfg:Sph3Config)->f64 {
    match cfg.eos {
        Sph3Eos::Polytrope=>
            cfg.polytropic_k*rho.powf(cfg.gamma-1.0)/(cfg.gamma-1.0),
        Sph3Eos::ColdDegenerate{mu_e}=>
            DegenerateElectrons::new(mu_e).unwrap().internal_energy_j_kg(rho),
    }
}
#[derive(Debug)]
pub struct Sph3Forces {
    pub central_accel_m_s2:V3,
    pub parcel_accel_m_s2:Vec<V3>,
    pub densities_kg_m3:Vec<f64>,
}
pub fn accelerations(
    central:&Massive3,particles:&[Fluid3],cfg:Sph3Config,
)->Sph3Forces {
    assert!(cfg.valid());
    let n=particles.len();
    let rho=densities_kg_m3(particles,cfg.smoothing_m);
    let pressure_over_rho2:Vec<f64>=rho.iter()
        .map(|density|pressure_pa(*density,cfg)/(density*density)).collect();
    let mut a=vec![V3::default();n];
    let mut central_a=V3::default();
    let eps2=cfg.gravity_softening_m*cfg.gravity_softening_m;
    for (i,p) in particles.iter().enumerate() {
        let delta=p.pos_m-central.pos_m;
        let r2=delta.len_sq()+eps2;
        let field=delta*(G_SI/(r2*r2.sqrt()));
        a[i]-=field*central.mass_kg;
        central_a+=field*p.mass_kg;
    }
    for i in 0..n {
        for j in (i+1)..n {
            let d=particles[j].pos_m-particles[i].pos_m;
            let r2=d.len_sq()+eps2;
            let grav=d*(G_SI/(r2*r2.sqrt()));
            a[i]+=grav*particles[j].mass_kg;
            a[j]-=grav*particles[i].mass_kg;
            let r=d.len();
            if r>0.0 && r<2.0*cfg.smoothing_m {
                let dw=spline_radial_derivative_3d(r,cfg.smoothing_m);
                let pressure_factor=
                    -(pressure_over_rho2[i]+pressure_over_rho2[j])*dw;
                let push=(particles[i].pos_m-particles[j].pos_m)*
                    (pressure_factor/r);
                a[i]+=push*particles[j].mass_kg;
                a[j]-=push*particles[i].mass_kg;
            }
        }
    }
    Sph3Forces{central_accel_m_s2:central_a,
        parcel_accel_m_s2:a,densities_kg_m3:rho}
}
pub fn total_mass_kg(central:&Massive3,particles:&[Fluid3])->f64 {
    central.mass_kg+particles.iter().map(|p|p.mass_kg).sum::<f64>()
}
pub fn total_momentum(central:&Massive3,particles:&[Fluid3])->V3 {
    particles.iter().fold(central.vel_m_s*central.mass_kg,
        |p,x|p+x.vel_m_s*x.mass_kg)
}
pub fn center_of_mass(central:&Massive3,particles:&[Fluid3])->V3 {
    let total=total_mass_kg(central,particles);
    (particles.iter().fold(central.pos_m*central.mass_kg,
        |p,x|p+x.pos_m*x.mass_kg))/total
}
/// Pair-symmetric Verlet. Energy drift is a timestep convergence diagnostic.
pub fn advance(
    central:&mut Massive3,particles:&mut [Fluid3],
    cfg:Sph3Config,dt_s:f64,
) {
    assert!(cfg.valid()&&dt_s.is_finite()&&dt_s>=0.0);
    if dt_s==0.0 {return;}
    let initial=accelerations(central,particles,cfg);
    central.pos_m+=central.vel_m_s*dt_s+
        initial.central_accel_m_s2*(0.5*dt_s*dt_s);
    central.vel_m_s+=initial.central_accel_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(initial.parcel_accel_m_s2) {
        p.pos_m+=p.vel_m_s*dt_s+a*(0.5*dt_s*dt_s);
        p.vel_m_s+=a*(0.5*dt_s);
    }
    let final_acc=accelerations(central,particles,cfg);
    central.vel_m_s+=final_acc.central_accel_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(final_acc.parcel_accel_m_s2) {
        p.vel_m_s+=a*(0.5*dt_s);
    }
}
/// Newtonian total energy with polytropic internal energy.
pub fn energy_j(central:&Massive3,particles:&[Fluid3],cfg:Sph3Config)->f64 {
    let rho=densities_kg_m3(particles,cfg.smoothing_m);
    let mut e=0.5*central.mass_kg*central.vel_m_s.len_sq();
    for (i,p) in particles.iter().enumerate() {
        e+=p.mass_kg*(0.5*p.vel_m_s.len_sq()+
            internal_energy_j_per_kg(rho[i],cfg));
        let d=(central.pos_m-p.pos_m).len_sq()+
            cfg.gravity_softening_m.powi(2);
        e-=G_SI*central.mass_kg*p.mass_kg/d.sqrt();
    }
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let d=(particles[i].pos_m-particles[j].pos_m).len_sq()+
                cfg.gravity_softening_m.powi(2);
            e-=G_SI*particles[i].mass_kg*particles[j].mass_kg/d.sqrt();
        }
    }
    e
}
/// Seed a 3D star with symmetric pairs of samples (Fibonacci directions).
/// Preserves the progenitor's center of mass and net momentum at initialization.
pub fn seed_from_progenitor(
    source:SiBody,center_z_m:f64,velocity_z_m_s:f64,count:usize,
)->(Vec<Fluid3>,Massive3) {
    assert!(count>=4 && count%2==0);
    let mut out=Vec::with_capacity(count);
    let half=count/2;
    let golden=PI*(3.0-5.0f64.sqrt());
    let pos0=V3::new(source.pos_m.x,source.pos_m.y,center_z_m);
    let vel0=V3::new(source.vel_m_s.x,source.vel_m_s.y,velocity_z_m_s);
    for i in 0..half {
        let z=1.0-2.0*(i as f64+0.5)/half as f64;
        let phi=i as f64*golden;
        let radial_xy=(1.0-z*z).sqrt();
        let dir=V3::new(radial_xy*phi.cos(),radial_xy*phi.sin(),z);
        let fraction=((i as f64+0.5)/half as f64).cbrt();
        let offset=dir*source.radius_m*fraction;
        let kick=dir*22_000.0;
        let mass=source.mass_kg/count as f64;
        out.push(Fluid3 {pos_m:pos0+offset,
            vel_m_s:vel0+kick,mass_kg:mass});
        out.push(Fluid3 {pos_m:pos0-offset,
            vel_m_s:vel0-kick,mass_kg:mass});
    }
    // The returned synthetic placeholder preserves metadata for source
    // comparison; caller supplies the actual central black hole separately.
    let placeholder=Massive3 {
        pos_m:pos0,vel_m_s:vel0,
        mass_kg:source.mass_kg,radius_m:source.radius_m,
    };
    (out,placeholder)
}
/// Equal-mass, antipodal samples with radial locations from the hydrostatic
/// cumulative mass profile. This gives a stratified progenitor distribution.
/// SPH particle sampling remains an approximation to continuum equilibrium.
pub fn seed_from_hydrostatic(
    profile:&WhiteDwarfProfile,center_m:V3,velocity_m_s:V3,
    count:usize,
)->Vec<Fluid3> {
    assert!(count>=4 && count%2==0);
    let half=count/2;
    let mut out=Vec::with_capacity(count);
    let golden=PI*(3.0-5.0f64.sqrt());
    let mass_each=profile.mass_kg/count as f64;
    for i in 0..half {
        let z=1.0-2.0*(i as f64+0.5)/half as f64;
        let angle=golden*i as f64;
        let xy=(1.0-z*z).sqrt();
        let dir=V3::new(xy*angle.cos(),xy*angle.sin(),z);
        let radius=profile.radius_at_mass_fraction((i as f64+0.5)/half as f64);
        let displacement=dir*radius;
        out.push(Fluid3{pos_m:center_m+displacement,
            vel_m_s:velocity_m_s,mass_kg:mass_each});
        out.push(Fluid3{pos_m:center_m-displacement,
            vel_m_s:velocity_m_s,mass_kg:mass_each});
    }
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    fn setup()->(Massive3,Vec<Fluid3>,Sph3Config) {
        let star=SiBody::new(0.65*SOLAR_MASS_KG,9.0e6,
            crate::V2::new(1.0e8,0.0),
            crate::V2::new(-3.0e5,2.0e4));
        let (p,_)=seed_from_progenitor(star,0.0,0.0,48);
        (Massive3{pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:10.0*SOLAR_MASS_KG,radius_m:3.0e4},p,
            Sph3Config::white_dwarf_proxy())
    }
    #[test]
    fn cubic_spline_integrates_to_unit_volume() {
        let h=3.0e6;
        let mut integral=0.0;
        for i in 0..10_000 {
            let r=(i as f64+0.5)*2.0*h/10_000.0;
            integral+=4.0*PI*r*r*cubic_spline_3d(r,h)*2.0*h/10_000.0;
        }
        assert!((integral-1.0).abs()<1e-6);
    }
    #[test]
    fn derivative_is_negative_and_compact() {
        let h=1e6;
        assert!(spline_radial_derivative_3d(h,h)<0.0);
        assert_eq!(spline_radial_derivative_3d(2.0*h,h),0.0);
    }
    #[test]
    fn initial_symmetric_pairs_preserve_mass_and_center() {
        let (central,p,_)=setup();
        let expected=central.mass_kg+0.65*SOLAR_MASS_KG;
        let sum=total_mass_kg(&central,&p);
        assert!((sum-expected).abs()/expected<1e-13);
        let m: f64=p.iter().map(|x|x.mass_kg).sum();
        let cm=p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/m;
        let v=p.iter().fold(V3::default(),|s,x|s+x.vel_m_s*x.mass_kg)/m;
        assert!((cm-V3::new(1.0e8,0.0,0.0)).len()<1e-5);
        assert!((v-V3::new(-3.0e5,2.0e4,0.0)).len()<1e-6);
    }
    #[test]
    fn total_pair_force_matches_newtons_third_law() {
        let (bh,p,cfg)=setup();
        let accelerations=accelerations(&bh,&p,cfg);
        let net=p.iter().zip(accelerations.parcel_accel_m_s2)
            .fold(accelerations.central_accel_m_s2*bh.mass_kg,
                |s,(p,a)|s+a*p.mass_kg);
        let denom=p[0].mass_kg*200_000.0;
        assert!(net.len()/denom<1e-9,"net force ratio={}",net.len()/denom);
    }
    #[test]
    fn three_dimensional_motion_departs_from_equatorial_plane() {
        let (mut bh,mut p,cfg)=setup();
        assert!(p.iter().any(|x|x.pos_m.z.abs()>1.0e5));
        let initial_z=p[2].pos_m.z;
        for _ in 0..20 {advance(&mut bh,&mut p,cfg,0.01);}
        assert!((p[2].pos_m.z-initial_z).abs()>0.01);
    }
    #[test]
    fn short_run_preserves_linear_momentum_and_mass() {
        let (mut bh,mut p,cfg)=setup();
        let initial=total_momentum(&bh,&p);
        let m=total_mass_kg(&bh,&p);
        for _ in 0..100 {advance(&mut bh,&mut p,cfg,0.005);}
        let drift=(total_momentum(&bh,&p)-initial).len()/
            (m*3e5);
        assert!(drift<1e-10,"momentum drift ratio {drift}");
        assert!((total_mass_kg(&bh,&p)-m).abs()/m<1e-13);
    }
    #[test]
    fn total_energy_short_interval_stays_near_initial() {
        let (mut bh,mut p,cfg)=setup();
        let e=energy_j(&bh,&p,cfg);
        for _ in 0..40 {advance(&mut bh,&mut p,cfg,0.002);}
        let drift=(energy_j(&bh,&p,cfg)-e).abs()/e.abs();
        assert!(drift<0.001,"energy drift ratio={drift}");
    }
    #[test]
    fn cold_degenerate_eos_is_wired_to_exact_reference() {
        let cfg=Sph3Config::cold_white_dwarf(2.0).unwrap();
        let e=DegenerateElectrons::new(2.0).unwrap();
        for rho in [1.0e7,1.0e9,1.0e11] {
            let a=pressure_pa(rho,cfg);
            assert!((a-e.pressure_pa(rho)).abs()/a<1e-13);
            let u=internal_energy_j_per_kg(rho,cfg);
            assert!((u-e.internal_energy_j_kg(rho)).abs()/u<1e-13);
        }
        assert!(!Sph3Config::cold_white_dwarf(0.0).is_some());
    }
    #[test]
    fn hydrostatic_seed_preserves_bulk_mass_center_and_momentum() {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=crate::white_dwarf::profile_for_mass(
            0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let center=V3::new(1e8,-2e7,3e7);
        let velocity=V3::new(1e5,2e4,-3e4);
        let samples=seed_from_hydrostatic(&star,center,velocity,128);
        let mass: f64=samples.iter().map(|p|p.mass_kg).sum();
        let com=samples.iter().fold(V3::default(),
            |s,p|s+p.pos_m*p.mass_kg)/mass;
        let momentum=samples.iter().fold(V3::default(),
            |s,p|s+p.vel_m_s*p.mass_kg)/mass;
        assert!((mass-star.mass_kg).abs()/mass<1e-12);
        assert!((com-center).len()<1e-5);
        assert!((momentum-velocity).len()<1e-6);
        assert!(samples.iter().all(|p|(p.pos_m-center).len()<=star.radius_m));
    }
    #[test]
    fn cold_eos_pair_forces_conserve_total_linear_momentum() {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=crate::white_dwarf::profile_for_mass(
            0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let mut samples=seed_from_hydrostatic(&star,
            V3::new(1e8,0.0,0.0),V3::new(-3e5,0.0,0.0),64);
        let mut central=Massive3{pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:10.0*SOLAR_MASS_KG,radius_m:3.0e4};
        let cfg=Sph3Config::cold_white_dwarf(2.0).unwrap();
        let p0=total_momentum(&central,&samples);
        for _ in 0..8 {advance(&mut central,&mut samples,cfg,0.001);}
        let p1=total_momentum(&central,&samples);
        assert!((p1-p0).len()/(star.mass_kg*3.0e5)<1e-10);
        assert!(samples.iter().all(|x|x.vel_m_s.finite()));
    }
    #[test]
    fn zero_step_identity() {
        let (mut bh,mut p,cfg)=setup();
        let original=bh.pos_m;
        let v=p[0].vel_m_s;
        advance(&mut bh,&mut p,cfg,0.0);
        assert_eq!(bh.pos_m,original);
        assert_eq!(p[0].vel_m_s,v);
    }
}
