//! Variable-h 3D smoothed-particle hydro for cold-degenerate stellar tests.
//!
//! The SPH pressure force includes a grad-h/Omega correction derived from
//! h_i = eta*(m_i/rho_i)^(1/3). Density and smoothing length are iterated
//! together. Gravity remains *fixed Plummer-softened* with a separate
//! epsilon by default. The optional symmetric variable-Plummer model
//! derives its gravitational grad-h (zeta) force from the same pair
//! potential used by the energy diagnostic. This is variational adaptive
//! Plummer softening, not a Poisson kernel matched to the cubic spline.
//!
//! The central massive body, where present, recoils from each tracer.
//! Pairwise pressure and gravitational forces conserve linear momentum
//! to rounding precision. Computational complexity is O(N²) per force pass.
use crate::sph3::{
    self,Fluid3,Massive3,Sph3Config,Sph3Eos,V3,
};
use std::f64::consts::PI;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum GravityModel {
    FixedPlummer,
    /// Pair potential -G m_i m_j /2 * [1/sqrt(r²+h_i²)+1/sqrt(r²+h_j²)].
    /// Requires the corresponding zeta/Omega grad-h force.
    MatchedAdaptivePlummer,
}
#[derive(Clone,Copy,Debug)]
pub struct AdaptiveConfig {
    pub eta:f64,
    pub min_h_m:f64,
    pub max_h_m:f64,
    pub gravity_softening_m:f64,
    pub gravity_model:GravityModel,
    pub density_iterations:usize,
}
impl AdaptiveConfig {
    pub fn from_star_radius(radius_m:f64)->Self {
        assert!(radius_m.is_finite()&&radius_m>0.0);
        Self{
            eta:1.2,min_h_m:0.08*radius_m,
            max_h_m:0.9*radius_m,
            gravity_softening_m:0.10*radius_m,
            gravity_model:GravityModel::FixedPlummer,
            density_iterations:9,
        }
    }
    pub fn with_adaptive_gravity(mut self)->Self {
        self.gravity_model=GravityModel::MatchedAdaptivePlummer;
        // Accurate coupled h/rho closure is particularly important for
        // variational gravitational-force consistency.
        self.density_iterations=24;
        self
    }
    pub fn valid(self)->bool {
        self.eta.is_finite()&&self.eta>0.5&&self.eta<3.0&&
        self.min_h_m.is_finite()&&self.min_h_m>0.0&&
        self.max_h_m.is_finite()&&self.max_h_m>self.min_h_m&&
        self.gravity_softening_m.is_finite()&&self.gravity_softening_m>0.0&&
        self.density_iterations>=2&&self.density_iterations<=32
    }
}

#[derive(Clone,Debug)]
pub struct AdaptiveState {
    pub smoothing_m:Vec<f64>,
    pub densities_kg_m3:Vec<f64>,
    /// grad-h factors Ω_i=1+(h_i/(3ρ_i)) ∂ρ_i/∂h_i.
    pub omega:Vec<f64>,
}
impl AdaptiveState {
    pub fn new(particles:&[Fluid3],initial_radius_m:f64,cfg:AdaptiveConfig)->Self {
        assert!(cfg.valid()&&initial_radius_m>0.0&&particles.len()>=2);
        let nominal_spacing=(4.0*PI*initial_radius_m.powi(3)/
            (3.0*particles.len() as f64)).cbrt();
        let initial_h=(cfg.eta*nominal_spacing)
            .clamp(cfg.min_h_m,cfg.max_h_m);
        Self{
            smoothing_m:vec![initial_h;particles.len()],
            densities_kg_m3:vec![0.0;particles.len()],
            omega:vec![1.0;particles.len()],
        }
    }
    pub fn valid_for(&self,n:usize)->bool {
        self.smoothing_m.len()==n&&
        self.densities_kg_m3.len()==n&&self.omega.len()==n
    }
}

fn kernel_and_derivative_h(distance_m:f64,h:f64)->(f64,f64) {
    let q=distance_m/h;
    if q>=2.0 {return (0.0,0.0);}
    let (f,dfdq)=if q<1.0 {
        (1.0-1.5*q*q+0.75*q.powi(3),-3.0*q+2.25*q*q)
    } else {
        (0.25*(2.0-q).powi(3),-0.75*(2.0-q).powi(2))
    };
    let norm=1.0/(PI*h.powi(3));
    let w=norm*f;
    let dwdh=-norm/h*(3.0*f+q*dfdq);
    (w,dwdh)
}
fn density_of_particle(i:usize,particles:&[Fluid3],h:f64)->(f64,f64) {
    let mut rho=0.0;
    let mut derivative=0.0;
    for x in particles {
        let (w,dwdh)=kernel_and_derivative_h(
            (particles[i].pos_m-x.pos_m).len(),h);
        rho+=x.mass_kg*w;
        derivative+=x.mass_kg*dwdh;
    }
    (rho,derivative)
}

/// Compute ρ_i and h_i from their coupled relation by fixed-point
/// iteration, with an underrelaxed h update. The explicit caps on h
/// are recorded, and Ω=1 at a hard cap where dh/dρ must vanish.
pub fn update_adaptive_densities(
    particles:&[Fluid3],state:&mut AdaptiveState,cfg:AdaptiveConfig,
) {
    assert!(cfg.valid()&&state.valid_for(particles.len()));
    for _ in 0..cfg.density_iterations {
        for i in 0..particles.len() {
            let h=state.smoothing_m[i];
            let (rho,_)=density_of_particle(i,particles,h);
            assert!(rho>0.0&&rho.is_finite());
            let target=cfg.eta*
                (particles[i].mass_kg/rho).cbrt();
            state.smoothing_m[i]=
                (0.45*h+0.55*target).clamp(cfg.min_h_m,cfg.max_h_m);
        }
    }
    for i in 0..particles.len() {
        let h=state.smoothing_m[i];
        let (rho,drho_dh)=density_of_particle(i,particles,h);
        state.densities_kg_m3[i]=rho;
        let at_cap=(h-cfg.min_h_m).abs()<cfg.min_h_m*1e-8||
            (h-cfg.max_h_m).abs()<cfg.max_h_m*1e-8;
        let omega=if at_cap {1.0}
            else {1.0+h*drho_dh/(3.0*rho)};
        // The floor applies only when the closure becomes numerically
        // singular; record min Ω in benchmark diagnostics.
        state.omega[i]=omega.max(0.02);
    }
}


/// Plummer potential per unit Gm, defined with a negative sign.
/// phi(r,h)=-1/sqrt(r²+h²).
pub fn adaptive_plummer_potential(r:f64,h:f64)->f64 {
    assert!(r>=0.0 && r.is_finite() && h>0.0 && h.is_finite());
    -1.0/(r*r+h*h).sqrt()
}
/// Derivative of the negative Plummer potential with respect to h.
/// dphi/dh = h/(r²+h²)^(3/2).
pub fn adaptive_plummer_dphi_dh(r:f64,h:f64)->f64 {
    let sq=r*r+h*h;
    h/(sq*sq.sqrt())
}
/// zeta_i= (dh_i/drho_i) sum_{j!=i} m_j dphi_ij(h_i)/dh_i.
/// With h∝rho^(-1/3), dh/drho=-h/(3rho) away from a hard h cap.
/// This is a negative coefficient for the Plummer potential convention.
pub fn adaptive_gravity_zeta(
    particles:&[Fluid3],state:&AdaptiveState,cfg:AdaptiveConfig,
)->Vec<f64> {
    assert!(state.valid_for(particles.len()));
    let mut zeta=vec![0.0;particles.len()];
    for (i,p) in particles.iter().enumerate() {
        let h=state.smoothing_m[i];
        let capped=(h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8 ||
            (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8;
        if capped {continue;}
        let mut derivative_sum=0.0;
        for (j,q) in particles.iter().enumerate() {
            if i==j {continue;}
            derivative_sum+=q.mass_kg*
                adaptive_plummer_dphi_dh((p.pos_m-q.pos_m).len(),h);
        }
        zeta[i]=-h/(3.0*state.densities_kg_m3[i])*derivative_sum;
    }
    zeta
}

#[derive(Debug)]
pub struct AdaptiveForces {
    pub central_acceleration_m_s2:V3,
    pub accelerations_m_s2:Vec<V3>,
    pub minimum_omega:f64,
    pub minimum_h_m:f64,
    pub maximum_h_m:f64,
}
pub fn forces(
    central:&Massive3,particles:&[Fluid3],
    state:&mut AdaptiveState,physics:Sph3Config,cfg:AdaptiveConfig,
)->AdaptiveForces {
    assert!(cfg.valid()&&physics.valid()&&state.valid_for(particles.len()));
    update_adaptive_densities(particles,state,cfg);
    let n=particles.len();
    let mut acceleration=vec![V3::default();n];
    let mut central_acc=V3::default();
    let inv_eps2=cfg.gravity_softening_m.powi(2);
    let zeta=if cfg.gravity_model==GravityModel::MatchedAdaptivePlummer {
        adaptive_gravity_zeta(particles,state,cfg)
    } else {vec![0.0;particles.len()]};
    let pressure_factor=state.densities_kg_m3.iter()
        .enumerate().map(|(i,rho)|{
            let pressure=sph3::pressure_pa(*rho,physics);
            pressure/(state.omega[i]*rho*rho)
        }).collect::<Vec<_>>();

    for (i,p) in particles.iter().enumerate() {
        let d=p.pos_m-central.pos_m;
        let r2=d.len_sq()+inv_eps2;
        let g=d*(crate::physics::G_SI/(r2*r2.sqrt()));
        acceleration[i]-=g*central.mass_kg;
        central_acc+=g*p.mass_kg;
    }
    for i in 0..n {
        for j in i+1..n {
            let d=particles[j].pos_m-particles[i].pos_m;
            let r=d.len();
            let r2=d.len_sq();
            let grav_scalar=match cfg.gravity_model {
                GravityModel::FixedPlummer=>
                    1.0/((r2+inv_eps2)*(r2+inv_eps2).sqrt()),
                GravityModel::MatchedAdaptivePlummer=>{
                    let hi2=state.smoothing_m[i].powi(2);
                    let hj2=state.smoothing_m[j].powi(2);
                    0.5*((r2+hi2).powf(-1.5)+(r2+hj2).powf(-1.5))
                }
            };
            let grav=d*(crate::physics::G_SI*grav_scalar);
            acceleration[i]+=grav*particles[j].mass_kg;
            acceleration[j]-=grav*particles[i].mass_kg;
            if r==0.0 {continue;}
            let grad_i=sph3::spline_radial_derivative_3d(
                r,state.smoothing_m[i]);
            let grad_j=sph3::spline_radial_derivative_3d(
                r,state.smoothing_m[j]);
            let coeff=pressure_factor[i]*grad_i+
                pressure_factor[j]*grad_j;
            // ∇_i W_ij=(x_i-x_j)/r * dW/dr; symmetric pair.
            let push=(particles[i].pos_m-particles[j].pos_m)*
                (-coeff/r);
            acceleration[i]+=push*particles[j].mass_kg;
            acceleration[j]-=push*particles[i].mass_kg;
            if cfg.gravity_model==GravityModel::MatchedAdaptivePlummer {
                // Price-Monaghan zeta correction derived by differentiating
                // the pair gravitational potential through h_i(rho_i).
                let gcoeff=0.5*crate::physics::G_SI*
                    (zeta[i]/state.omega[i]*grad_i+
                     zeta[j]/state.omega[j]*grad_j)/r;
                let correction=d*gcoeff;
                acceleration[i]+=correction*particles[j].mass_kg;
                acceleration[j]-=correction*particles[i].mass_kg;
            }
        }
    }
    AdaptiveForces {
        central_acceleration_m_s2:central_acc,
        accelerations_m_s2:acceleration,
        minimum_omega:state.omega.iter().copied().fold(f64::INFINITY,f64::min),
        minimum_h_m:state.smoothing_m.iter().copied().fold(f64::INFINITY,f64::min),
        maximum_h_m:state.smoothing_m.iter().copied().fold(0.0_f64,f64::max),
    }
}

/// Time-symmetric velocity Verlet with variable-h pressure and fixed
/// gravitational softening. Momentum is conserved by pairwise forces.
pub fn advance(
    central:&mut Massive3,particles:&mut [Fluid3],
    state:&mut AdaptiveState,physics:Sph3Config,
    cfg:AdaptiveConfig,dt_s:f64,
) {
    assert!(dt_s.is_finite()&&dt_s>=0.0);
    if dt_s==0.0 {return;}
    let before=forces(central,particles,state,physics,cfg);
    central.pos_m+=central.vel_m_s*dt_s+
        before.central_acceleration_m_s2*(0.5*dt_s*dt_s);
    central.vel_m_s+=before.central_acceleration_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(before.accelerations_m_s2) {
        p.pos_m+=p.vel_m_s*dt_s+a*(0.5*dt_s*dt_s);
        p.vel_m_s+=a*(0.5*dt_s);
    }
    let after=forces(central,particles,state,physics,cfg);
    central.vel_m_s+=after.central_acceleration_m_s2*(0.5*dt_s);
    for (p,a) in particles.iter_mut().zip(after.accelerations_m_s2) {
        p.vel_m_s+=a*(0.5*dt_s);
    }
}

pub fn isolated_force_rms_m_s2(
    particles:&[Fluid3],state:&mut AdaptiveState,
    physics:Sph3Config,cfg:AdaptiveConfig,
)->f64 {
    let empty=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let f=forces(&empty,particles,state,physics,cfg);
    let mass=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    (particles.iter().zip(f.accelerations_m_s2)
        .map(|(p,a)|p.mass_kg*a.len_sq()).sum::<f64>()/mass).sqrt()
}
pub fn adaptive_energy_j(
    central:&Massive3,particles:&[Fluid3],
    state:&mut AdaptiveState,
    physics:Sph3Config,cfg:AdaptiveConfig,
)->f64 {
    update_adaptive_densities(particles,state,cfg);
    energy_j_cached(central,particles,state,physics,cfg)
}
/// Energy with the smoothing-length/density snapshot already updated by the
/// most recent force call, for low-overhead per-frame telemetry.
pub fn energy_j_cached(
    central:&Massive3,particles:&[Fluid3],
    state:&AdaptiveState,
    physics:Sph3Config,cfg:AdaptiveConfig,
)->f64 {
    assert!(state.valid_for(particles.len()));
    let mut e=central.mass_kg*central.vel_m_s.len_sq()*0.5;
    for (i,p) in particles.iter().enumerate() {
        e+=p.mass_kg*(0.5*p.vel_m_s.len_sq()+
            sph3::internal_energy_j_per_kg(state.densities_kg_m3[i],physics));
        let r2=(p.pos_m-central.pos_m).len_sq()+
            cfg.gravity_softening_m.powi(2);
        e-=crate::physics::G_SI*p.mass_kg*central.mass_kg/r2.sqrt();
    }
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let r2=(particles[i].pos_m-particles[j].pos_m).len_sq();
            let pair_inverse_r=match cfg.gravity_model {
                GravityModel::FixedPlummer=>
                    1.0/(r2+cfg.gravity_softening_m.powi(2)).sqrt(),
                GravityModel::MatchedAdaptivePlummer=>
                    0.5*(1.0/(r2+state.smoothing_m[i].powi(2)).sqrt()+
                         1.0/(r2+state.smoothing_m[j].powi(2)).sqrt()),
            };
            e-=crate::physics::G_SI*particles[i].mass_kg*
                particles[j].mass_kg*pair_inverse_r;
        }
    }
    e
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn setup(n:usize)->(Vec<Fluid3>,Sph3Config,AdaptiveConfig,AdaptiveState) {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let samples=sph3::seed_from_hydrostatic(
            &star,V3::new(5e7,0.0,0.0),V3::default(),n);
        let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m);
        let state=AdaptiveState::new(&samples,star.radius_m,cfg);
        (samples,physics,cfg,state)
    }
    #[test]
    fn kernel_h_derivative_matches_finite_difference() {
        for q in [0.0,0.15,0.6,1.15,1.9] {
            let h=1.7e6;
            let r=q*h;
            let (_,dwdh)=kernel_and_derivative_h(r,h);
            let dh=1.0;
            let finite=(sph3::cubic_spline_3d(r,h+dh)-
                sph3::cubic_spline_3d(r,h-dh))/(2.0*dh);
            assert!((dwdh-finite).abs()<=
                3e-6*dwdh.abs().max(1e-30));
        }
    }
    #[test]
    fn density_smoothing_and_omega_stay_finite() {
        let (particles,phys,cfg,mut state)=setup(48);
        update_adaptive_densities(&particles,&mut state,cfg);
        assert!(state.smoothing_m.iter().all(|h|
            *h>=cfg.min_h_m&&*h<=cfg.max_h_m));
        assert!(state.densities_kg_m3.iter().all(|rho|
            rho.is_finite()&&*rho>0.0));
        assert!(state.omega.iter().all(|o|o.is_finite()&&*o>=0.02));
        let f=forces(&Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0
        },&particles,&mut state,phys,cfg);
        assert!(f.minimum_omega>0.0);
    }
    #[test]
    fn symmetric_adaptive_forces_preserve_total_momentum() {
        let (particles,phys,cfg,mut state)=setup(64);
        let hole=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:10.0*SOLAR_MASS_KG,radius_m:3e4
        };
        let force=forces(&hole,&particles,&mut state,phys,cfg);
        let net=particles.iter().zip(force.accelerations_m_s2)
            .fold(force.central_acceleration_m_s2*hole.mass_kg,
                |sum,(p,a)|sum+a*p.mass_kg);
        let scale=particles[0].mass_kg*2.0e6;
        assert!(net.len()/scale<1e-8,"net {}",net.len()/scale);
    }
    #[test]
    fn plummer_h_derivative_matches_finite_difference() {
        let r=3.0e6;let h=2.0e6;let dh=1.0;
        let finite=(adaptive_plummer_potential(r,h+dh)-
            adaptive_plummer_potential(r,h-dh))/(2.0*dh);
        let analytic=adaptive_plummer_dphi_dh(r,h);
        assert!((analytic-finite).abs()/analytic<1e-7);
    }
    #[test]
    fn matched_softening_strengthens_gravity_via_grad_h_and_conserves_momentum() {
        let (particles,phys,old,_) =setup(64);
        let cfg=old.with_adaptive_gravity();
        let mut state=AdaptiveState::new(&particles,8.4e6,cfg);
        let hole=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let f=forces(&hole,&particles,&mut state,phys,cfg);
        let net=particles.iter().zip(&f.accelerations_m_s2)
            .fold(V3::default(),|sum,(p,a)|sum+*a*p.mass_kg);
        let scale=particles[0].mass_kg*f.accelerations_m_s2[0].len();
        assert!(net.len()/scale<1e-9);
        let zeta=adaptive_gravity_zeta(&particles,&state,cfg);
        assert!(zeta.iter().all(|z|z.is_finite()&&*z<=0.0));
        assert!(zeta.iter().any(|z|*z<0.0));
    }
    #[test]
    fn matched_pair_potential_and_force_drift_remain_finite() {
        let (mut p,phys,old,_) =setup(48);
        let cfg=old.with_adaptive_gravity();
        let mut state=AdaptiveState::new(&p,8.4e6,cfg);
        let mut hole=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let e0=adaptive_energy_j(&hole,&p,&mut state,phys,cfg);
        let p0=sph3::total_momentum(&hole,&p);
        for _ in 0..60 {
            advance(&mut hole,&mut p,&mut state,phys,cfg,0.002);
        }
        let e1=adaptive_energy_j(&hole,&p,&mut state,phys,cfg);
        let momentum=(sph3::total_momentum(&hole,&p)-p0).len();
        assert!((e1-e0).abs()/e0.abs()<1e-3,
            "matched gravitational energy drift = {}",(e1-e0).abs()/e0.abs());
        assert!(momentum/(p[0].mass_kg*3e5)<1e-8);
    }
    #[test]
    fn adaptive_zero_step_and_short_drift() {
        let (mut particles,phys,cfg,mut state)=setup(48);
        let mut hole=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:10.0*SOLAR_MASS_KG,radius_m:3e4
        };
        let initial=particles[0].pos_m;
        advance(&mut hole,&mut particles,&mut state,phys,cfg,0.0);
        assert_eq!(particles[0].pos_m,initial);
        let p0=sph3::total_momentum(&hole,&particles);
        let mass=sph3::total_mass_kg(&hole,&particles);
        let e0=adaptive_energy_j(&hole,&particles,&mut state,phys,cfg);
        for _ in 0..20 {
            advance(&mut hole,&mut particles,&mut state,phys,cfg,0.001);
        }
        let p1=sph3::total_momentum(&hole,&particles);
        let e1=adaptive_energy_j(&hole,&particles,&mut state,phys,cfg);
        assert!((p1-p0).len()/(mass*3.0e5)<1e-9);
        assert!(e1.is_finite()&&e0.is_finite());
    }
}
