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
use crate::sph_gravity_kernel;
use crate::sph_hroot;
use crate::sph_neighbors::CellNeighbors;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum GravityModel {
    FixedPlummer,
    /// Pair potential -G m_i m_j /2 * [1/sqrt(r²+h_i²)+1/sqrt(r²+h_j²)].
    /// Requires the corresponding zeta/Omega grad-h force.
    MatchedAdaptivePlummer,
    /// Poisson potential of the same compact 3D cubic spline used for SPH
    /// density and pressure; exact Newtonian interaction for r >= 2h.
    MatchedCubicSpline,
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum DensityMethod {
    /// Existing BG-PHY-08 underrelaxed h/rho fixed-point closure.
    Underrelaxed,
    /// Bracketed Newton root of rho(h)h³=eta³m, with explicit h caps.
    SafeguardedNewton,
}
#[derive(Clone,Copy,Debug)]
pub struct AdaptiveConfig {
    pub eta:f64,
    pub min_h_m:f64,
    pub max_h_m:f64,
    pub gravity_softening_m:f64,
    pub gravity_model:GravityModel,
    pub density_iterations:usize,
    pub density_method:DensityMethod,
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
            density_method:DensityMethod::Underrelaxed,
        }
    }
    pub fn with_adaptive_gravity(mut self)->Self {
        self.gravity_model=GravityModel::MatchedAdaptivePlummer;
        self.density_iterations=24;
        self
    }
    pub fn with_cubic_gravity(mut self)->Self {
        self.gravity_model=GravityModel::MatchedCubicSpline;
        self.density_iterations=24;
        self
    }
    pub fn with_newton_density(mut self)->Self {
        self.density_method=DensityMethod::SafeguardedNewton;
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
fn density_of_particle(
    i:usize,particles:&[Fluid3],h:f64,candidates:&[usize],
)->(f64,f64) {
    let mut rho=0.0;
    let mut derivative=0.0;
    for &j in candidates {
        let x=&particles[j];
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
    if cfg.density_method==DensityMethod::SafeguardedNewton {
        // The complete admissible interval must be covered, including
        // potential cap solutions. The expensive root backend is opt-in
        // until its physical stability and runtime are separately measured.
        let neighbors=CellNeighbors::new(particles,cfg.max_h_m)
            .all_candidates();
        for i in 0..particles.len() {
            let root=sph_hroot::solve_particle(
                i,particles,&neighbors[i],cfg.eta,
                cfg.min_h_m,cfg.max_h_m,state.smoothing_m[i],
                1e-9,64);
            state.smoothing_m[i]=root.h_m;
            state.densities_kg_m3[i]=root.density_kg_m3;
            state.omega[i]=root.omega;
        }
        return;
    }
    // Build a conservative support index once per force evaluation, with
    // headroom for h changes during the closure. If an iteration exceeds
    // the cached support, rebuild before the *next* density evaluation.
    // This preserves exact compact support and avoids 24 separate cell
    // hash constructions on most force calls.
    let mut supported_h=(state.smoothing_m.iter().copied()
        .fold(0.0_f64,f64::max)*1.25)
        .clamp(cfg.min_h_m,cfg.max_h_m);
    let mut candidates=CellNeighbors::new(particles,supported_h)
        .all_candidates();
    for _ in 0..cfg.density_iterations {
        for i in 0..particles.len() {
            let h=state.smoothing_m[i];
            let (rho,_)=density_of_particle(i,particles,h,&candidates[i]);
            assert!(rho>0.0&&rho.is_finite());
            let target=cfg.eta*
                (particles[i].mass_kg/rho).cbrt();
            state.smoothing_m[i]=
                (0.45*h+0.55*target).clamp(cfg.min_h_m,cfg.max_h_m);
        }
        let evolved_h=state.smoothing_m.iter().copied()
            .fold(0.0_f64,f64::max);
        if evolved_h>supported_h {
            supported_h=(evolved_h*1.25).min(cfg.max_h_m);
            candidates=CellNeighbors::new(particles,supported_h)
                .all_candidates();
        }
    }
    // Final smoothing lengths remain inside the cached support radius.
    // A sorted candidate list yields the same summation order as all-pairs.
    for i in 0..particles.len() {
        let h=state.smoothing_m[i];
        let (rho,drho_dh)=density_of_particle(i,particles,h,&candidates[i]);
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
/// Both the Plummer and Poisson-consistent cubic potential use a negative
/// phi convention and produce a negative zeta coefficient.
pub fn adaptive_gravity_zeta(
    particles:&[Fluid3],state:&AdaptiveState,cfg:AdaptiveConfig,
)->Vec<f64> {
    assert!(state.valid_for(particles.len()));
    let mut zeta=vec![0.0;particles.len()];
    // Plummer's dphi/dh has infinite support. Only the compact M4 model
    // can safely prune the gravitational zeta summation at distance 2*h.
    let compact=if cfg.gravity_model==GravityModel::MatchedCubicSpline {
        Some(CellNeighbors::new(particles,
            state.smoothing_m.iter().copied().fold(0.0_f64,f64::max))
            .all_candidates())
    } else {None};
    for (i,p) in particles.iter().enumerate() {
        let h=state.smoothing_m[i];
        let capped=(h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8 ||
            (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8;
        if capped {continue;}
        let mut derivative_sum=0.0;
        let candidates=compact.as_ref().map(|x|x[i].clone())
            .unwrap_or_else(||(0..particles.len()).collect());
        for j in candidates {
            if i==j {continue;}
            let q=&particles[j];
            let r=(p.pos_m-q.pos_m).len();
            let derivative=match cfg.gravity_model {
                GravityModel::FixedPlummer=>0.0,
                GravityModel::MatchedAdaptivePlummer=>
                    adaptive_plummer_dphi_dh(r,h),
                GravityModel::MatchedCubicSpline=>
                    sph_gravity_kernel::smoothing_derivative(r,h),
            };
            derivative_sum+=q.mass_kg*derivative;
        }
        zeta[i]=-h/(3.0*state.densities_kg_m3[i])*derivative_sum;
    }
    zeta
}

#[derive(Debug)]
pub struct AdaptiveForces {
    pub central_acceleration_m_s2:V3,
    pub accelerations_m_s2:Vec<V3>,
    /// Pressure acceleration, including the hydrodynamic grad-h correction.
    /// Kept separate for radial hydrostatic-force diagnostics.
    pub pressure_accelerations_m_s2:Vec<V3>,
    /// Self-gravity, central field, and gravitational grad-h corrections.
    pub gravity_accelerations_m_s2:Vec<V3>,
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
    let mut pressure_acc=vec![V3::default();n];
    let mut gravity_acc=vec![V3::default();n];
    let mut central_acc=V3::default();
    let inv_eps2=cfg.gravity_softening_m.powi(2);
    let zeta=if cfg.gravity_model!=GravityModel::FixedPlummer {
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
        gravity_acc[i]-=g*central.mass_kg;
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
                GravityModel::MatchedCubicSpline=>{
                    if r>0.0 {
                        0.5*(sph_gravity_kernel::radial_derivative(
                            r,state.smoothing_m[i])+
                            sph_gravity_kernel::radial_derivative(
                            r,state.smoothing_m[j]))/r
                    } else {0.0}
                }
            };
            let grav=d*(crate::physics::G_SI*grav_scalar);
            acceleration[i]+=grav*particles[j].mass_kg;
            acceleration[j]-=grav*particles[i].mass_kg;
            gravity_acc[i]+=grav*particles[j].mass_kg;
            gravity_acc[j]-=grav*particles[i].mass_kg;
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
            pressure_acc[i]+=push*particles[j].mass_kg;
            pressure_acc[j]-=push*particles[i].mass_kg;
            if cfg.gravity_model!=GravityModel::FixedPlummer {
                // Price-Monaghan zeta correction derived by differentiating
                // the pair gravitational potential through h_i(rho_i).
                let gcoeff=0.5*crate::physics::G_SI*
                    (zeta[i]/state.omega[i]*grad_i+
                     zeta[j]/state.omega[j]*grad_j)/r;
                let correction=d*gcoeff;
                acceleration[i]+=correction*particles[j].mass_kg;
                acceleration[j]-=correction*particles[i].mass_kg;
                gravity_acc[i]+=correction*particles[j].mass_kg;
                gravity_acc[j]-=correction*particles[i].mass_kg;
            }
        }
    }
    AdaptiveForces {
        central_acceleration_m_s2:central_acc,
        accelerations_m_s2:acceleration,
        pressure_accelerations_m_s2:pressure_acc,
        gravity_accelerations_m_s2:gravity_acc,
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
                GravityModel::MatchedCubicSpline=>{
                    let radius=r2.sqrt();
                    -0.5*(sph_gravity_kernel::potential(
                        radius,state.smoothing_m[i])+
                        sph_gravity_kernel::potential(
                        radius,state.smoothing_m[j]))
                },
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
    fn cell_list_density_matches_full_pairwise_sum() {
        let (particles,_,cfg,mut state)=setup(128);
        update_adaptive_densities(&particles,&mut state,cfg);
        for i in 0..particles.len() {
            let h=state.smoothing_m[i];
            let mut full=0.0;
            let mut deriv=0.0;
            for j in 0..particles.len() {
                let r=(particles[i].pos_m-particles[j].pos_m).len();
                let (w,dwdh)=kernel_and_derivative_h(r,h);
                full+=particles[j].mass_kg*w;
                deriv+=particles[j].mass_kg*dwdh;
            }
            let diff=(full-state.densities_kg_m3[i]).abs()/full;
            assert!(diff<1e-13,"cell density mismatch particle {i}: {diff}");
            let expected_omega=if (h-cfg.min_h_m).abs()<cfg.min_h_m*1e-8||
                (h-cfg.max_h_m).abs()<cfg.max_h_m*1e-8 {
                1.0
            } else {(1.0+h*deriv/(3.0*full)).max(0.02)};
            assert!((state.omega[i]-expected_omega).abs()<1e-11);
        }
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
    fn matched_energy_gradient_agrees_with_variational_force() {
        let (mut particles,physics,old,_)=setup(32);
        let cfg=old.with_adaptive_gravity();
        let center=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let mut state=AdaptiveState::new(&particles,8.4e6,cfg);
        let force=forces(&center,&particles,&mut state,physics,cfg);
        let predicted=force.accelerations_m_s2[0].x;
        let step_m=100.0;
        particles[0].pos_m.x+=step_m;
        let mut plus=AdaptiveState::new(&particles,8.4e6,cfg);
        let e_plus=adaptive_energy_j(&center,&particles,&mut plus,physics,cfg);
        particles[0].pos_m.x-=2.0*step_m;
        let mut minus=AdaptiveState::new(&particles,8.4e6,cfg);
        let e_minus=adaptive_energy_j(&center,&particles,&mut minus,physics,cfg);
        let from_energy=-(e_plus-e_minus)/(2.0*step_m*particles[0].mass_kg);
        let scale=force.accelerations_m_s2[0].len().max(1.0e3);
        let disagreement=(predicted-from_energy).abs()/scale;
        assert!(disagreement<0.04,
            "energy gradient force mismatch relative {disagreement}, analytic {predicted}, finite {from_energy}");
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
    fn safeguarded_newton_closure_matches_density_and_caps() {
        let (p,phys,cfg_old,mut state)=setup(64);
        let cfg=cfg_old.with_cubic_gravity().with_newton_density();
        update_adaptive_densities(&p,&mut state,cfg);
        assert!(state.densities_kg_m3.iter().all(|x|x.is_finite()&&*x>0.0));
        assert!(state.omega.iter().all(|x|x.is_finite()&&*x>=0.02));
        let full=(0..p.len()).collect::<Vec<_>>();
        for i in 0..p.len() {
            let h=state.smoothing_m[i];
            let is_cap=(h-cfg.min_h_m).abs()<cfg.min_h_m*1e-8 ||
                (h-cfg.max_h_m).abs()<cfg.max_h_m*1e-8;
            let residual=sph_hroot::closure_residual(
                i,&p,&full,cfg.eta,h);
            if !is_cap {
                assert!(residual<1e-8,
                    "particle {i} h-rho residual {residual}");
            } else {
                assert_eq!(state.omega[i],1.0);
            }
            let direct=p.iter().map(|j|
                j.mass_kg*sph3::cubic_spline_3d(
                    (p[i].pos_m-j.pos_m).len(),h)).sum::<f64>();
            assert!((direct-state.densities_kg_m3[i]).abs()/direct<1e-12);
        }
        // A physical force evaluation must also use the new method.
        let mut new_state=state.clone();
        let central=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let f=forces(&central,&p,&mut new_state,phys,cfg);
        assert!(f.accelerations_m_s2.iter().all(|a|a.finite()));
    }
    #[test]
    fn cubic_kernel_gravity_is_newtonian_beyond_compact_support() {
        let (mut particles,phys,old,_)=setup(32);
        let cfg=old.with_cubic_gravity();
        let mut state=AdaptiveState::new(&particles,8.4e6,cfg);
        let mut central=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        // Separate just two equal masses by 12h in otherwise empty space.
        particles.truncate(2);
        particles[0].pos_m=V3::new(0.0,0.0,0.0);
        particles[1].pos_m=V3::new(40.0e6,0.0,0.0);
        state=AdaptiveState::new(&particles,8.4e6,cfg);
        let force=forces(&central,&particles,&mut state,phys,cfg);
        let r=40.0e6;
        let expected=crate::physics::G_SI*particles[1].mass_kg/(r*r);
        // Pressure has compact support too, so the pair responds to
        // exactly Newtonian self-gravity at this separation.
        let got=force.accelerations_m_s2[0].x;
        assert!((got-expected).abs()/expected<1e-12);
        let zeta=adaptive_gravity_zeta(&particles,&state,cfg);
        assert!(zeta.iter().all(|z|z.abs()<1e-14));
        central.pos_m=V3::new(2e8,0.0,0.0);
    }
    #[test]
    fn cubic_kernel_gravity_force_matches_energy_gradient() {
        let (mut particles,phys,old,_)=setup(48);
        let cfg=old.with_cubic_gravity();
        let central=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let mut state=AdaptiveState::new(&particles,8.4e6,cfg);
        let a=forces(&central,&particles,&mut state,phys,cfg)
            .accelerations_m_s2[2].x;
        let delta=30.0;
        particles[2].pos_m.x+=delta;
        let mut positive=state.clone();
        let e_plus=adaptive_energy_j(&central,&particles,
            &mut positive,phys,cfg);
        particles[2].pos_m.x-=2.0*delta;
        let mut negative=state.clone();
        let e_minus=adaptive_energy_j(&central,&particles,
            &mut negative,phys,cfg);
        let numeric=-(e_plus-e_minus)/(2.0*delta*particles[2].mass_kg);
        let relative=(a-numeric).abs()/a.abs().max(1.0e3);
        assert!(relative<0.03,
            "cubic gravity variational energy gradient error {relative}, analytic={a}, numeric={numeric}");
    }
    #[test]
    fn cubic_gravity_preserves_momentum_and_short_step_energy() {
        let (mut particles,phys,old,_)=setup(64);
        let cfg=old.with_cubic_gravity();
        let mut state=AdaptiveState::new(&particles,8.4e6,cfg);
        let mut central=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let p0=sph3::total_momentum(&central,&particles);
        let e0=adaptive_energy_j(&central,&particles,&mut state,phys,cfg);
        for _ in 0..60 {
            advance(&mut central,&mut particles,&mut state,phys,cfg,0.002);
        }
        let p1=sph3::total_momentum(&central,&particles);
        let e1=adaptive_energy_j(&central,&particles,&mut state,phys,cfg);
        assert!((p1-p0).len()/(particles[0].mass_kg*3e5)<1e-8);
        assert!((e1-e0).abs()/e0.abs()<1e-3);
    }
    #[test]
    fn force_decomposition_reconstructs_total_and_pair_momentum() {
        for model in [
            GravityModel::FixedPlummer,
            GravityModel::MatchedAdaptivePlummer,
            GravityModel::MatchedCubicSpline,
        ] {
            let (p,phys,mut cfg,_)=setup(48);
            cfg.gravity_model=model;
            cfg.density_iterations=24;
            let mut state=AdaptiveState::new(&p,8.4e6,cfg);
            let central=Massive3{
                pos_m:V3::default(),vel_m_s:V3::default(),
                mass_kg:0.0,radius_m:0.0,
            };
            let f=forces(&central,&p,&mut state,phys,cfg);
            let mut net_pressure=V3::default();
            let mut net_gravity=V3::default();
            for (i,part) in p.iter().enumerate() {
                let reconstructed=f.pressure_accelerations_m_s2[i]+
                    f.gravity_accelerations_m_s2[i];
                let err=(reconstructed-f.accelerations_m_s2[i]).len();
                assert!(err/f.accelerations_m_s2[i].len().max(1.0)<1e-11);
                net_pressure+=f.pressure_accelerations_m_s2[i]*part.mass_kg;
                net_gravity+=f.gravity_accelerations_m_s2[i]*part.mass_kg;
            }
            let scale=p[0].mass_kg*f.accelerations_m_s2[0].len().max(1.0);
            assert!(net_pressure.len()/scale<1e-9);
            assert!(net_gravity.len()/scale<1e-9);
        }
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
