//! Safeguarded Newton-bisection for the implicit SPH density/h relation.
 //!
 //! For each particle i, solve rho_i(h) h^3 = eta^3 m_i.
 //! The derivative is 3 h^2 rho_i Omega_i, with the same grad-h Omega
 //! used by the pressure and gravitational zeta corrections.
 //!
 //! A candidate list must cover the whole allowed interval [h_min,h_max].
 //! If the root lies outside it, return the appropriate cap explicitly.
 //! All root and residual diagnostics are independent of the fixed-iteration
 //! underrelaxed method preserved in sph_adaptive.rs.
use crate::sph3::{Fluid3,V3};
use crate::sph_kernel::{self,DensityKernel};

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum RootLocation { Interior, LowerBound, UpperBound }

#[derive(Clone,Copy,Debug)]
pub struct RootSolve {
    pub h_m:f64,
    pub density_kg_m3:f64,
    pub density_derivative_h:f64,
    pub omega:f64,
    pub iterations:usize,
    /// |rho h^3 / (eta^3 m) -1|, dimensionless. A constrained solution
    /// can retain a large residual when its unconstrained root is outside
    /// the configured smoothing bounds.
    pub residual:f64,
    pub location:RootLocation,
}
fn density_and_dh(
    i:usize,p:&[Fluid3],indices:&[usize],h:f64,kernel:DensityKernel
)->(f64,f64) {
    let mut density=0.0;
    let mut derivative=0.0;
    for &j in indices {
        let distance=(p[i].pos_m-p[j].pos_m).len();
        let (w,dwdh)=sph_kernel::kernel_and_dh(distance,h,kernel);
        density+=p[j].mass_kg*w;
        derivative+=p[j].mass_kg*dwdh;
    }
    (density,derivative)
}
fn evaluate(i:usize,p:&[Fluid3],candidates:&[usize],eta:f64,h:f64,
    kernel:DensityKernel)
    ->(f64,f64,f64,f64) {
    let (rho,drho)=density_and_dh(i,p,candidates,h,kernel);
    assert!(rho.is_finite()&&rho>0.0);
    let normalization=(eta.powi(3)*p[i].mass_kg).max(1e-50);
    let f=h.powi(3)*rho/normalization-1.0;
    let fp=h*h*(3.0*rho+h*drho)/normalization;
    (f,fp,rho,drho)
}
fn finish(i:usize,p:&[Fluid3],candidates:&[usize],
    eta:f64,h:f64,n:usize,location:RootLocation,
    kernel:DensityKernel)->RootSolve {
    let (f,_,rho,drho)=evaluate(i,p,candidates,eta,h,kernel);
    let omega=if location==RootLocation::Interior {
        (1.0+h*drho/(3.0*rho)).max(0.02)
    } else {1.0};
    RootSolve {
        h_m:h,density_kg_m3:rho,density_derivative_h:drho,
        omega,iterations:n,residual:f.abs(),location,
    }
}

/// All supplied candidates must contain every particle within 2*h_max
/// of the queried particle. Sorted candidates preserve reproducible
/// density summation. This function does not update particle velocities.
pub fn solve_particle(
    i:usize,p:&[Fluid3],candidates:&[usize],
    eta:f64,h_min:f64,h_max:f64,h_initial:f64,
    relative_tolerance:f64,max_iter:usize,
)->RootSolve {
    solve_particle_with_kernel(i,p,candidates,eta,h_min,h_max,h_initial,
        relative_tolerance,max_iter,DensityKernel::CubicM4)
}
/// Independent kernel selection while retaining the M4 entry point.
pub fn solve_particle_with_kernel(
    i:usize,p:&[Fluid3],candidates:&[usize],
    eta:f64,h_min:f64,h_max:f64,h_initial:f64,
    relative_tolerance:f64,max_iter:usize,kernel:DensityKernel,
)->RootSolve {
    assert!(i<p.len()&&p[i].mass_kg>0.0&&
        p[i].pos_m.finite()&&h_min>0.0&&h_max>h_min&&
        eta.is_finite()&&eta>0.5&&
        relative_tolerance.is_finite()&&relative_tolerance>0.0&&
        relative_tolerance<=1e-3&&max_iter>0);
    let (flo,_,_,_)=evaluate(i,p,candidates,eta,h_min,kernel);
    if flo>=0.0 {
        return finish(i,p,candidates,eta,h_min,0,RootLocation::LowerBound,kernel);
    }
    let (fhi,_,_,_)=evaluate(i,p,candidates,eta,h_max,kernel);
    if fhi<=0.0 {
        return finish(i,p,candidates,eta,h_max,0,RootLocation::UpperBound,kernel);
    }
    let mut lo=h_min;let mut hi=h_max;
    let mut h=h_initial.clamp(lo,hi);
    let mut iterations=0;
    for step in 0..max_iter {
        let (f,fp,_,_)=evaluate(i,p,candidates,eta,h,kernel);
        iterations=step+1;
        if f.abs()<=relative_tolerance {
            return finish(i,p,candidates,eta,h,iterations,RootLocation::Interior,kernel);
        }
        if f>0.0 {hi=h;} else {lo=h;}
        let newton=if fp.is_finite()&&fp>0.0 {h-f/fp}
            else {f64::NAN};
        // A Newton step very close to a bracket endpoint is often
        // inefficient; bisect when the proposed step leaves the safe
        // central 98% of the bracket.
        h=if newton.is_finite()&&
            newton>lo+0.01*(hi-lo)&&
            newton<hi-0.01*(hi-lo) {newton}
            else {0.5*(lo+hi)};
    }
    finish(i,p,candidates,eta,h,iterations,RootLocation::Interior,kernel)
}

/// Residual of the h-rho constraint at a preexisting smoothing length.
pub fn closure_residual(
    i:usize,p:&[Fluid3],candidates:&[usize],eta:f64,h:f64
)->f64 {
    closure_residual_with_kernel(i,p,candidates,eta,h,DensityKernel::CubicM4)
}
pub fn closure_residual_with_kernel(
    i:usize,p:&[Fluid3],candidates:&[usize],eta:f64,h:f64,
    kernel:DensityKernel,
)->f64 {
    evaluate(i,p,candidates,eta,h,kernel).0.abs()
}

#[derive(Clone,Copy,Debug)]
pub struct RootAudit {
    pub particle_count:usize,
    pub rms_equation_residual:f64,
    pub max_equation_residual:f64,
    pub max_unconstrained_residual:f64,
    pub rms_iterations:f64,
    pub capped_fraction:f64,
    pub max_omega_delta:f64,
    pub min_omega:f64,
}
pub fn summarize(solves:&[RootSolve])->RootAudit {
    assert!(!solves.is_empty());
    let n=solves.len();
    let mut residual_sq=0.0;
    let mut iterations_sq=0.0;
    let mut max_residual=0.0_f64;
    let mut max_free=0.0_f64;
    let mut min_omega=f64::INFINITY;
    let mut caps=0;
    for s in solves {
        residual_sq+=s.residual*s.residual;
        iterations_sq+=(s.iterations as f64).powi(2);
        max_residual=max_residual.max(s.residual);
        if s.location==RootLocation::Interior {
            max_free=max_free.max(s.residual);
        } else {caps+=1;}
        min_omega=min_omega.min(s.omega);
    }
    RootAudit{
        particle_count:n,
        rms_equation_residual:(residual_sq/n as f64).sqrt(),
        max_equation_residual:max_residual,
        max_unconstrained_residual:max_free,
        rms_iterations:(iterations_sq/n as f64).sqrt(),
        capped_fraction:caps as f64/n as f64,
        max_omega_delta:0.0,
        min_omega,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn positions()->Vec<Fluid3> {
        // Dense, roughly uniform 5^3 lattice with true positive self-density.
        let mut p=Vec::new();
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    p.push(Fluid3{
                        pos_m:V3::new(x as f64*1e6,y as f64*1e6,z as f64*1e6),
                        vel_m_s:V3::default(),mass_kg:1e27,
                    });
                }
            }
        }
        p
    }
    #[test]
    fn exact_lattice_root_is_small_and_density_matches_source() {
        let p=positions();let candidates=(0..p.len()).collect::<Vec<_>>();
        let central=p.len()/2;
        let root=solve_particle(
            central,&p,&candidates,1.2,0.3e6,5.0e6,1.6e6,1e-10,80);
        assert_eq!(root.location,RootLocation::Interior);
        assert!(root.residual<1e-10);
        assert!(root.omega>0.0&&root.omega<=2.0);
        let full=p.iter().map(|q|q.mass_kg*
            crate::sph3::cubic_spline_3d(
                (p[central].pos_m-q.pos_m).len(),root.h_m))
            .sum::<f64>();
        assert!((full-root.density_kg_m3).abs()/full<1e-13);
        assert!(closure_residual(central,&p,&candidates,1.2,
            root.h_m)<1e-10);
    }
    #[test]
    fn constrained_root_reports_boundary_and_equation_error() {
        let p=positions();let indices=(0..p.len()).collect::<Vec<_>>();
        let i=62;
        let lo=solve_particle(i,&p,&indices,1.2,
            2.8e6,4e6,3e6,1e-10,80);
        assert_eq!(lo.location,RootLocation::LowerBound);
        let hi=solve_particle(i,&p,&indices,1.2,
            0.2e6,0.65e6,0.4e6,1e-10,80);
        assert_eq!(hi.location,RootLocation::UpperBound);
        assert!(lo.residual>0.0&&hi.residual>0.0);
        assert_eq!(lo.omega,1.0);
    }
    #[test]
    fn root_solver_is_deterministic_and_start_invariant() {
        let p=positions();let candidates=(0..p.len()).collect::<Vec<_>>();
        let a=solve_particle(62,&p,&candidates,1.2,
            0.4e6,4e6,0.45e6,1e-11,80);
        let b=solve_particle(62,&p,&candidates,1.2,
            0.4e6,4e6,3.9e6,1e-11,80);
        assert_eq!(a.location,RootLocation::Interior);
        assert_eq!(b.location,RootLocation::Interior);
        assert!((a.h_m-b.h_m).abs()/a.h_m<1e-10);
    }
    #[test]
    fn implicit_root_can_cover_seed_with_zero_velocity() {
        let p=positions();
        let candidates=(0..p.len()).collect::<Vec<_>>();
        let solutions=(0..p.len()).map(|i|
            solve_particle(i,&p,&candidates,1.2,0.3e6,5e6,
                1.6e6,1e-8,64)).collect::<Vec<_>>();
        let audit=summarize(&solutions);
        assert!(audit.max_unconstrained_residual<1e-8);
        assert!(audit.min_omega>0.0);
    }
}
