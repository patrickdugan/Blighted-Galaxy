//! BG-PHY-22: analytic total derivative of the cold EOS energy with
//! BG-PHY-20's partition-gated affine MLS density, including implicit h.
//!
//! This is a small-N O(N^2) reverse/adjoint calculation of
//!     U(x) = sum_i m_i u_EOS(rho_hat_i(x)).
//! Every SPH density, smoothing length h, source volume, partition
//! S, MLS coefficient and smoothstep exposure participates in dU/dx.
//! The derivative is checked against BG-PHY-21B's *independent* complete
//! finite-difference re-solve, not against a differentiable approximation
//! of the old pair pressure force.
//!
//! The kernel is compact and C1 at its support edge. Away from changes
//! in the h-cap/MLS-fallback active set the derivative is smooth.
//! A hard cap has dh/dx=0; a floored Omega or MLS rank change is flagged.
//!
//! This module remains an opt-in, gravity-free barotropic research pilot.
//! It does not substitute forces into the shipped gameplay integrator.
use crate::sph3::{self,Fluid3,Sph3Config,V3};
use crate::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use crate::sph_density_audit::{self,DensityEstimator};
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_rk_compatible::DensityChoice;

type M4=[[f64;4];4];
fn dot4(a:[f64;4],b:[f64;4])->f64 {
    (0..4).map(|i|a[i]*b[i]).sum()
}
fn phi(q:V3)->[f64;4]{[1.0,q.x,q.y,q.z]}
fn cross(a:V3,b:V3)->V3 {
    V3::new(a.y*b.z-a.z*b.y,
        a.z*b.x-a.x*b.z,
        a.x*b.y-a.y*b.x)
}
fn kernel_terms(d:V3,h:f64,k:DensityKernel)->(f64,f64,V3){
    let r=d.len();
    let (w,dh)=sph_kernel::kernel_and_dh(r,h,k);
    let g=if r>0.0 {
        d*(sph_kernel::radial_derivative(r,h,k)/r)
    } else {V3::default()};
    (w,dh,g)
}
#[derive(Clone,Copy)]
struct Term {
    j:usize,
    q:V3,
    w:f64,
    kernel_w:f64,
    kernel_dh:f64,
    grad_w:V3,
    volume:f64,
}
/// Diagonally scaled Cholesky of the same affine MLS normal matrix as
/// sph_density_audit::weighted_intercept. Return the intercept fit beta
/// and z=M^{-1}e0, required for d(beta[0])=z·(db-dM beta).
fn affine_solve(terms:&[Term],rho:&[f64])
    ->Option<([f64;4],[f64;4],f64)> {
    if terms.len()<6{return None;}
    let mut matrix:M4=[[0.0;4];4];
    let mut rhs=[0.0;4];
    for t in terms {
        let x=phi(t.q);
        for a in 0..4 {
            rhs[a]+=t.w*x[a]*rho[t.j];
            for b in 0..4 {matrix[a][b]+=t.w*x[a]*x[b];}
        }
    }
    let mut scale=[0.0;4];
    for a in 0..4 {
        if !matrix[a][a].is_finite()||matrix[a][a]<=0.0 {
            return None;
        }
        scale[a]=matrix[a][a].sqrt();
    }
    let mut l:M4=[[0.0;4];4];
    let (mut pmin,mut pmax)=(f64::INFINITY,0.0_f64);
    for a in 0..4 {
        for b in 0..=a{
            let mut t=matrix[a][b]/(scale[a]*scale[b]);
            for k in 0..b {t-=l[a][k]*l[b][k];}
            if a==b {
                if !t.is_finite()||t<=1e-10{return None;}
                pmin=pmin.min(t);
                pmax=pmax.max(t);
                l[a][a]=t.sqrt();
            } else {l[a][b]=t/l[b][b];}
        }
    }
    let spread=pmax/pmin;
    if !spread.is_finite()||spread>1e8{return None;}
    let solve=|b:[f64;4]| ->[f64;4] {
        let mut y=[0.0;4];
        for a in 0..4{
            let mut t=b[a]/scale[a];
            for k in 0..a{t-=l[a][k]*y[k];}
            y[a]=t/l[a][a];
        }
        let mut x=[0.0;4];
        for a in (0..4).rev(){
            let mut t=y[a];
            for k in a+1..4{t-=l[k][a]*x[k];}
            x[a]=t/l[a][a];
        }
        for a in 0..4{x[a]/=scale[a];}
        x
    };
    let beta=solve(rhs);
    if !beta[0].is_finite()||beta[0]<=0.0{return None;}
    let z=solve([1.0,0.0,0.0,0.0]);
    Some((beta,z,spread))
}
#[derive(Clone,Debug)]
pub struct AnalyticEnergyGradient {
    pub forces_n:Vec<V3>,
    pub potential_j:f64,
    pub density_kg_m3:Vec<f64>,
    pub net_force_ratio:f64,
    pub net_torque_ratio:f64,
    /// RMS magnitude of direct kernel and basis displacement terms.
    pub direct_force_rms_n:f64,
    /// RMS magnitude of dependence on the baseline SPH rho closure.
    pub density_force_rms_n:f64,
    /// RMS magnitude of implicit smoothing-length terms.
    pub smoothing_force_rms_n:f64,
    pub mls_fallback_fraction:f64,
    pub max_mls_relative_disagreement:f64,
    pub minimum_omega:f64,
    pub capped_fraction:f64,
    pub omega_floor_fraction:f64,
}
fn force_rms(a:&[V3])->f64 {
    (a.iter().map(|x|x.len_sq()).sum::<f64>()/a.len() as f64).sqrt()
}
/// Compute the true analytic -dU/dx for the *specific* BG-PHY-20
/// surface-gated affine MLS construction, with all global particle
/// dependencies and h_i implicit derivatives included.
pub fn total_energy_gradient(
    p:&[Fluid3],physics:Sph3Config,cfg:AdaptiveConfig,
    reference_radius_m:f64,mode:DensityChoice,
)->AnalyticEnergyGradient {
    let n=p.len();
    // BG-PHY-23 extends the validated analytic adjoint to the N256
    // research-only stellar pilot. Gameplay paths are unchanged.
    assert!((16..=256).contains(&n)&&cfg.valid()&&
        physics.valid()&&reference_radius_m>0.0);
    let exact_cfg=cfg.with_newton_density();
    let mut state=AdaptiveState::new(p,reference_radius_m,exact_cfg);
    sph_adaptive::update_adaptive_densities(p,&mut state,exact_cfg);
    let rho=&state.densities_kg_m3;
    let mut field=rho.to_vec();
    let mut reference_blend=None;
    if mode==DensityChoice::SurfaceBlended {
        let report=sph_density_audit::estimate(p,&state,exact_cfg.kernel);
        let v=report.fields.iter().find(|x|
            x.method==DensityEstimator::SurfaceBlendedMls)
            .expect("surface blend must be present").rho.clone();
        field=v;
        reference_blend=Some(report);
    }
    let u=field.iter().zip(p).map(|(rho,x)|x.mass_kg*
        sph3::internal_energy_j_per_kg(*rho,physics)).sum::<f64>();
    let l=field.iter().zip(p).map(|(rho,x)|
        x.mass_kg*sph3::pressure_pa(*rho,physics)/(rho*rho))
        .collect::<Vec<_>>();
    let mut adj_rho=vec![0.0_f64;n];
    let mut adj_h=vec![0.0_f64;n];
    let mut grad_x=vec![V3::default();n];
    let mut invalid=0_usize;
    let mut max_disagreement=0.0_f64;
    if mode==DensityChoice::Original {
        adj_rho.copy_from_slice(&l);
    } else {
        let batch=reference_blend.as_ref().unwrap();
        for i in 0..n{
            let h=state.smoothing_m[i];
            let mut terms=Vec::new();
            let mut partition=0.0;
            for j in 0..n {
                let d=p[i].pos_m-p[j].pos_m;
                let (w,dh,grad_w)=kernel_terms(
                    d,h,exact_cfg.kernel);
                if w<=0.0 {continue;}
                let volume=p[j].mass_kg/rho[j];
                let weighted=volume*w;
                partition+=weighted;
                terms.push(Term{
                    j,q:(p[j].pos_m-p[i].pos_m)*(1.0/h),
                    w:weighted,kernel_w:w,
                    kernel_dh:dh,grad_w,volume,
                });
            }
            let fitted=affine_solve(&terms,rho);
            if let Some((beta,z,_))=fitted {
                let t=((1.0-partition)/0.5).clamp(0.0,1.0);
                let alpha=t*t*(3.0-2.0*t);
                let alpha_prime=if t>0.0&&t<1.0 {
                    -12.0*t*(1.0-t)
                } else {0.0};
                let expected=(1.0-alpha)*rho[i]+alpha*beta[0];
                let discrepancy=(expected-field[i]).abs()/
                    field[i].max(1e-100);
                max_disagreement=max_disagreement.max(discrepancy);
                // The independent evaluator must agree with the identical
                // BG-PHY-20 arithmetic and fallback branch.
                assert!(discrepancy<1e-9,
                    "MLS scalar mismatch i={i}, rel={discrepancy}");
                assert!((batch.partition[i]-partition).abs()<1e-9);
                adj_rho[i]+=l[i]*(1.0-alpha);
                let dbeta=l[i]*alpha;
                let ds=l[i]*(beta[0]-rho[i])*alpha_prime;
                for t in &terms {
                    let j=t.j;
                    let x=phi(t.q);
                    let zx=dot4(z,x);
                    let bx=dot4(beta,x);
                    let residual=rho[j]-bx;
                    // d beta0 = z^T(db-dM beta) includes dM^{-1}.
                    let dw_coeff=ds+dbeta*zx*residual;
                    let drho_coeff=dbeta*t.w*zx-
                        dw_coeff*t.w/rho[j];
                    adj_rho[j]+=drho_coeff;
                    let dq=V3::new(z[1],z[2],z[3])*
                        (dbeta*t.w*residual)-
                        V3::new(beta[1],beta[2],beta[3])*
                        (dbeta*t.w*zx);
                    // dW = grad_i W ·(dx_i-dx_j)+dW/dh dh_i
                    // dq = (dx_j-dx_i)/h_i-q dh_i/h_i.
                    let direct=t.grad_w*(dw_coeff*t.volume)-
                        dq*(1.0/h);
                    grad_x[i]+=direct;
                    grad_x[j]-=direct;
                    adj_h[i]+=dw_coeff*t.volume*t.kernel_dh-
                        dq.dot(t.q)/h;
                }
            } else {
                invalid+=1;
                // An invalid/singular/nonpositive affine fit is a
                // piecewise baseline branch in sph_density_audit.
                assert!((field[i]-rho[i]).abs()/rho[i]<1e-12);
                adj_rho[i]+=l[i];
            }
        }
    }
    let mut grad_density=vec![V3::default();n];
    let mut grad_smoothing=vec![V3::default();n];
    let mut caps=0_usize;
    let mut floors=0_usize;
    for i in 0..n {
        let h=state.smoothing_m[i];
        let r=rho[i];
        let omega=state.omega[i];
        let capped=(h-exact_cfg.min_h_m).abs()<=exact_cfg.min_h_m*1e-8||
            (h-exact_cfg.max_h_m).abs()<=exact_cfg.max_h_m*1e-8;
        if capped{caps+=1;}
        if omega<=0.0200001 {floors+=1;}
        // Drho_total = drho_fixed_h/Omega away from h caps.
        // dh = -h/(3 rho Omega)*drho_fixed_h.
        // The hard cap uses dh=0 and Omega=1 by construction.
        let rho_coeff=adj_rho[i]/omega;
        let h_coeff=if capped{0.0} else {
            -adj_h[i]*h/(3.0*r*omega)
        };
        for j in 0..n{
            if i==j{continue;}
            let d=p[i].pos_m-p[j].pos_m;
            let r=d.len();
            if r>=2.0*h||r==0.0{continue;}
            let g=d*(sph_kernel::radial_derivative(
                r,h,exact_cfg.kernel)/r*p[j].mass_kg);
            let dg=g*rho_coeff;
            let hg=g*h_coeff;
            grad_density[i]+=dg;
            grad_density[j]-=dg;
            grad_smoothing[i]+=hg;
            grad_smoothing[j]-=hg;
        }
    }
    let mut forces=Vec::with_capacity(n);
    let mut net=V3::default();
    let mut torque=V3::default();
    let mut sum_f=0.0;
    let mut sum_t=0.0;
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),
        |a,x|a+x.pos_m*x.mass_kg)/mass;
    for i in 0..n{
        let force=(grad_x[i]+grad_density[i]+grad_smoothing[i])*(-1.0);
        forces.push(force);
        net+=force;
        let arm=p[i].pos_m-center;
        torque+=cross(arm,force);
        sum_f+=force.len();
        sum_t+=arm.len()*force.len();
    }
    AnalyticEnergyGradient{
        forces_n:forces,potential_j:u,density_kg_m3:field,
        net_force_ratio:net.len()/sum_f.max(1.0),
        net_torque_ratio:torque.len()/sum_t.max(1.0),
        direct_force_rms_n:force_rms(&grad_x),
        density_force_rms_n:force_rms(&grad_density),
        smoothing_force_rms_n:force_rms(&grad_smoothing),
        mls_fallback_fraction:invalid as f64/n as f64,
        max_mls_relative_disagreement:max_disagreement,
        minimum_omega:state.omega.iter().copied().fold(f64::INFINITY,f64::min),
        capped_fraction:caps as f64/n as f64,
        omega_floor_fraction:floors as f64/n as f64,
    }
}
/// Pressure-only Hamiltonian velocity Verlet using the analytic
/// -grad_x sum m u_EOS(rho_hat(x)). No independent compatible-U reservoir.
#[derive(Clone)]
pub struct AnalyticPilot {
    pub particles:Vec<Fluid3>,
    pub physics:Sph3Config,
    pub cfg:AdaptiveConfig,
    pub reference_radius_m:f64,
    pub density:DensityChoice,
}
impl AnalyticPilot {
    pub fn gradient(&self)->AnalyticEnergyGradient{
        total_energy_gradient(&self.particles,self.physics,self.cfg,
            self.reference_radius_m,self.density)
    }
    pub fn energy_j(&self)->f64{
        let thermal=self.gradient().potential_j;
        let kinetic=self.particles.iter().map(|x|
            0.5*x.mass_kg*x.vel_m_s.len_sq()).sum::<f64>();
        thermal+kinetic
    }
    pub fn advance(&mut self,dt:f64){
        assert!(dt.is_finite()&&dt.abs()<=0.02);
        if dt==0.0{return;}
        let f0=self.gradient();
        for (p,f) in self.particles.iter_mut().zip(f0.forces_n){
            p.vel_m_s+=f*(0.5*dt/p.mass_kg);
            p.pos_m+=p.vel_m_s*dt;
        }
        let f1=self.gradient();
        for (p,f) in self.particles.iter_mut().zip(f1.forces_n){
            p.vel_m_s+=f*(0.5*dt/p.mass_kg);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_rk_variational;
    use crate::white_dwarf::{profile_for_mass,DegenerateElectrons};
    fn fixture(n:usize,kernel:DensityKernel)
        ->(Vec<Fluid3>,Sph3Config,AdaptiveConfig,f64){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,
            V3::new(1.0e8,0.0,0.0),V3::default(),n,17).0;
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==DensityKernel::WendlandC4{cfg=cfg.with_wendland_c4();}
        cfg.eta=1.7;
        (p,Sph3Config::cold_white_dwarf(2.0).unwrap(),cfg,star.radius_m)
    }
    #[test]
    fn analytic_original_recovers_independent_total_energy_force(){
        for k in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let (p,phys,cfg,radius)=fixture(16,k);
            let a=total_energy_gradient(&p,phys,cfg,radius,DensityChoice::Original);
            let b=sph_rk_variational::total_gradient(&p,phys,cfg,radius,
                DensityChoice::Original,1e-4);
            let mut num=0.0;
            let mut denom=0.0;
            for (fa,fb) in a.forces_n.iter().zip(&b.forces_n){
                num+=(*fa-*fb).len_sq();denom+=fb.len_sq();
            }
            let err=(num/denom).sqrt();
            assert!(err<2e-5,"original {:?}: {err:e}",k);
            assert!(a.net_force_ratio<1e-10);
            assert!(a.net_torque_ratio<1e-10);
        }
    }
    #[test]
    fn analytic_surface_blend_matches_complete_perturbation_oracle(){
        for k in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let (p,phys,cfg,radius)=fixture(16,k);
            let a=total_energy_gradient(&p,phys,cfg,radius,DensityChoice::SurfaceBlended);
            let b=sph_rk_variational::total_gradient(&p,phys,cfg,radius,
                DensityChoice::SurfaceBlended,1e-4);
            let mut num=0.0;
            let mut denom=0.0;
            for (fa,fb) in a.forces_n.iter().zip(&b.forces_n){
                num+=(*fa-*fb).len_sq();denom+=fb.len_sq();
            }
            let err=(num/denom).sqrt();
            assert!(err<5e-5,"surface-blend {:?}: {err:e}",k);
            assert!(a.max_mls_relative_disagreement<1e-9);
            assert!(a.net_force_ratio<1e-9);
            assert!(a.net_torque_ratio<1e-9);
        }
    }
}
