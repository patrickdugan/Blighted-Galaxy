//! BG-PHY-17: kernel zeroth/first-order SPH consistency at frozen positions.
//!
//! E0_i = h_i * sum_j (m_j/rho_j) grad_i W_ij(h_i)
//! measures the spurious constant-field gradient, independent of EOS.
//! The first-moment trace samples the linear-gradient consistency:
//! T_i = sum_j (m_j/rho_j) ((r_j-r_i) dot grad_i W_ij)/3,
//! expected to approach unity for an isotropic well-sampled interior.
//! Boundary bias and angular noise are reported separately by radial quartile.
use crate::sph3::{Fluid3,V3};
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_neighbors::CellNeighbors;
#[derive(Clone,Copy,Debug)]
pub struct ConsistencyMetrics {
    pub e0_rms:f64,
    pub e0_inner_rms:f64,
    pub e0_outer_rms:f64,
    pub e0_outer_mean_radial:f64,
    pub e0_inner_tangent_rms:f64,
    pub e1_trace_mean:f64,
    pub e1_trace_outer_mean:f64,
}
pub fn audit(
    p:&[Fluid3],state:&AdaptiveState,kernel:DensityKernel
)->ConsistencyMetrics {
    assert!(!p.is_empty()&&state.valid_for(p.len()));
    let n=p.len();
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/mass;
    let mut sorted=(0..n).collect::<Vec<_>>();
    sorted.sort_by(|&i,&j|
        (p[i].pos_m-center).len_sq()
            .total_cmp(&(p[j].pos_m-center).len_sq()).then(i.cmp(&j)));
    let mut e0=vec![V3::default();n];
    let mut traces=vec![0.0;n];
    let hmax=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let candidates=CellNeighbors::new(p,hmax).all_candidates();
    for i in 0..n {
        let h=state.smoothing_m[i];
        let mut sum=V3::default();
        let mut trace=0.0;
        for &j in &candidates[i] {
            if j==i {continue;}
            let disp=p[i].pos_m-p[j].pos_m;
            let r=disp.len();
            if r<=0.0||r>=2.0*h {continue;}
            let dwd_r=sph_kernel::radial_derivative(r,h,kernel);
            let grad=disp*(dwd_r/r);
            let volume=p[j].mass_kg/state.densities_kg_m3[j];
            sum+=grad*volume;
            trace+=volume*((p[j].pos_m-p[i].pos_m).dot(grad));
        }
        e0[i]=sum*h;
        traces[i]=trace/3.0;
    }
    let quarter=(n/4).max(1);
    let e0_rms=(e0.iter().map(|v|v.len_sq()).sum::<f64>()/n as f64).sqrt();
    let e0_inner_rms=(sorted.iter().take(quarter)
        .map(|&i|e0[i].len_sq()).sum::<f64>()/quarter as f64).sqrt();
    let e0_outer_rms=(sorted.iter().rev().take(quarter)
        .map(|&i|e0[i].len_sq()).sum::<f64>()/quarter as f64).sqrt();
    let mut outer_radial=0.0;
    let mut inner_tangent_sq=0.0;
    for &i in sorted.iter().rev().take(quarter) {
        let unit=(p[i].pos_m-center).normalized();
        outer_radial+=e0[i].dot(unit);
    }
    for &i in sorted.iter().take(quarter) {
        let unit=(p[i].pos_m-center).normalized();
        let t=e0[i]-unit*e0[i].dot(unit);
        inner_tangent_sq+=t.len_sq();
    }
    ConsistencyMetrics{
        e0_rms,e0_inner_rms,e0_outer_rms,
        e0_outer_mean_radial:outer_radial/quarter as f64,
        e0_inner_tangent_rms:(inner_tangent_sq/quarter as f64).sqrt(),
        e1_trace_mean:traces.iter().sum::<f64>()/n as f64,
        e1_trace_outer_mean:sorted.iter().rev().take(quarter)
            .map(|&i|traces[i]).sum::<f64>()/quarter as f64,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{profile_for_mass,DegenerateElectrons};
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_adaptive::{AdaptiveConfig,update_adaptive_densities};
    #[test]
    fn consistency_is_finite_and_invariant_under_translation() {
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let (mut p,_)=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),96,17);
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        let mut astate=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut astate,cfg);
        let a=audit(&p,&astate,cfg.kernel);
        assert!(a.e0_rms.is_finite()&&a.e0_rms>0.0);
        assert!(a.e0_outer_rms.is_finite());
        assert!(a.e1_trace_mean.is_finite());
        for x in &mut p{x.pos_m+=V3::new(-2.0e8,3.0e8,1.0e8);}
        let mut bstate=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut bstate,cfg);
        let b=audit(&p,&bstate,cfg.kernel);
        assert!((a.e0_rms-b.e0_rms).abs()<1e-10);
        assert!((a.e1_trace_mean-b.e1_trace_mean).abs()<1e-10);
    }
}
