//! BG-PHY-18: frozen pressure-gradient reconstruction with exact moments.
//!
//! Linear-gradient reconstruction uses the Bonet–Lok matrix
//! M_i = sum_j (m_j/rho_j) grad_i W_ij ⊗ (x_j-x_i), and
//! grad P = M_i^{-1} sum_j V_j(P_j-P_i) grad_i W_ij.
//! The direct difference gradient is also recorded. Both annihilate a
//! constant-pressure field. A second "oracle" field evaluates the SAME
//! discrete operators using the 1D hydrostatic P_ref(r), rho_ref(r).
//!
//! This module NEVER supplies a force to the dynamics. Replacing a
//! conservative symmetric SPH pair force by a directed corrected gradient
//! without a new variational formulation would alter conservation laws.
use crate::physics::G_SI;
use crate::sph3::{Fluid3,V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_neighbors::CellNeighbors;
use crate::white_dwarf::WhiteDwarfProfile;

type Mat3=[[f64;3];3];

fn vec3(v:V3)->[f64;3]{[v.x,v.y,v.z]}
fn from3(v:[f64;3])->V3{V3::new(v[0],v[1],v[2])}
fn multiply(a:&Mat3,v:V3)->V3{
    let p=vec3(v);
    from3(std::array::from_fn(|i|a[i].iter().zip(p).map(|(a,b)|a*b).sum()))
}
fn inverse(a:&Mat3)->Option<(Mat3,f64)>{
    let (x,y,z)=(a[0][0],a[0][1],a[0][2]);
    let (u,v,w)=(a[1][1],a[1][2],a[2][2]);
    let co=[
        [u*w-v*v, z*v-y*w, y*v-z*u],
        [z*v-y*w, x*w-z*z, y*z-x*v],
        [y*v-z*u, y*z-x*v, x*u-y*y]
    ];
    let det=x*co[0][0]+y*co[0][1]+z*co[0][2];
    let tr=(x+u+w)/3.0;
    if !det.is_finite()||!tr.is_finite()||tr<=0.0||
       det<=1e-12*tr.powi(3){return None;}
    let inv=co.map(|r|r.map(|v|v/det));
    let norm=|m:&Mat3|m.iter().map(|r|r.iter()
        .map(|x|x.abs()).sum::<f64>())
        .fold(0.0_f64,f64::max);
    let cond=norm(a)*norm(&inv);
    if !cond.is_finite()||cond>1e7{return None;}
    Some((inv,cond))
}
#[derive(Debug,Clone,Copy)]
pub struct GradientMetrics {
    pub count:usize,
    pub usable_fraction:f64,
    pub reference_rms_m_s2:f64,
    /// Error of -(sum V_j (P_j-P_i) grad W_ij) / rho_i.
    pub direct_num_error:f64,
    /// Same numerical pressure and density with first-moment inversion.
    pub corrected_num_error:f64,
    /// Replace P,rho at both positions by independently sampled 1D reference.
    pub direct_ref_error:f64,
    pub corrected_ref_error:f64,
    pub direct_num_radial_ratio:f64,
    pub corrected_num_radial_ratio:f64,
    pub direct_ref_radial_ratio:f64,
    pub corrected_ref_radial_ratio:f64,
    pub corrected_num_tangent_rms:f64,
    pub corrected_ref_tangent_rms:f64,
    pub density_ln_ratio_mean:f64,
    pub density_ln_ratio_rms:f64,
    pub moment_trace_mean:f64,
    pub moment_condition_max:f64,
}
#[derive(Debug,Clone)]
pub struct GradientReport {
    pub global:GradientMetrics,
    pub inner:GradientMetrics,
    pub outer:GradientMetrics,
}
#[derive(Clone,Copy,Debug)]
struct Sample {
    m:f64,
    g:f64,
    radial:V3,
    density_ln_ratio:f64,
    moment_trace:f64,
    moment_condition:Option<f64>,
    direct_num:V3,
    corrected_num:V3,
    direct_ref:V3,
    corrected_ref:V3,
}
#[derive(Default)]
struct Sums{
    total_count:usize,
    valid_count:usize,
    mass:f64,
    g2:f64,
    gmass:f64,
    err:[f64;4],
    radial:[f64;4],
    tangent_corrected_num:f64,
    tangent_corrected_ref:f64,
    log_ratio:f64,
    log_ratio2:f64,
    trace:f64,
    condition_max:f64,
}
impl Sums {
    fn add(&mut self,p:Sample){
        self.total_count+=1;
        if let Some(cond)=p.moment_condition{
            self.valid_count+=1;
            self.condition_max=self.condition_max.max(cond);
            self.mass+=p.m;
            self.g2+=p.m*p.g*p.g;
            self.gmass+=p.m*p.g;
            self.log_ratio+=p.m*p.density_ln_ratio;
            self.log_ratio2+=p.m*p.density_ln_ratio*p.density_ln_ratio;
            self.trace+=p.m*p.moment_trace;
            for (k,v) in [p.direct_num,p.corrected_num,
                          p.direct_ref,p.corrected_ref].iter().enumerate(){
                self.err[k]+=p.m*(*v-p.radial*p.g).len_sq();
                self.radial[k]+=p.m*v.dot(p.radial);
            }
            for (v,sum) in [
                (p.corrected_num,&mut self.tangent_corrected_num),
                (p.corrected_ref,&mut self.tangent_corrected_ref),
            ] {
                let t=v-p.radial*v.dot(p.radial);
                *sum+=p.m*t.len_sq();
            }
        }
    }
    fn finish(&self)->GradientMetrics {
        let scale=self.g2.max(1e-50);
        let radial_norm=self.gmass.max(1e-50);
        let m=self.mass.max(1e-50);
        GradientMetrics{
            count:self.total_count,
            usable_fraction:self.valid_count as f64/
                self.total_count.max(1) as f64,
            reference_rms_m_s2:(self.g2/m).sqrt(),
            direct_num_error:(self.err[0]/scale).sqrt(),
            corrected_num_error:(self.err[1]/scale).sqrt(),
            direct_ref_error:(self.err[2]/scale).sqrt(),
            corrected_ref_error:(self.err[3]/scale).sqrt(),
            direct_num_radial_ratio:self.radial[0]/radial_norm,
            corrected_num_radial_ratio:self.radial[1]/radial_norm,
            direct_ref_radial_ratio:self.radial[2]/radial_norm,
            corrected_ref_radial_ratio:self.radial[3]/radial_norm,
            corrected_num_tangent_rms:
                (self.tangent_corrected_num/scale).sqrt(),
            corrected_ref_tangent_rms:
                (self.tangent_corrected_ref/scale).sqrt(),
            density_ln_ratio_mean:self.log_ratio/m,
            density_ln_ratio_rms:(self.log_ratio2/m).sqrt(),
            moment_trace_mean:self.trace/m,
            moment_condition_max:self.condition_max,
        }
    }
}
/// Matrix-consistency correction audit at FROZEN particle positions.
/// Smoothing lengths and densities must already be initialized in state.
pub fn audit(
    p:&[Fluid3],state:&AdaptiveState,
    kernel:DensityKernel,star:&WhiteDwarfProfile
)->GradientReport {
    let n=p.len();
    assert!(n>=8&&state.valid_for(n));
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),
        |s,x|s+x.pos_m*x.mass_kg)/mass;
    let mut radii=vec![0.0_f64;n];
    let mut numerical_p=vec![0.0_f64;n];
    let mut ref_p=vec![0.0_f64;n];
    let mut ref_rho=vec![0.0_f64;n];
    let mut reference_g=vec![0.0_f64;n];
    for (i,x) in p.iter().enumerate(){
        let r=(x.pos_m-center).len();
        assert!(r>0.0&&r<star.radius_m,
            "BG-PHY-18 measures seeded interior particles only");
        radii[i]=r;
        ref_rho[i]=star.density_at(r);
        assert!(ref_rho[i]>0.0&&state.densities_kg_m3[i]>0.0);
        numerical_p[i]=star.eos.pressure_pa(state.densities_kg_m3[i]);
        ref_p[i]=star.eos.pressure_pa(ref_rho[i]);
        reference_g[i]=G_SI*star.enclosed_mass_at(r)/(r*r);
    }
    let mut rank=(0..n).collect::<Vec<_>>();
    rank.sort_by(|&i,&j|radii[i].total_cmp(&radii[j]).then(i.cmp(&j)));
    let mut samples=vec![None;n];
    let hmax=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let lists=CellNeighbors::new(p,hmax).all_candidates();
    for i in 0..n{
        let h=state.smoothing_m[i];
        let mut mat=[[0.0_f64;3];3];
        let mut dnum=V3::default();
        let mut dref=V3::default();
        for &j in &lists[i]{
            if i==j{continue;}
            let dr=p[j].pos_m-p[i].pos_m;
            let radius=dr.len();
            if radius<=0.0||radius>=2.0*h {continue;}
            let grad=(p[i].pos_m-p[j].pos_m)*
                (sph_kernel::radial_derivative(radius,h,kernel)/radius);
            let volume=p[j].mass_kg/state.densities_kg_m3[j];
            let u=vec3(grad);
            let d=vec3(dr);
            for a in 0..3{
                for b in 0..3{mat[a][b]+=volume*u[a]*d[b];}
            }
            dnum+=grad*(volume*(numerical_p[j]-numerical_p[i]));
            dref+=grad*(volume*(ref_p[j]-ref_p[i]));
        }
        let correction=inverse(&mat);
        let normal=(p[i].pos_m-center).normalized();
        let raw_num=dnum*(-1.0/state.densities_kg_m3[i]);
        let raw_ref=dref*(-1.0/ref_rho[i]);
        let (corrected_num,corrected_ref,condition)=
            if let Some((inv,cond))=correction{
                (multiply(&inv,dnum)*(-1.0/state.densities_kg_m3[i]),
                 multiply(&inv,dref)*(-1.0/ref_rho[i]),Some(cond))
            } else {(V3::default(),V3::default(),None)};
        samples[i]=Some(Sample{
            m:p[i].mass_kg,
            g:reference_g[i],radial:normal,
            density_ln_ratio:
                (state.densities_kg_m3[i]/ref_rho[i]).ln(),
            moment_trace:(mat[0][0]+mat[1][1]+mat[2][2])/3.0,
            moment_condition:condition,
            direct_num:raw_num,corrected_num,
            direct_ref:raw_ref,corrected_ref,
        });
    }
    let mut all=Sums::default();
    let mut inner=Sums::default();
    let mut outer=Sums::default();
    let quarter=(n/4).max(1);
    for (index,&i) in rank.iter().enumerate(){
        let s=samples[i].expect("every particle sampled");
        all.add(s);
        if index<quarter {inner.add(s);}
        if index>=n-quarter {outer.add(s);}
    }
    GradientReport{
        global:all.finish(),inner:inner.finish(),outer:outer.finish(),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_adaptive::{AdaptiveConfig,update_adaptive_densities};
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    #[test]
    fn inversion_reconstructs_linear_gradient_and_detects_singular() {
        let a=[[2.0,0.2,0.1],[0.2,1.8,0.3],[0.1,0.3,1.4]];
        let g=V3::new(7.0,-3.0,5.0);
        let b=multiply(&a,g);
        let (inverted,condition)=inverse(&a).unwrap();
        let found=multiply(&inverted,b);
        assert!((found-g).len()<1e-12);
        assert!(condition>1.0&&condition<10.0);
        assert!(inverse(&[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,0.0]]).is_none());
    }
    #[test]
    fn seeded_gradients_finite_and_translation_invariant() {
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let (mut p,_)=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),64,17);
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut state,cfg);
        let a=audit(&p,&state,cfg.kernel,&star);
        assert!(a.global.usable_fraction>0.8);
        assert!(a.global.corrected_num_error.is_finite());
        assert!(a.outer.corrected_ref_error.is_finite());
        for x in &mut p {x.pos_m+=V3::new(3.0e8,-4.0e8,5.0e8);}
        let mut shifted=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut shifted,cfg);
        let b=audit(&p,&shifted,cfg.kernel,&star);
        assert!((a.global.corrected_ref_error-
            b.global.corrected_ref_error).abs()<1e-8);
        assert!((a.outer.direct_ref_error-
            b.outer.direct_ref_error).abs()<1e-8);
    }
}
