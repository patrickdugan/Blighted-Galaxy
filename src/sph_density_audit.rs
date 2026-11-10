//! BG-PHY-20: frozen SPH density estimators at a stellar free surface.
//!
//! All methods consume only particle positions, masses, the *existing* SPH
//! density snapshot, and the selected density kernel. Hydrostatic rho_ref
//! never enters an estimator; it is reserved for out-of-sample scoring.
//!
//! The independent free-boundary diagnostics are:
//! - original rho_i = sum_j m_j W_ij (control);
//! - Shepard partition-of-unity density: sum mW / sum (m/rho)W;
//! - local affine and quadratic weighted least-squares intercepts fitted
//!   to the existing SPH density samples;
//! - spherical enclosed-mass finite-volume secant (diagnostic only).
//!
//! All fields are FROZEN diagnostics. The live, conservative SPH density
//! closure, h, Omega, gravity, and momentum/energy forces stay unchanged.
use crate::sph3::{Fluid3,V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_neighbors::CellNeighbors;
use std::f64::consts::PI;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum DensityEstimator {
    Original,
    Shepard,
    MlsLinear,
    MlsQuadratic,
    /// Weight affine MLS by the missing kernel partition, preserving the
    /// original density in fully supported interior regions.
    SurfaceBlendedMls,
    RadialMassSecant,
}
impl DensityEstimator {
    pub fn name(self)->&'static str {
        match self {
            Self::Original=>"original",
            Self::Shepard=>"shepard",
            Self::MlsLinear=>"mls_linear",
            Self::MlsQuadratic=>"mls_quadratic",
            Self::SurfaceBlendedMls=>"surface_blended_mls",
            Self::RadialMassSecant=>"radial_mass_secant",
        }
    }
}
#[derive(Clone,Debug)]
pub struct DensityField {
    pub method:DensityEstimator,
    pub rho:Vec<f64>,
    /// When a fit is singular or predicts nonpositive rho, use baseline
    /// density explicitly and report the fallback fraction.
    pub fallback_fraction:f64,
    /// Maximum scaled-normal-matrix Cholesky pivot spread; a condition proxy.
    pub maximum_pivot_spread:f64,
}
#[derive(Clone,Debug)]
pub struct DensityBatch {
    pub fields:Vec<DensityField>,
    /// Discrete kernel partition unity sum_j (m_j/rho_SPH,j) W_ij.
    pub partition:Vec<f64>,
}
type Matrix=[[f64;10];10];
fn basis(q:V3)->[f64;10]{
    [1.0,q.x,q.y,q.z,0.5*q.x*q.x,0.5*q.y*q.y,
     0.5*q.z*q.z,q.x*q.y,q.x*q.z,q.y*q.z]
}
/// Return fitted local intercept and pivot-spread diagnostics.
/// The 4-component fit reproduces affine fields; the 10-component fit
/// reproduces quadratic fields at arbitrary irregular supported points.
fn weighted_intercept(samples:&[([f64;10],f64,f64)],n:usize)
    ->Option<(f64,f64)> {
    assert!(n==4||n==10);
    if samples.len()<n+2{return None;}
    let mut a:Matrix=[[0.0;10];10];
    let mut b=[0.0_f64;10];
    for (x,w,y) in samples {
        if !w.is_finite()||*w<=0.0 {continue;}
        for i in 0..n {
            b[i]+=w*x[i]*y;
            for j in 0..n{a[i][j]+=w*x[i]*x[j];}
        }
    }
    let mut scale=[0.0;10];
    for i in 0..n {
        if !a[i][i].is_finite()||a[i][i]<=0.0{return None;}
        scale[i]=a[i][i].sqrt();
    }
    let mut l:Matrix=[[0.0;10];10];
    let mut pivot_min=f64::INFINITY;
    let mut pivot_max=0.0_f64;
    for i in 0..n {
        for j in 0..=i {
            let mut value=a[i][j]/(scale[i]*scale[j]);
            for k in 0..j {value-=l[i][k]*l[j][k];}
            if i==j {
                if !value.is_finite()||value<=1e-10{return None;}
                pivot_min=pivot_min.min(value);
                pivot_max=pivot_max.max(value);
                l[i][i]=value.sqrt();
            } else {l[i][j]=value/l[j][j];}
        }
    }
    let spread=pivot_max/pivot_min;
    if !spread.is_finite()||spread>1e8{return None;}
    let mut z=[0.0;10];
    for i in 0..n{
        let mut y=b[i]/scale[i];
        for j in 0..i{y-=l[i][j]*z[j];}
        z[i]=y/l[i][i];
    }
    let mut x=[0.0;10];
    for i in (0..n).rev(){
        let mut value=z[i];
        for j in i+1..n{value-=l[j][i]*x[j];}
        x[i]=value/l[i][i];
    }
    let intercept=x[0]/scale[0];
    if intercept.is_finite()&&intercept>0.0{
        Some((intercept,spread))
    } else {None}
}
/// A spherical-symmetry diagnostic from ranked enclosed particle mass and
/// radial volume. Equal-mass particles approximate dM/dV via a local secant.
/// It is intentionally unsuitable for an aspherical tidal disruption.
fn radial_secant(p:&[Fluid3],center:V3,half_window:usize)->Vec<Option<f64>> {
    let n=p.len();
    let mut ranks=(0..n).map(|i|{
        let r=(p[i].pos_m-center).len();
        (r,i)
    }).collect::<Vec<_>>();
    ranks.sort_by(|a,b|a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut output=vec![None;n];
    for (index,&(_,i)) in ranks.iter().enumerate() {
        let mut lo=index.saturating_sub(half_window);
        let mut hi=(index+half_window).min(n-1);
        if hi-lo<2 {continue;}
        // Broaden a one-sided window at domain ends without injecting
        // an external surface boundary or reference mass-CDF.
        if lo==0 {hi=(2*half_window).min(n-1);}
        if hi==n-1 {lo=n.saturating_sub(2*half_window+1);}
        let (rlo,rhi)=(ranks[lo].0,ranks[hi].0);
        let volume=(4.0*PI/3.0)*(rhi.powi(3)-rlo.powi(3));
        let shell_mass=ranks[lo..hi].iter()
            .map(|x|p[x.1].mass_kg).sum::<f64>();
        let density=shell_mass/volume;
        if density.is_finite()&&density>0.0 {
            output[i]=Some(density);
        }
    }
    output
}
pub fn estimate(
    particles:&[Fluid3],state:&AdaptiveState,kernel:DensityKernel,
)->DensityBatch {
    let n=particles.len();
    assert!(n>=16&&state.valid_for(n));
    let baseline=&state.densities_kg_m3;
    assert!(baseline.iter().all(|x|x.is_finite()&&*x>0.0));
    let total_mass=particles.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=particles.iter().fold(V3::default(),|s,p|
        s+p.pos_m*p.mass_kg)/total_mass;
    let maxh=state.smoothing_m.iter().copied()
        .fold(0.0_f64,f64::max);
    let candidates=CellNeighbors::new(particles,maxh).all_candidates();
    let mut partition=vec![0.0;n];
    let mut shepard=baseline.to_vec();
    let mut affine=baseline.to_vec();
    let mut quadratic=baseline.to_vec();
    let mut failed1=0_usize;
    let mut failed2=0_usize;
    let mut max1=0.0_f64;
    let mut max2=0.0_f64;
    for i in 0..n {
        let h=state.smoothing_m[i];
        let mut samples=Vec::with_capacity(candidates[i].len());
        let mut sum_weight=0.0;
        let mut sum_rho_weight=0.0;
        for &j in &candidates[i] {
            let disp=particles[j].pos_m-particles[i].pos_m;
            let r=disp.len();
            if r>=2.0*h {continue;}
            let w=particles[j].mass_kg/baseline[j]*
                sph_kernel::kernel_and_dh(r,h,kernel).0;
            if w<=0.0 {continue;}
            sum_weight+=w;
            sum_rho_weight+=w*baseline[j];
            samples.push((basis(disp*(1.0/h)),w,baseline[j]));
        }
        assert!(sum_weight>0.0);
        partition[i]=sum_weight;
        shepard[i]=sum_rho_weight/sum_weight;
        if let Some((rho,spread))=weighted_intercept(&samples,4){
            affine[i]=rho;max1=max1.max(spread);
        } else {failed1+=1;}
        if let Some((rho,spread))=weighted_intercept(&samples,10){
            quadratic[i]=rho;max2=max2.max(spread);
        } else {failed2+=1;}
    }
    let half_window=((n as f64).sqrt()/2.0).round() as usize;
    // A half-space through a smooth free surface contributes roughly half
    // the compact kernel's volume. Interpolate the correction by the
    // normalized missing partition fraction, with a C1 smoothstep to
    // avoid a sharp transition in the diagnostic scalar field.
    let mut blended=baseline.to_vec();
    for i in 0..n {
        let x=((1.0-partition[i])/0.5).clamp(0.0,1.0);
        let weight=x*x*(3.0-2.0*x);
        blended[i]=baseline[i]*(1.0-weight)+affine[i]*weight;
    }
    let radial=radial_secant(particles,center,half_window.max(4));
    let failed_radial=radial.iter().filter(|x|x.is_none()).count();
    let radial_values=radial.into_iter().enumerate()
        .map(|(i,v)|v.unwrap_or(baseline[i])).collect::<Vec<_>>();
    DensityBatch{
        fields:vec![
            DensityField {
                method:DensityEstimator::Original,rho:baseline.to_vec(),
                fallback_fraction:0.0,maximum_pivot_spread:0.0,
            },
            DensityField {
                method:DensityEstimator::Shepard,rho:shepard,
                fallback_fraction:0.0,maximum_pivot_spread:0.0,
            },
            DensityField {
                method:DensityEstimator::MlsLinear,rho:affine,
                fallback_fraction:failed1 as f64/n as f64,
                maximum_pivot_spread:max1,
            },
            DensityField {
                method:DensityEstimator::MlsQuadratic,rho:quadratic,
                fallback_fraction:failed2 as f64/n as f64,
                maximum_pivot_spread:max2,
            },
            DensityField {
                method:DensityEstimator::SurfaceBlendedMls,rho:blended,
                fallback_fraction:failed1 as f64/n as f64,
                maximum_pivot_spread:max1,
            },
            DensityField {
                method:DensityEstimator::RadialMassSecant,rho:radial_values,
                fallback_fraction:failed_radial as f64/n as f64,
                maximum_pivot_spread:0.0,
            },
        ],
        partition,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_adaptive::{AdaptiveConfig,update_adaptive_densities};
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    #[test]
    fn mls_intercept_reproduces_affine_and_quadratic_fields() {
        let mut points=Vec::new();
        for x in [-0.9,-0.3,0.2,0.7]{
            for y in [-0.8,0.1,0.8]{
                for z in [-0.7,0.25,0.9]{
                    let q=V3::new(x,y,z);
                    let b=basis(q);
                    let linear=100.0+5.0*x-3.0*y+2.0*z;
                    let curved=linear+10.0*b[4]+8.0*b[5]+
                        6.0*b[6]-3.0*b[7]+2.0*b[8]+7.0*b[9];
                    points.push((b,0.5+(-x*x-y*y-z*z).exp(),curved));
                }
            }
        }
        let (q,_) =weighted_intercept(&points,10).unwrap();
        assert!((q-100.0).abs()<1e-9);
        let affine=points.iter().map(|(b,w,_)|
            (*b,*w,100.0+5.0*b[1]-3.0*b[2]+2.0*b[3]))
            .collect::<Vec<_>>();
        let (a,_) =weighted_intercept(&affine,4).unwrap();
        assert!((a-100.0).abs()<1e-10);
        let constant=points.iter().map(|(b,w,_)|(*b,*w,31.0))
            .collect::<Vec<_>>();
        let (c,_) =weighted_intercept(&constant,10).unwrap();
        assert!((c-31.0).abs()<1e-9);
    }
    #[test]
    fn seeded_fields_are_positive_and_have_exact_baseline(){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,V3::new(1e8,0.0,0.0),
            V3::default(),96,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut state,cfg);
        let report=estimate(&p,&state,cfg.kernel);
        assert_eq!(report.fields.len(),6);
        assert_eq!(report.fields[0].rho,state.densities_kg_m3);
        assert!(report.partition.iter().all(|v|v.is_finite()&&*v>0.0));
        for field in &report.fields {
            assert!(field.rho.iter().all(|rho|rho.is_finite()&&*rho>0.0));
            assert!((0.0..=1.0).contains(&field.fallback_fraction));
        }
    }
}
