//! BG-PHY-19: local polynomial reconstructions at the white-dwarf surface.
//!
//! A frozen-particle diagnostic comparing weighted least-squares linear
//! and quadratic fits to four scalar fields: SPH pressure, 1D-reference
//! pressure, SPH specific enthalpy, and 1D-reference specific enthalpy.
//!
//! For the barotropic cold electron EOS, a_pressure = -grad(H), where
//! H=int dP/rho. H remains smooth enough near the free surface that
//! its gradient can be better conditioned than -grad(P)/rho.
//!
//! No reconstructed gradient is used by live forces or dynamics.
use crate::physics::G_SI;
use crate::sph3::{Fluid3,V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_neighbors::CellNeighbors;
use crate::white_dwarf::WhiteDwarfProfile;

pub const FIELD_LABELS:[&str;8]=[
    "pressure_num_linear","pressure_num_quadratic",
    "pressure_ref_linear","pressure_ref_quadratic",
    "enthalpy_num_linear","enthalpy_num_quadratic",
    "enthalpy_ref_linear","enthalpy_ref_quadratic",
];
type Vec9=[f64;9];
type Mat9=[[f64;9];9];
type Rhs=[[f64;4];9];
/// Relative normalized errors use the same independent 1D gravity field.
#[derive(Clone,Copy,Debug)]
pub struct PolynomialMetrics {
    pub count:usize,
    pub usable_fraction:f64,
    pub vector_error:[f64;8],
    pub radial_ratio:[f64;8],
    pub tangent_rms:[f64;8],
    /// Squared Cholesky pivot ratio; a conditioning *proxy*.
    pub max_linear_pivot_spread:f64,
    pub max_quadratic_pivot_spread:f64,
    pub density_ln_bias_mean:f64,
}
#[derive(Clone,Debug)]
pub struct PolynomialAudit {
    pub global:PolynomialMetrics,
    pub inner:PolynomialMetrics,
    pub outer:PolynomialMetrics,
}
/// Displacement divided by local h. Quadratic coefficients encode the
/// Hessian using 0.5 for each diagonal monomial.
fn basis(q:V3)->Vec9 {
    [q.x,q.y,q.z,0.5*q.x*q.x,0.5*q.y*q.y,0.5*q.z*q.z,
        q.x*q.y,q.x*q.z,q.y*q.z]
}
/// Solve four weighted scalar-field regressions with a shared Cholesky
/// factorization of diagonally preconditioned normal equations.
/// Returns 9 coefficients per field (higher coefficients zero for linear),
/// plus a Cholesky pivot spread diagnostic.
fn fit(samples:&[(Vec9,f64,[f64;4])],n:usize)->Option<(Rhs,f64)> {
    assert!(n==3||n==9);
    if samples.len()<n+2 {return None;}
    let mut matrix:Mat9=[[0.0;9];9];
    let mut right:Rhs=[[0.0;4];9];
    for (q,w,fields) in samples {
        if *w<=0.0||!w.is_finite() {continue;}
        for i in 0..n {
            for j in 0..n {matrix[i][j]+=w*q[i]*q[j];}
            for k in 0..4 {right[i][k]+=w*q[i]*fields[k];}
        }
    }
    let mut scales=[0.0;9];
    for i in 0..n {
        if matrix[i][i]<=0.0||!matrix[i][i].is_finite() {return None;}
        scales[i]=matrix[i][i].sqrt();
    }
    let mut l:Mat9=[[0.0;9];9];
    let mut max_pivot=0.0_f64;
    let mut min_pivot=f64::INFINITY;
    for i in 0..n {
        for j in 0..=i {
            let mut a=matrix[i][j]/(scales[i]*scales[j]);
            for k in 0..j {a-=l[i][k]*l[j][k];}
            if i==j {
                if a<=1e-10||!a.is_finite() {return None;}
                max_pivot=max_pivot.max(a);
                min_pivot=min_pivot.min(a);
                l[i][i]=a.sqrt();
            } else {l[i][j]=a/l[j][j];}
        }
    }
    let spread=max_pivot/min_pivot;
    if !spread.is_finite()||spread>1e8{return None;}
    let mut result:Rhs=[[0.0;4];9];
    for col in 0..4 {
        let mut y=[0.0_f64;9];
        for i in 0..n {
            let mut v=right[i][col]/scales[i];
            for j in 0..i {v-=l[i][j]*y[j];}
            y[i]=v/l[i][i];
        }
        let mut x=[0.0_f64;9];
        for i in (0..n).rev() {
            let mut v=y[i];
            for j in i+1..n {v-=l[j][i]*x[j];}
            x[i]=v/l[i][i];
        }
        for i in 0..n{result[i][col]=x[i]/scales[i];}
    }
    Some((result,spread))
}
fn acceleration(beta:&Rhs,field:usize,h:f64,rho:f64)->V3 {
    let divisor=if field<2 {h*rho}else{h};
    V3::new(-beta[0][field]/divisor,
            -beta[1][field]/divisor,
            -beta[2][field]/divisor)
}
#[derive(Clone,Copy)]
struct Sample{
    mass:f64,
    gravity:f64,
    radial:V3,
    predictions:[V3;8],
    linear_spread:f64,
    quadratic_spread:f64,
    density_ln_bias:f64,
    valid:bool,
}
#[derive(Default)]
struct Sums {
    count:usize,
    valid:usize,
    m:f64,
    reference2:f64,
    reference_linear:f64,
    error:[f64;8],
    radial:[f64;8],
    tangent:[f64;8],
    linear_max:f64,
    quadratic_max:f64,
    density_ln_sum:f64,
}
impl Sums {
    fn add(&mut self,p:Sample){
        self.count+=1;
        if !p.valid{return;}
        self.valid+=1;
        self.m+=p.mass;
        self.reference2+=p.mass*p.gravity*p.gravity;
        self.reference_linear+=p.mass*p.gravity;
        self.linear_max=self.linear_max.max(p.linear_spread);
        self.quadratic_max=self.quadratic_max.max(p.quadratic_spread);
        self.density_ln_sum+=p.mass*p.density_ln_bias;
        for (k,a) in p.predictions.iter().enumerate(){
            self.error[k]+=p.mass*(*a-p.radial*p.gravity).len_sq();
            self.radial[k]+=p.mass*a.dot(p.radial);
            self.tangent[k]+=p.mass*
                (*a-p.radial*a.dot(p.radial)).len_sq();
        }
    }
    fn finish(&self)->PolynomialMetrics{
        let divisor=self.reference2.max(1e-50);
        PolynomialMetrics{
            count:self.count,
            usable_fraction:self.valid as f64/self.count.max(1) as f64,
            vector_error:self.error.map(|s|(s/divisor).sqrt()),
            radial_ratio:self.radial.map(|v|v/self.reference_linear.max(1e-50)),
            tangent_rms:self.tangent.map(|s|(s/divisor).sqrt()),
            max_linear_pivot_spread:self.linear_max,
            max_quadratic_pivot_spread:self.quadratic_max,
            density_ln_bias_mean:self.density_ln_sum/self.m.max(1e-50),
        }
    }
}
/// Weighted polynomial fits over SPH-supported neighbors, at frozen
/// positions with exact-root densities already in state.
pub fn audit(
    p:&[Fluid3],state:&AdaptiveState,kernel:DensityKernel,
    star:&WhiteDwarfProfile,
)->PolynomialAudit {
    let n=p.len();
    assert!(n>=16&&state.valid_for(n));
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),
        |s,x|s+x.pos_m*x.mass_kg)/mass;
    let mut radius=vec![0.0;n];
    let mut fields=vec![[0.0_f64;4];n];
    for (i,x) in p.iter().enumerate(){
        let r=(x.pos_m-center).len();
        assert!(r>0.0&&r<star.radius_m);
        radius[i]=r;
        let rho_num=state.densities_kg_m3[i];
        let rho_ref=star.density_at(r);
        assert!(rho_num>0.0&&rho_ref>0.0);
        fields[i]=[
            star.eos.pressure_pa(rho_num),
            star.eos.pressure_pa(rho_ref),
            star.eos.enthalpy_m2_s2(rho_num),
            star.eos.enthalpy_m2_s2(rho_ref),
        ];
    }
    let maxh=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let candidates=CellNeighbors::new(p,maxh).all_candidates();
    let mut samples=vec![None;n];
    for i in 0..n {
        let h=state.smoothing_m[i];
        let mut neighbors=Vec::new();
        for &j in &candidates[i] {
            if i==j {continue;}
            let d=p[j].pos_m-p[i].pos_m;
            let r=d.len();
            if r>=2.0*h||r==0.0{continue;}
            let w=p[j].mass_kg/state.densities_kg_m3[j]*
                sph_kernel::kernel_and_dh(r,h,kernel).0;
            let diff=std::array::from_fn(|k|fields[j][k]-fields[i][k]);
            neighbors.push((basis(d*(1.0/h)),w,diff));
        }
        let linear=fit(&neighbors,3);
        let quadratic=fit(&neighbors,9);
        let valid=linear.is_some()&&quadratic.is_some();
        let mut predictions=[V3::default();8];
        let mut spreads=(0.0,0.0);
        if let (Some((a,sa)),Some((b,sb)))=(linear,quadratic) {
            spreads=(sa,sb);
            let rho_num=state.densities_kg_m3[i];
            let rho_ref=star.density_at(radius[i]);
            for k in 0..4{
                let denom=if k%2==0{rho_num}else{rho_ref};
                predictions[2*k]=acceleration(&a,k,h,denom);
                predictions[2*k+1]=acceleration(&b,k,h,denom);
            }
        }
        let radial=(p[i].pos_m-center).normalized();
        samples[i]=Some(Sample{
            mass:p[i].mass_kg,
            gravity:G_SI*star.enclosed_mass_at(radius[i])/radius[i].powi(2),
            radial,predictions,
            linear_spread:spreads.0,quadratic_spread:spreads.1,
            density_ln_bias:(state.densities_kg_m3[i]/
                star.density_at(radius[i])).ln(),valid,
        });
    }
    let mut rank=(0..n).collect::<Vec<_>>();
    rank.sort_by(|&i,&j|radius[i].total_cmp(&radius[j]).then(i.cmp(&j)));
    let quarter=(n/4).max(1);
    let mut global=Sums::default();
    let mut inner=Sums::default();
    let mut outer=Sums::default();
    for (k,&i) in rank.iter().enumerate(){
        let sample=samples[i].expect("particle audit missing");
        global.add(sample);
        if k<quarter {inner.add(sample);}
        if k>=n-quarter {outer.add(sample);}
    }
    PolynomialAudit{
        global:global.finish(),inner:inner.finish(),outer:outer.finish(),
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_adaptive::{AdaptiveConfig,update_adaptive_densities};
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    #[test]
    fn reproduces_quadratic_and_linear_fields_with_irregular_neighbors(){
        let coefficients=[1.25,-0.7,0.9,0.6,-0.3,0.8,0.2,0.4,-0.5];
        let mut s=Vec::new();
        for x in [-0.8,-0.3,0.2,0.7]{
            for y in [-0.7,0.1,0.8] {
                for z in [-0.6,0.25,0.9]{
                    let p=basis(V3::new(x,y,z));
                    let delta=p.iter().zip(coefficients)
                        .map(|(a,b)|a*b).sum::<f64>();
                    s.push((p,0.4+(x*x+y*y+z*z).exp().recip(),
                        [delta,2.0*delta,-0.5*delta,3.0*delta]));
                }
            }
        }
        let (quad,qcond)=fit(&s,9).unwrap();
        let (lin,_) =fit(&s,3).unwrap();
        assert!(qcond.is_finite()&&qcond<1e8);
        for i in 0..9{
            assert!((quad[i][0]-coefficients[i]).abs()<1e-9);
            assert!((quad[i][1]-2.0*coefficients[i]).abs()<1e-9);
        }
        assert!((quad[0][0]-lin[0][0]).abs()>1e-3);
        let constant=s.iter().map(|(q,w,_)|(*q,*w,[0.0;4]))
            .collect::<Vec<_>>();
        let (flat,_)=fit(&constant,9).unwrap();
        assert!(flat.iter().all(|row|row.iter().all(|v|v.abs()<1e-10)));
    }
    #[test]
    fn real_star_fits_are_finite_and_translation_invariant(){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let (mut p,_)=seed_cubic_shells_oriented(
            &star,V3::new(1e8,0.0,0.0),V3::default(),96,17);
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        let mut s0=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut s0,cfg);
        let a=audit(&p,&s0,cfg.kernel,&star);
        assert!(a.global.usable_fraction>0.75);
        assert!(a.global.vector_error.iter().all(|x|x.is_finite()));
        for x in &mut p{x.pos_m+=V3::new(-2e8,3e8,5e8);}
        let mut s1=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut s1,cfg);
        let b=audit(&p,&s1,cfg.kernel,&star);
        assert!((a.global.vector_error[7]-b.global.vector_error[7]).abs()<1e-8);
    }
}
