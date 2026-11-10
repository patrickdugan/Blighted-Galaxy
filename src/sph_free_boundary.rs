//! BG-PHY-25: free-surface enthalpy reconstructions with geometry-only
//! one-sided vacuum support. Frozen scalar diagnostics only.
//!
//! A smooth barotropic enthalpy field is locally affine near the vacuum
//! surface even when the pressure is strongly curved: grad H = -g.
//! The local exterior normal is inferred from the missing first moment
//! Q_i = sum_j (m_j/rho_j) W_ij (x_j-x_i)/h_i; -Q points outward.
//! The distance to a locally planar support boundary is inferred by
//! inverting the normalized half-space mass fraction of the chosen kernel.
//! A missing-support ghost datum H=0 at that plane can then enter a
//! weighted affine regression. No hydrostatic reference labels are used.
//!
//! We retain BG-PHY-20's smooth partition exposure for the *diagnostic*
//! blended scalar; full energy derivatives must be derived and verified
//! separately before this scalar may feed live gameplay dynamics.
use crate::sph3::{Fluid3,V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self,DensityKernel};
use crate::sph_neighbors::CellNeighbors;
use crate::white_dwarf::DegenerateElectrons;
use std::f64::consts::PI;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum EnthalpyCandidate {
    MlsH,
    BlendedH,
    GhostH,
    GhostBlendedH,
    /// Compensate common volume bias via |Q|/S before placing H=0 ghost.
    MomentRatioGhostBlend,
}
impl EnthalpyCandidate {
    pub fn label(self)->&'static str {
        match self{
            Self::MlsH=>"enthalpy_affine",
            Self::BlendedH=>"enthalpy_affine_blend",
            Self::GhostH=>"vacuum_ghost_enthalpy",
            Self::GhostBlendedH=>"vacuum_ghost_blend",
            Self::MomentRatioGhostBlend=>"moment_ratio_ghost_blend",
        }
    }
}
#[derive(Clone,Debug)]
pub struct EnthalpyDensity {
    pub mode:EnthalpyCandidate,
    pub rho:Vec<f64>,
    pub fallback_fraction:f64,
    pub maximum_pivot_spread:f64,
}
#[derive(Clone,Copy,Debug)]
pub struct BoundaryGeometry{
    /// Raw volume partition sum sum_j V_j W_ij.
    pub support:f64,
    /// Dimensionless missing first moment.
    pub first_moment:V3,
    /// Outward normal inferred from the first moment alone.
    pub outward:V3,
    /// Estimated distance to the locally planar vacuum in units h.
    pub plane_distance_h:f64,
    pub ghost_weight:f64,
    pub first_to_zeroth_ratio:f64,
    pub ratio_plane_distance_h:f64,
    pub ratio_ghost_weight:f64,
    pub blend_weight:f64,
}
#[derive(Clone,Debug)]
pub struct BoundaryBatch{
    pub fields:Vec<EnthalpyDensity>,
    pub geometry:Vec<BoundaryGeometry>,
}
#[derive(Clone,Copy,Debug)]
struct Sample{
    q:V3,
    w:f64,
    h:f64,
}
fn basis(q:V3)->[f64;4]{[1.0,q.x,q.y,q.z]}
fn weighted_affine(samples:&[Sample],ghost:Option<(V3,f64)>)
    ->Option<([f64;4],f64)> {
    if samples.len()<6{return None;}
    let mut a=[[0.0_f64;4];4];
    let mut b=[0.0_f64;4];
    for s in samples {
        let v=basis(s.q);
        for i in 0..4{
            b[i]+=s.w*v[i]*s.h;
            for j in 0..4{a[i][j]+=s.w*v[i]*v[j];}
        }
    }
    if let Some((q,w))=ghost {
        let v=basis(q);
        for i in 0..4{
            for j in 0..4{a[i][j]+=w*v[i]*v[j];}
        }
    }
    let mut scale=[0.0;4];
    for i in 0..4{
        if a[i][i]<=0.0||!a[i][i].is_finite(){return None;}
        scale[i]=a[i][i].sqrt();
    }
    let mut l=[[0.0_f64;4];4];
    let (mut min_pivot,mut max_pivot)=(f64::INFINITY,0.0_f64);
    for i in 0..4{
        for j in 0..=i {
            let mut value=a[i][j]/(scale[i]*scale[j]);
            for k in 0..j {value-=l[i][k]*l[j][k];}
            if i==j {
                if !value.is_finite()||value<=1e-10{return None;}
                min_pivot=min_pivot.min(value);
                max_pivot=max_pivot.max(value);
                l[i][i]=value.sqrt();
            }else{l[i][j]=value/l[j][j];}
        }
    }
    let spread=max_pivot/min_pivot;
    if !spread.is_finite()||spread>1e8{return None;}
    let mut y=[0.0;4];
    for i in 0..4{
        let mut value=b[i]/scale[i];
        for j in 0..i{value-=l[i][j]*y[j];}
        y[i]=value/l[i][i];
    }
    let mut beta=[0.0;4];
    for i in (0..4).rev(){
        let mut value=y[i];
        for j in i+1..4{value-=l[j][i]*beta[j];}
        beta[i]=value/l[i][i];
    }
    for i in 0..4{beta[i]/=scale[i];}
    if !beta.iter().all(|v|v.is_finite()){return None;}
    Some((beta,spread))
}
fn blend_weight(s:f64)->f64{
    let t=(2.0*(1.0-s)).clamp(0.0,1.0);
    t*t*(3.0-2.0*t)
}
/// Geometric 3D compact-support mass inside a planar half-space.
/// The point is a distance a*h inside a flat boundary (a in [0,2]).
fn planar_support_fraction(a:f64,k:DensityKernel)->f64 {
    assert!((0.0..=2.0).contains(&a));
    let n=128;
    let step=2.0/n as f64;
    let mut mass=0.0;
    for i in 0..n {
        let q=(i as f64+0.5)*step;
        let angle=if q<=a {1.0}else{0.5*(1.0+a/q)};
        let w=sph_kernel::kernel_and_dh(q,1.0,k).0;
        mass+=4.0*PI*q*q*w*angle*step;
    }
    mass
}
fn plane_lookup(k:DensityKernel)->[f64;65]{
    let mut out=[0.0;65];
    let norm=planar_support_fraction(2.0,k);
    for i in 0..=64{
        out[i]=planar_support_fraction(2.0*i as f64/64.0,k)/norm;
    }
    out
}
/// For q>a, the fluid-side angular first-moment coefficient
/// is [(a/q)^2-1]/4. Return its positive magnitude.
fn planar_first_moment(a:f64,k:DensityKernel)->f64{
    assert!((0.0..=2.0).contains(&a));
    let n=128;
    let step=2.0/n as f64;
    let mut magnitude=0.0;
    for i in 0..n{
        let q=(i as f64+0.5)*step;
        if q<=a {continue;}
        let w=sph_kernel::kernel_and_dh(q,1.0,k).0;
        magnitude+=PI*q.powi(3)*w*(1.0-(a/q).powi(2))*step;
    }
    magnitude
}
fn ratio_lookup(k:DensityKernel)->[f64;65]{
    let mut out=[0.0;65];
    for i in 0..=64{
        let a=2.0*i as f64/64.0;
        out[i]=planar_first_moment(a,k)/
            planar_support_fraction(a,k).max(1e-12);
    }
    out
}
fn invert_moment_ratio(r:f64,lut:&[f64;65])->f64{
    if r>=lut[0]{return 0.0;}
    if r<=0.0{return 2.0;}
    for i in 1..=64{
        if lut[i]<=r{
            let denom=(lut[i-1]-lut[i]).max(1e-12);
            return 2.0/64.0*((i-1) as f64+
                ((lut[i-1]-r)/denom).clamp(0.0,1.0));
        }
    }
    2.0
}
fn ideal_fraction_at_depth(a:f64,lut:&[f64;65])->f64{
    let x=(a*32.0).clamp(0.0,64.0);
    let i=(x.floor() as usize).min(63);
    lut[i]+(x-i as f64)*(lut[i+1]-lut[i])
}
fn inverse_planar_fraction(s:f64,lut:&[f64;65])->f64{
    if s<=0.5{return 0.0;}
    if s>=1.0{return 2.0;}
    for i in 1..=64 {
        if lut[i]>=s{
            let denom=(lut[i]-lut[i-1]).max(1e-12);
            return 2.0/64.0*((i-1) as f64+
                ((s-lut[i-1])/denom).clamp(0.0,1.0));
        }
    }
    2.0
}
pub fn reconstruct(
    p:&[Fluid3],state:&AdaptiveState,
    kernel:DensityKernel,eos:DegenerateElectrons,
)->BoundaryBatch {
    let n=p.len();
    assert!(n>=16&&state.valid_for(n));
    let rho=&state.densities_kg_m3;
    assert!(rho.iter().all(|x|x.is_finite()&&*x>0.0));
    let enthalpy=rho.iter().map(|x|eos.enthalpy_m2_s2(*x))
        .collect::<Vec<_>>();
    let lut=plane_lookup(kernel);
    let ratio_lut=ratio_lookup(kernel);
    let max_h=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let candidates=CellNeighbors::new(p,max_h).all_candidates();
    let mut geometry=Vec::with_capacity(n);
    let mut plain=Vec::with_capacity(n);
    let mut blended=Vec::with_capacity(n);
    let mut ghost_plain=Vec::with_capacity(n);
    let mut ghost_blended=Vec::with_capacity(n);
    let mut ratio_ghost_blended=Vec::with_capacity(n);
    let (mut failed_plain,mut failed_ghost,mut failed_ratio)=
        (0usize,0usize,0usize);
    let (mut worst_plain,mut worst_ghost,mut worst_ratio)=
        (0.0_f64,0.0_f64,0.0_f64);
    for i in 0..n{
        let h=state.smoothing_m[i];
        let mut samples=Vec::new();
        let mut partition=0.0;
        let mut first=V3::default();
        for &j in &candidates[i]{
            let d=p[j].pos_m-p[i].pos_m;
            let r=d.len();
            if r>=2.0*h{continue;}
            let w=p[j].mass_kg/rho[j]*
                sph_kernel::kernel_and_dh(r,h,kernel).0;
            if w<=0.0 {continue;}
            let q=d*(1.0/h);
            first+=q*w;
            partition+=w;
            samples.push(Sample{q,w,h:enthalpy[j]});
        }
        assert!(partition>0.0);
        let outward=first*(-1.0/first.len().max(1e-30));
        let distance=inverse_planar_fraction(partition,&lut);
        let ghost_weight=(1.0-partition).clamp(0.0,0.5);
        let moment_ratio=first.len()/partition;
        let ratio_depth=invert_moment_ratio(moment_ratio,&ratio_lut);
        // Use the ideal missing *geometric* volume at this depth,
        // so an overall biased source volume does not over-weight
        // the vacuum ghost in the affine normal equations.
        let ratio_ghost_weight=
            (1.0-ideal_fraction_at_depth(ratio_depth,&lut))
                .clamp(0.0,0.5);
        let alpha=blend_weight(partition);
        geometry.push(BoundaryGeometry{
            support:partition,first_moment:first,outward,
            plane_distance_h:distance,ghost_weight,
            first_to_zeroth_ratio:moment_ratio,
            ratio_plane_distance_h:ratio_depth,
            ratio_ghost_weight,
            blend_weight:alpha,
        });
        let ordinary=weighted_affine(&samples,None);
        let with_ghost=if ghost_weight>0.0&&first.len()>1e-8{
            weighted_affine(&samples,Some((outward*distance,ghost_weight)))
        }else{ordinary};
        let with_ratio_ghost=if ratio_ghost_weight>0.0&&first.len()>1e-8{
            weighted_affine(&samples,Some((
                outward*ratio_depth,ratio_ghost_weight)))
        }else{ordinary};
        let to_rho=|fit:Option<([f64;4],f64)>|
            fit.filter(|(beta,_)|beta[0]>0.0)
                .map(|(beta,spread)|(eos.density_from_enthalpy(beta[0]),spread))
                .filter(|(value,_)|value.is_finite()&&*value>0.0);
        let hplain=to_rho(ordinary);
        let hghost=to_rho(with_ghost);
        let hratio=to_rho(with_ratio_ghost);
        if hplain.is_none(){failed_plain+=1;}
        if hghost.is_none(){failed_ghost+=1;}
        if hratio.is_none(){failed_ratio+=1;}
        if let Some((_,c))=hplain{worst_plain=worst_plain.max(c);}
        if let Some((_,c))=hghost{worst_ghost=worst_ghost.max(c);}
        if let Some((_,c))=hratio{worst_ratio=worst_ratio.max(c);}
        let pvalue=hplain.map(|x|x.0).unwrap_or(rho[i]);
        let gvalue=hghost.map(|x|x.0).unwrap_or(rho[i]);
        let rvalue=hratio.map(|x|x.0).unwrap_or(rho[i]);
        plain.push(pvalue);
        ghost_plain.push(gvalue);
        blended.push((1.0-alpha)*rho[i]+alpha*pvalue);
        ghost_blended.push((1.0-alpha)*rho[i]+alpha*gvalue);
        ratio_ghost_blended.push((1.0-alpha)*rho[i]+alpha*rvalue);
    }
    BoundaryBatch{
        fields:vec![
            EnthalpyDensity{
                mode:EnthalpyCandidate::MlsH,rho:plain,
                fallback_fraction:failed_plain as f64/n as f64,
                maximum_pivot_spread:worst_plain,
            },
            EnthalpyDensity{
                mode:EnthalpyCandidate::BlendedH,rho:blended,
                fallback_fraction:failed_plain as f64/n as f64,
                maximum_pivot_spread:worst_plain,
            },
            EnthalpyDensity{
                mode:EnthalpyCandidate::GhostH,rho:ghost_plain,
                fallback_fraction:failed_ghost as f64/n as f64,
                maximum_pivot_spread:worst_ghost,
            },
            EnthalpyDensity{
                mode:EnthalpyCandidate::GhostBlendedH,rho:ghost_blended,
                fallback_fraction:failed_ghost as f64/n as f64,
                maximum_pivot_spread:worst_ghost,
            },
            EnthalpyDensity{
                mode:EnthalpyCandidate::MomentRatioGhostBlend,
                rho:ratio_ghost_blended,
                fallback_fraction:failed_ratio as f64/n as f64,
                maximum_pivot_spread:worst_ratio,
            },
        ],
        geometry,
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::profile_for_mass;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::sph_adaptive::{AdaptiveConfig,update_adaptive_densities};
    #[test]
    fn planar_fraction_half_at_surface_and_unity_inside(){
        for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let lut=plane_lookup(kernel);
            assert!((lut[0]-0.5).abs()<2e-12);
            assert!((lut[64]-1.0).abs()<2e-12);
            for i in 1..=64{assert!(lut[i]>=lut[i-1]);}
            for d in [0.2,0.5,1.0,1.5]{
                let s=planar_support_fraction(d,kernel)/
                    planar_support_fraction(2.0,kernel);
                let recovered=inverse_planar_fraction(s,&lut);
                assert!((d-recovered).abs()<0.005,
                    "kernel {:?}, distance {d}, recovered {recovered}",kernel);
            }
        }
    }
    #[test]
    fn moment_ratio_recovers_boundary_depth_under_volume_bias(){
        for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let ratio=ratio_lookup(kernel);
            let fraction=plane_lookup(kernel);
            assert!(ratio[0]>0.0&&ratio[64]==0.0);
            for i in 1..=64{
                assert!(ratio[i]<=ratio[i-1]+1e-13);
            }
            for depth in [0.1,0.2,0.5,0.8,1.2,1.6]{
                let expected_s=planar_support_fraction(depth,kernel);
                let expected_q=planar_first_moment(depth,kernel);
                for bias in [0.63,0.81,1.0]{
                    // V_j -> bias*V_j multiplies both S and Q
                    // by the same factor, and cancels in |Q|/S.
                    let found=invert_moment_ratio(
                        bias*expected_q/(bias*expected_s),&ratio);
                    assert!((found-depth).abs()<0.006,
                        "kernel {:?}, depth {depth}, found {found}",
                        kernel);
                    let inferred_s=ideal_fraction_at_depth(found,&fraction);
                    assert!((inferred_s-expected_s).abs()<0.005);
                }
            }
        }
    }
    #[test]
    fn nonuniform_volume_bias_confounds_the_first_to_zeroth_ratio(){
        for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let depth=0.6;
            let normal_ratio=planar_first_moment(depth,kernel)/
                planar_support_fraction(depth,kernel);
            // A local source-volume error 1+b*q_z creates a new
            // first moment Q'=-Q+b*M_zz, so Q'/S' no longer
            // cancels the volume error. The flat half-space
            // second normal moment follows the exact angular integral.
            let n=128;
            let step=2.0/n as f64;
            let mut mzz=0.0;
            for i in 0..n {
                let q=(i as f64+0.5)*step;
                let fraction=if q<=depth {
                    1.0/3.0
                }else{
                    (1.0+(depth/q).powi(3))/6.0
                };
                mzz+=4.0*PI*q.powi(4)*
                    sph_kernel::kernel_and_dh(q,1.0,kernel).0*
                    fraction*step;
            }
            let b=-0.3;
            let s=planar_support_fraction(depth,kernel);
            let q=planar_first_moment(depth,kernel);
            let tilted_ratio=(q-b*mzz)/(s-b*q);
            let inferred=invert_moment_ratio(tilted_ratio,&ratio_lookup(kernel));
            assert!((inferred-depth).abs()>0.05,
                "nonuniform volume error should bias boundary depth");
            assert!((normal_ratio-q/s).abs()<1e-14);
        }
    }
    #[test]
    fn synthetic_one_sided_affine_vacuum_ghost_reproduces_linear_h(){
        let mut samples=Vec::new();
        // z <= 0 is fluid, z=0 is known vacuum H=0.
        // At a query point z=-0.6, H=10*(0-z).
        for x in [-1.0,-0.3,0.5,1.0]{
            for y in [-0.8,0.2,0.9]{
                for z in [-1.8,-1.1,-0.4,0.0]{
                    let dz=z+0.6;
                    samples.push(Sample{
                        q:V3::new(x,y,dz),w:(-x*x-y*y-dz*dz).exp(),
                        h:10.0*(-z),
                    });
                }
            }
        }
        let (fit,_)=weighted_affine(&samples,Some((V3::new(0.0,0.0,0.6),1.0))).unwrap();
        assert!((fit[0]-6.0).abs()<1e-10);
        assert!((fit[3]+10.0).abs()<1e-10);
        let (unconstrained,_)=weighted_affine(&samples,None).unwrap();
        assert!((unconstrained[0]-fit[0]).abs()<1e-10);
        // A constant field is reproduced exactly by ordinary MLS.
        for s in &mut samples{s.h=13.0;}
        let (flat,_)=weighted_affine(&samples,None).unwrap();
        assert!((flat[0]-13.0).abs()<1e-11);
        assert!(flat[1..4].iter().all(|x|x.abs()<1e-10));
    }
    #[test]
    fn half_space_normal_and_curved_enthalpy_are_finite(){
        let mut samples=Vec::new();
        let mut qmoment=V3::default();
        let mut total=0.0;
        for x in [-1.0,-0.4,0.4,1.0] {
            for y in [-0.8,0.1,0.9] {
                for z in [-1.6,-1.0,-0.4,0.0] {
                    let q=V3::new(x,y,z+0.6);
                    let depth=-z;
                    // Physical vacuum H(z=0)=0 with finite near-surface
                    // curvature. The affine fit is a first-order model.
                    let h=10.0*depth+3.0*depth*depth;
                    let w=(-q.len_sq()).exp();
                    qmoment+=q*w;
                    total+=w;
                    samples.push(Sample{q,w,h});
                }
            }
        }
        let outward=qmoment*(-1.0/qmoment.len());
        assert!(outward.z>0.3,
            "the geometry-only outward normal should face the vacuum");
        let (regular,_)=weighted_affine(&samples,None).unwrap();
        let (ghost,_)=weighted_affine(
            &samples,Some((outward*0.6,0.4*total))).unwrap();
        let true_center=10.0*0.6+3.0*0.6*0.6;
        // For finite curvature the two extrapolants may differ.
        // This synthetic test requires physical finite positive
        // enthalpy; the white dwarf experiment measures error.
        assert!(regular[0].is_finite()&&regular[0]>0.0);
        assert!(ghost[0].is_finite()&&ghost[0]>0.0);
        assert!((ghost[0]-true_center).abs()<true_center);
        assert!((regular[0]-true_center).abs()<true_center);
    }
    #[test]
    fn real_star_surface_candidates_positive_and_finite(){
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(&star,V3::new(1e8,0.0,0.0),
            V3::default(),96,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4().with_newton_density();
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        update_adaptive_densities(&p,&mut state,cfg);
        let batch=reconstruct(&p,&state,cfg.kernel,star.eos);
        assert_eq!(batch.fields.len(),5);
        for field in &batch.fields{
            assert_eq!(field.rho.len(),p.len());
            assert!(field.rho.iter().all(|v|v.is_finite()&&*v>0.0));
        }
        assert!(batch.geometry.iter().all(|g|
            g.support.is_finite()&&g.outward.finite()&&
            (0.0..=2.0).contains(&g.plane_distance_h)));
    }
}
