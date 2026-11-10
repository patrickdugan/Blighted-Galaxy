//! BG-PHY-21: first-order 3D reproducing kernel with full QUERY gradients.
//!
//! At each particle i, q_ij=(x_i-x_j)/h_i, V_j=m_j/rho_SPH,j.
//! m0=sum V_j W_ij; m1=sum V_j q_ij W_ij; m2=sum V_j q_ij q_ij^T W_ij.
//! B_i=-m2^{-1} m1, A_i=1/(m0+B_i dot m1).
//! W^R_ij=A_i(1+B_i dot q_ij)W_ij.
//!
//! dW^R/dx_i includes dA/dx_i and dB/dx_i from all moment
//! derivatives. This is a *query-coordinate* derivative at fixed
//! particle source coordinates, volumes V_j and smoothing h_i. It is
//! checked against independently re-evaluated finite differences.
//!
//! Full Lagrangian particle-coordinate differentiation ALSO acts on
//! source coordinates, volumes, MLS surface factors and variable h;
//! this module does not claim to be that total derivative. No legacy
//! SPH behavior is modified.
use crate::sph3::{Fluid3,V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self,DensityKernel};

type Mat3=[[f64;3];3];
type MatDeriv=[[[f64;3];3];3];

fn as_arr(v:V3)->[f64;3]{[v.x,v.y,v.z]}
fn from_arr(a:[f64;3])->V3{V3::new(a[0],a[1],a[2])}
fn dot(a:[f64;3],b:[f64;3])->f64{
    a.iter().zip(b).map(|(x,y)|x*y).sum()
}
fn mat_vec(a:&Mat3,v:[f64;3])->[f64;3]{
    std::array::from_fn(|i|dot(a[i],v))
}
fn invert_spd(a:&Mat3)->Option<Mat3>{
    let (x,y,z)=(a[0][0],a[0][1],a[0][2]);
    let (u,v,w)=(a[1][1],a[1][2],a[2][2]);
    let cof=[
        [u*w-v*v,z*v-y*w,y*v-z*u],
        [z*v-y*w,x*w-z*z,y*z-x*v],
        [y*v-z*u,y*z-x*v,x*u-y*y],
    ];
    let det=dot(a[0],cof[0]);
    let tr=(x+u+w)/3.0;
    if !tr.is_finite()||tr<=0.0||!det.is_finite()||
        det<=1e-12*tr.powi(3){return None;}
    let inv=cof.map(|row|row.map(|q|q/det));
    let row_norm=|a:&Mat3|a.iter()
        .map(|row|row.iter().map(|v|v.abs()).sum::<f64>())
        .fold(0.0_f64,f64::max);
    let condition=row_norm(a)*row_norm(&inv);
    if !condition.is_finite()||condition>1e8{return None;}
    Some(inv)
}
#[derive(Clone,Copy,Debug)]
pub struct RkCoefficient {
    pub a:f64,
    /// Dimensionless linear reproduction coefficient in q=(xi-xj)/h.
    pub b:V3,
    pub grad_a:V3,
    /// grad_b[alpha] is the physical query gradient of b_alpha.
    pub grad_b:[V3;3],
    pub valid:bool,
    pub raw_partition:f64,
}
impl RkCoefficient {
    fn fallback(partition:f64)->Self {
        Self {a:1.0,b:V3::default(),grad_a:V3::default(),
            grad_b:[V3::default();3],valid:false,
            raw_partition:partition}
    }
}
/// Rebuild the correction for an arbitrary query coordinate. In this
/// operation all source particle coordinates and V_j are FIXED.
pub fn query_coefficient(
    query:V3,p:&[Fluid3],volumes:&[f64],h:f64,kernel:DensityKernel,
)->RkCoefficient {
    assert!(!p.is_empty()&&p.len()==volumes.len()&&h>0.0);
    let mut m0=0.0_f64;
    let mut m1=[0.0;3];
    let mut m2:Mat3=[[0.0;3];3];
    let mut dm0=[0.0;3];
    let mut dm1:Mat3=[[0.0;3];3]; // [alpha][gamma]
    let mut dm2:MatDeriv=[[[0.0;3];3];3]; // [alpha][beta][gamma]
    for (j,src) in p.iter().enumerate(){
        let d=query-src.pos_m;
        let r=d.len();
        if r>=2.0*h {continue;}
        let q=as_arr(d*(1.0/h));
        let volume=volumes[j];
        assert!(volume>0.0&&volume.is_finite());
        let w=volume*sph_kernel::kernel_and_dh(r,h,kernel).0;
        let gradw=if r>0.0 {
            as_arr(d*(volume*sph_kernel::radial_derivative(r,h,kernel)/r))
        } else {[0.0;3]};
        m0+=w;
        for alpha in 0..3 {
            m1[alpha]+=w*q[alpha];
            dm0[alpha]+=gradw[alpha];
            for gamma in 0..3 {
                dm1[alpha][gamma]+=gradw[gamma]*q[alpha]+
                    if alpha==gamma{w/h}else{0.0};
            }
            for beta in 0..3{
                m2[alpha][beta]+=w*q[alpha]*q[beta];
                for gamma in 0..3{
                    dm2[alpha][beta][gamma]+=
                        gradw[gamma]*q[alpha]*q[beta]+
                        (w/h)*((if alpha==gamma{q[beta]}else{0.0})+
                            (if beta==gamma{q[alpha]}else{0.0}));
                }
            }
        }
    }
    let inv=match invert_spd(&m2){
        Some(x)=>x,None=>return RkCoefficient::fallback(m0),
    };
    let mut b=mat_vec(&inv,m1);
    for x in &mut b {*x= -*x;}
    let denom=m0+dot(b,m1);
    if !denom.is_finite()||denom<=1e-12 {
        return RkCoefficient::fallback(m0);
    }
    let a=1.0/denom;
    let mut grad_a=[0.0;3];
    let mut grad_b=[[0.0;3];3];
    for gamma in 0..3{
        let mut rhs=[0.0;3];
        for alpha in 0..3{
            rhs[alpha]=dm1[alpha][gamma]+
                (0..3).map(|beta|dm2[alpha][beta][gamma]*b[beta])
                    .sum::<f64>();
        }
        let mut db=mat_vec(&inv,rhs);
        for x in &mut db {*x= -*x;}
        let ddenom=dm0[gamma]+dot(db,m1)+
            (0..3).map(|alpha|b[alpha]*dm1[alpha][gamma]).sum::<f64>();
        grad_a[gamma]=-a*a*ddenom;
        for alpha in 0..3{grad_b[alpha][gamma]=db[alpha];}
    }
    let result=RkCoefficient {
        a,b:from_arr(b),grad_a:from_arr(grad_a),
        grad_b:grad_b.map(from_arr),valid:true,raw_partition:m0,
    };
    if !result.a.is_finite()||!result.b.finite()||
       !result.grad_a.finite()||
       !result.grad_b.iter().all(|g|g.finite()){
        return RkCoefficient::fallback(m0);
    }
    result
}
/// At a query x, return corrected scalar W^R_ij and complete query gradient.
pub fn weight_gradient(
    x:V3,source:V3,h:f64,kernel:DensityKernel,c:RkCoefficient
)->(f64,V3) {
    let d=x-source;
    let r=d.len();
    if r>=2.0*h{return (0.0,V3::default());}
    let (w,_)=sph_kernel::kernel_and_dh(r,h,kernel);
    let gradw=if r>0.0{
        d*(sph_kernel::radial_derivative(r,h,kernel)/r)
    } else {V3::default()};
    if !c.valid{return (w,gradw);}
    let q=d*(1.0/h);
    let factor=1.0+c.b.dot(q);
    let grad_factor=c.b*(1.0/h)+
        c.grad_b[0]*q.x+c.grad_b[1]*q.y+c.grad_b[2]*q.z;
    let grad=gradw*(c.a*factor)+
        c.grad_a*(factor*w)+grad_factor*(c.a*w);
    (c.a*factor*w,grad)
}
#[derive(Clone,Debug)]
pub struct RkGeometry {
    pub coefficients:Vec<RkCoefficient>,
    pub volumes:Vec<f64>,
    pub support_h:Vec<f64>,
    pub kernel:DensityKernel,
}
#[derive(Clone,Copy,Debug)]
pub struct RkIdentityAudit{
    pub valid_fraction:f64,
    pub maximum_constant_error:f64,
    pub maximum_linear_error:f64,
    pub maximum_query_gradient_error:f64,
}
impl RkGeometry {
    pub fn from_state(p:&[Fluid3],state:&AdaptiveState,kernel:DensityKernel)->Self{
        assert!(state.valid_for(p.len()));
        let volumes=p.iter().zip(&state.densities_kg_m3)
            .map(|(x,rho)|x.mass_kg/rho).collect::<Vec<_>>();
        let coefficients=p.iter().zip(&state.smoothing_m)
            .map(|(x,h)|query_coefficient(x.pos_m,p,&volumes,*h,kernel))
            .collect::<Vec<_>>();
        Self{coefficients,volumes,support_h:state.smoothing_m.clone(),kernel}
    }
    /// Corrected W^R and gradient wrt i's query point, sources fixed.
    pub fn query(&self,p:&[Fluid3],i:usize,j:usize)->(f64,V3){
        weight_gradient(p[i].pos_m,p[j].pos_m,self.support_h[i],
            self.kernel,self.coefficients[i])
    }
    /// dWtilde_ij = 0.5*(grad_i W^R_ij-grad_j W^R_ji).
    /// dWtilde_ji = -dWtilde_ij by construction.
    pub fn antisymmetric_pair_gradient(
        &self,p:&[Fluid3],i:usize,j:usize,
    )->V3 {
        let (_,left)=self.query(p,i,j);
        let (_,right)=self.query(p,j,i);
        (left-right)*0.5
    }
    pub fn identities(&self,p:&[Fluid3])->RkIdentityAudit{
        assert_eq!(p.len(),self.coefficients.len());
        let mut valid=0;
        let mut max0=0.0_f64;
        let mut max1=0.0_f64;
        let mut maxg=0.0_f64;
        for (i,c) in self.coefficients.iter().enumerate(){
            if !c.valid {continue;}
            valid+=1;
            let mut constant=0.0;
            let mut first=V3::default();
            let mut gradient_sum=V3::default();
            for j in 0..p.len(){
                let (w,g)=self.query(p,i,j);
                let volume=self.volumes[j];
                constant+=volume*w;
                first+=(p[i].pos_m-p[j].pos_m)*
                    (volume*w/self.support_h[i]);
                gradient_sum+=g*volume*self.support_h[i];
            }
            max0=max0.max((constant-1.0).abs());
            max1=max1.max(first.len());
            maxg=maxg.max(gradient_sum.len());
        }
        RkIdentityAudit{
            valid_fraction:valid as f64/p.len() as f64,
            maximum_constant_error:max0,
            maximum_linear_error:max1,
            maximum_query_gradient_error:maxg,
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn seeded_grid()->(Vec<Fluid3>,Vec<f64>){
        let mut p=Vec::new();
        for x in -2..=2{
            for y in -2..=2{
                for z in -2..=2 {
                    p.push(Fluid3{
                        pos_m:V3::new(
                            x as f64*1.1e6+y as f64*5e4,
                            y as f64*1.1e6+z as f64*7e4,
                            z as f64*1.1e6+x as f64*9e4),
                        vel_m_s:V3::default(),mass_kg:1e27,
                    });
                }
            }
        }
        let volume=1.1e6_f64.powi(3);
        let len=p.len();
        (p,vec![volume;len])
    }
    #[test]
    fn query_reproduces_constants_and_linear_fields_at_boundary(){
        let (p,volume)=seeded_grid();
        let h=2.0e6;
        for i in [0,7,62,90,124]{
            let c=query_coefficient(p[i].pos_m,&p,&volume,h,DensityKernel::WendlandC4);
            assert!(c.valid,"singular RK at particle {i}");
            let mut zeroth=0.0;
            let mut first=V3::default();
            for (j,src) in p.iter().enumerate() {
                let (w,_)=weight_gradient(p[i].pos_m,src.pos_m,h,
                    DensityKernel::WendlandC4,c);
                zeroth+=volume[j]*w;
                first+=(p[i].pos_m-src.pos_m)*(volume[j]*w/h);
            }
            assert!((zeroth-1.0).abs()<1e-10,"zeroth {i} {zeroth}");
            assert!(first.len()<1e-10,"linear {i}: {:?}",first);
        }
    }
    #[test]
    fn rotation_covariance_for_full_query_gradients(){
        let (p,volumes)=seeded_grid();
        let rotated=p.iter().map(|x|Fluid3 {
            pos_m:V3::new(-x.pos_m.y,x.pos_m.x,x.pos_m.z),
            vel_m_s:V3::default(),mass_kg:x.mass_kg,
        }).collect::<Vec<_>>();
        let i=0;
        let j=7;
        let h=2.0e6;
        for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
            let a=query_coefficient(p[i].pos_m,&p,&volumes,h,kernel);
            let b=query_coefficient(rotated[i].pos_m,&rotated,&volumes,h,kernel);
            assert!(a.valid&&b.valid);
            let (w,g)=weight_gradient(p[i].pos_m,p[j].pos_m,h,kernel,a);
            let (wr,gr)=weight_gradient(rotated[i].pos_m,rotated[j].pos_m,
                h,kernel,b);
            let rotated_g=V3::new(-g.y,g.x,g.z);
            let scale=g.len().max(1e-40);
            assert!((w-wr).abs()/w.abs().max(1e-40)<1e-10);
            assert!((gr-rotated_g).len()/scale<1e-9);
        }
    }
    #[test]
    fn full_query_derivative_matches_central_difference(){
        let (p,volume)=seeded_grid();
        let i=7;
        let j=9;
        let h=2.3e6;
        let x=p[i].pos_m;
        let c=query_coefficient(x,&p,&volume,h,DensityKernel::WendlandC4);
        assert!(c.valid);
        let (_,grad)=weight_gradient(x,p[j].pos_m,h,
            DensityKernel::WendlandC4,c);
        let dx=h*1e-4;
        for axis in 0..3 {
            let delta=match axis{
                0=>V3::new(dx,0.0,0.0),
                1=>V3::new(0.0,dx,0.0),
                _=>V3::new(0.0,0.0,dx),
            };
            let plus=query_coefficient(x+delta,&p,&volume,h,
                DensityKernel::WendlandC4);
            let minus=query_coefficient(x-delta,&p,&volume,h,
                DensityKernel::WendlandC4);
            let fp=weight_gradient(x+delta,p[j].pos_m,h,
                DensityKernel::WendlandC4,plus).0;
            let fm=weight_gradient(x-delta,p[j].pos_m,h,
                DensityKernel::WendlandC4,minus).0;
            let numeric=(fp-fm)/(2.0*dx);
            let g=[grad.x,grad.y,grad.z][axis];
            let scale=fp.abs().max(fm.abs())/h;
            assert!((numeric-g).abs()<5e-4*scale.max(1e-40),
                "axis {axis} analytic {g:e} numeric {numeric:e}");
        }
    }
}
