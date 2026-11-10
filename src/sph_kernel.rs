//! BG-PHY-16 density and Poisson-matched gravitational smoothing kernels.
//!
//! Both choices have compact support q=r/h<2. The M4 path forwards to
//! existing stable code. Wendland C4 is the normalized three-dimensional
//! positive-Fourier-transform kernel used by Dehnen & Aly (2012).
//!
//! W_C4 = 495/(256*pi*h^3)*(1-q/2)^6*(1+3q+35q^2/12).
//!
//! The gravitational kernel is derived by solving the spherical Poisson
//! equation for exactly this W, with phi=-1/r outside r>=2h. Its grad-h
//! derivative closes the same variational force terms as M4.
use std::f64::consts::PI;
use crate::sph3;
use crate::sph_gravity_kernel;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum DensityKernel {
    CubicM4,
    WendlandC4,
}
impl DensityKernel {
    pub fn label(self)->&'static str {
        match self {Self::CubicM4=>"cubic_m4", Self::WendlandC4=>"wendland_c4"}
    }
}
fn c4_shape(q:f64)->(f64,f64) {
    if q>=2.0 {return (0.0,0.0);}
    let s=1.0-0.5*q;
    let shape=s.powi(6)*(1.0+3.0*q+(35.0/12.0)*q*q);
    let dshape=-(7.0/3.0)*q*(2.0+5.0*q)*s.powi(5);
    (shape,dshape)
}
/// Return W and dW/dh at fixed physical separation.
pub fn kernel_and_dh(r:f64,h:f64,k:DensityKernel)->(f64,f64) {
    assert!(r>=0.0&&r.is_finite()&&h>0.0&&h.is_finite());
    let q=r/h;
    if q>=2.0 {return (0.0,0.0);}
    let (shape,derivative,normalization)=match k {
        DensityKernel::CubicM4=>{
            let (f,fp)=if q<1.0 {
                (1.0-1.5*q*q+0.75*q.powi(3),-3.0*q+2.25*q*q)
            } else {
                (0.25*(2.0-q).powi(3),-0.75*(2.0-q).powi(2))
            };
            (f,fp,1.0/PI)
        }
        DensityKernel::WendlandC4=>{
            let (f,fp)=c4_shape(q);
            (f,fp,495.0/(256.0*PI))
        }
    };
    let w=normalization*shape/h.powi(3);
    let dwdh=-normalization*(3.0*shape+q*derivative)/h.powi(4);
    (w,dwdh)
}
/// Scalar radial derivative dW/dr (always nonpositive for r>0).
pub fn radial_derivative(r:f64,h:f64,k:DensityKernel)->f64 {
    assert!(r>=0.0&&r.is_finite()&&h>0.0&&h.is_finite());
    match k {
        DensityKernel::CubicM4=>sph3::spline_radial_derivative_3d(r,h),
        DensityKernel::WendlandC4=>{
            let (_,df)=c4_shape(r/h);
            (495.0/(256.0*PI))*df/h.powi(4)
        }
    }
}
/// Exact enclosed-mass fraction M(<r) / M for a single C4 density kernel.
fn c4_enclosed(q:f64)->f64 {
    if q>=2.0 {return 1.0;}
    if q<=0.0 {return 0.0;}
    if q<=1.0 {
        let q2=q*q;
        q.powi(3)/16384.0 *
            (42240.0+q2*(-59136.0+q2*(79200.0+
                q*(-73920.0+q*(30800.0+q*(-6336.0+525.0*q))))))
    } else {
        let q2=q*q;
        1.0-(2.0-q).powi(7)*(525.0*q2*q2+
            1014.0*q2*q+896.0*q2+448.0*q+128.0)/16384.0
    }
}
/// phi(r,h)/Gm; negative with continuous Newtonian -1/r exterior.
pub fn potential(r:f64,h:f64,k:DensityKernel)->f64 {
    assert!(r>=0.0&&r.is_finite()&&h>0.0&&h.is_finite());
    match k {
        DensityKernel::CubicM4=>sph_gravity_kernel::potential(r,h),
        DensityKernel::WendlandC4=>{
            let q=r/h;
            if q>=2.0 {return -1.0/r;}
            let p=-55.0/32.0 + q*q*(
                165.0/128.0+q*q*(-231.0/256.0+q*q*(
                825.0/1024.0+q*(-165.0/256.0+
                q*(1925.0/8192.0+q*(-11.0/256.0+q*(105.0/32768.0)))))));
            p/h
        }
    }
}
/// Positive dphi/dr = M(<r) / r^2, with the exact exterior Newtonian limit.
pub fn gravity_radial_derivative(r:f64,h:f64,k:DensityKernel)->f64 {
    assert!(r>=0.0&&r.is_finite()&&h>0.0&&h.is_finite());
    match k {
        DensityKernel::CubicM4=>sph_gravity_kernel::radial_derivative(r,h),
        DensityKernel::WendlandC4=>{
            if r==0.0 {0.0} else {c4_enclosed(r/h)/(r*r)}
        }
    }
}
/// Positive dphi/dh at fixed r, vanishing identically for r>=2h.
pub fn gravity_smoothing_derivative(r:f64,h:f64,k:DensityKernel)->f64 {
    assert!(r>=0.0&&r.is_finite()&&h>0.0&&h.is_finite());
    match k {
        DensityKernel::CubicM4=>sph_gravity_kernel::smoothing_derivative(r,h),
        DensityKernel::WendlandC4=>{
            let q=r/h;
            if q>=2.0 {return 0.0;}
            55.0*(2.0-q).powi(7)*(3.0*q+2.0)*
                (7.0*q*q+8.0*q+4.0)/(32768.0*h*h)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cubic_path_preserves_old_density_force_and_potential() {
        let h=2.1e6;
        for q in [0.0,0.4,0.99,1.0,1.9,2.0,2.4] {
            let r=q*h;
            let (w,_)=kernel_and_dh(r,h,DensityKernel::CubicM4);
            assert!((w-sph3::cubic_spline_3d(r,h)).abs()<1e-24);
            assert_eq!(radial_derivative(r,h,DensityKernel::CubicM4),
                sph3::spline_radial_derivative_3d(r,h));
            assert_eq!(potential(r,h,DensityKernel::CubicM4),
                sph_gravity_kernel::potential(r,h));
        }
    }
    #[test]
    fn c4_spherical_normalization_and_poisson_mass() {
        let n=1600;
        let h=2.0;
        let delta=2.0*h/n as f64;
        let mut sum=0.0;
        for i in 0..=n {
            let r=i as f64*delta;
            let w=kernel_and_dh(r,h,DensityKernel::WendlandC4).0;
            let weight=if i==0||i==n{1.0}else if i%2==0{2.0}else{4.0};
            sum+=weight*4.0*PI*r*r*w;
        }
        let integral=sum*delta/3.0;
        assert!((integral-1.0).abs()<2e-10,"C4 normalization {integral}");
        for q in [0.05,0.2,0.7,1.0,1.5,1.9,2.0,2.5] {
            let r=q*h;
            let phi=potential(r,h,DensityKernel::WendlandC4);
            let radial=gravity_radial_derivative(r,h,DensityKernel::WendlandC4);
            let smooth=gravity_smoothing_derivative(r,h,DensityKernel::WendlandC4);
            assert!(phi<0.0&&radial>0.0&&smooth>=0.0);
            let dr=1e-5*h;
            let dphidr=(potential(r+dr,h,DensityKernel::WendlandC4)-
                potential(r-dr,h,DensityKernel::WendlandC4))/(2.0*dr);
            assert!((dphidr-radial).abs()<2e-6*radial.max(1e-12),
                "radial mismatch q={q}: {dphidr} vs {radial}");
            let dh=1e-5*h;
            let dphidh=(potential(r,h+dh,DensityKernel::WendlandC4)-
                potential(r,h-dh,DensityKernel::WendlandC4))/(2.0*dh);
            assert!((dphidh-smooth).abs()<5e-6*smooth.max(1e-12),
                "grad-h mismatch q={q}: {dphidh} vs {smooth}");
        }
    }
    #[test]
    fn c4_kernel_h_derivative_matches_finite_difference() {
        let h=1.8e6;
        for q in [0.0,0.2,0.7,1.15,1.8] {
            let r=q*h;let dh=2.0;
            let (w,analytic)=kernel_and_dh(r,h,DensityKernel::WendlandC4);
            let numeric=(kernel_and_dh(r,h+dh,DensityKernel::WendlandC4).0-
                kernel_and_dh(r,h-dh,DensityKernel::WendlandC4).0)/(2.0*dh);
            assert!(w>=0.0);
            assert!((analytic-numeric).abs()<1e-5*analytic.abs().max(1e-35),
                "dW/dh mismatch q={q}: {analytic} vs {numeric}");
        }
    }
}
