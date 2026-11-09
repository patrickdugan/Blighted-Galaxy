//! Compact-support cubic-spline gravitational softening matched to the
//! *density* M4 SPH kernel used by sph3.rs.
 //!
 //! Price & Monaghan (2007), MNRAS 374, Appendix A, equations A1–A3.
 //! W_3D=1/(pi h^3) [1 - 1.5q² + 0.75q³] for 0<=q<1,
 //! W_3D=1/(4 pi h^3) (2-q)^3 for 1<=q<2, q=r/h.
 //!
 //! phi(r,h) is the unit negative gravitational potential (-1/r at
 //! r>=2h); dphi/dr is positive, and dphi/dh controls adaptive softening
 //! with zeta/Omega. The kernel is exactly Newtonian at >=2h.
pub fn potential(r_m:f64,h_m:f64)->f64 {
    assert!(r_m>=0.0&&r_m.is_finite()&&h_m>0.0&&h_m.is_finite());
    let q=r_m/h_m;
    if q<1.0 {
        let q2=q*q;let q4=q2*q2;
        ((2.0/3.0)*q2-0.3*q4+0.1*q4*q-1.4)/h_m
    } else if q<2.0 {
        let q2=q*q;
        ((4.0/3.0)*q2-q2*q+0.3*q2*q2-
            (1.0/30.0)*q2*q2*q-1.6+1.0/(15.0*q))/h_m
    } else {
        -1.0/r_m
    }
}
/// Positive dphi/dr, SI 1/m². The vector gravitational acceleration
/// on i is +G*m_j*(r_j-r_i)/r * dphi/dr.
pub fn radial_derivative(r_m:f64,h_m:f64)->f64 {
    assert!(r_m>=0.0&&r_m.is_finite()&&h_m>0.0&&h_m.is_finite());
    let q=r_m/h_m;
    if q<1.0 {
        ((4.0/3.0)*q-1.2*q.powi(3)+0.5*q.powi(4))/h_m.powi(2)
    } else if q<2.0 {
        ((8.0/3.0)*q-3.0*q*q+1.2*q.powi(3)-
            q.powi(4)/6.0-1.0/(15.0*q*q))/h_m.powi(2)
    } else {
        1.0/(r_m*r_m)
    }
}
/// Positive partial derivative dphi/dh (SI 1/m²) for fixed r.
/// It vanishes outside the support of the mass distribution.
pub fn smoothing_derivative(r_m:f64,h_m:f64)->f64 {
    assert!(r_m>=0.0&&r_m.is_finite()&&h_m>0.0&&h_m.is_finite());
    let q=r_m/h_m;
    if q<1.0 {
        (-2.0*q*q+1.5*q.powi(4)-0.6*q.powi(5)+1.4)/h_m.powi(2)
    } else if q<2.0 {
        (-4.0*q*q+4.0*q.powi(3)-1.5*q.powi(4)+
            0.2*q.powi(5)+1.6)/h_m.powi(2)
    } else {
        0.0
    }
}
/// Newtonian enclosed-mass fraction of the unit-normalized M4 kernel.
pub fn enclosed_mass_fraction(r_m:f64,h_m:f64)->f64 {
    (r_m*r_m*radial_derivative(r_m,h_m)).clamp(0.0,1.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sph3;
    fn near(x:f64,y:f64,tol:f64) {
        assert!((x-y).abs()<=tol*y.abs().max(1e-30),
            "observed {x:.15e}, expected {y:.15e}");
    }
    #[test]
    fn analytic_softening_is_regular_and_matches_newton_at_two_h() {
        let h=2.0e6;
        near(potential(0.0,h),-1.4/h,1e-13);
        near(radial_derivative(0.0,h),0.0,1e-14);
        near(potential(2.0*h,h),-1.0/(2.0*h),1e-14);
        near(radial_derivative(2.0*h,h),1.0/(4.0*h*h),1e-14);
        assert_eq!(smoothing_derivative(2.0*h,h),0.0);
        for q in [0.0,0.1,0.7,1.0,1.5,1.99,2.0,2.5] {
            let r=q*h;
            assert!(potential(r,h)<0.0);
            assert!(radial_derivative(r,h)>=0.0);
            assert!(smoothing_derivative(r,h)>=-1e-25);
        }
    }
    #[test]
    fn piecewise_potential_and_force_are_continuous_at_breaks() {
        let h=3.0e6;
        for q in [1.0,2.0] {
            let r=q*h;let dr=1e-8*h;
            let reference=potential(r,h);
            near(potential(r-dr,h),reference,1e-7);
            near(potential(r+dr,h),reference,1e-7);
            near(radial_derivative(r-dr,h),radial_derivative(r,h),1e-7);
            near(radial_derivative(r+dr,h),radial_derivative(r,h),1e-7);
        }
    }
    #[test]
    fn potential_derivatives_match_finite_differences() {
        let h=2.4e6;
        for q in [0.03,0.4,0.92,1.2,1.82,2.3] {
            let r=q*h;let dr=h*1e-5;let dh=h*1e-5;
            let fd_r=(potential(r+dr,h)-potential(r-dr,h))/(2.0*dr);
            let fd_h=(potential(r,h+dh)-potential(r,h-dh))/(2.0*dh);
            near(radial_derivative(r,h),fd_r,2e-5);
            if q!=2.3 {
                near(smoothing_derivative(r,h),fd_h,5e-5);
            } else {
                near(fd_h,0.0,1e-12);
            }
        }
    }
    #[test]
    fn force_obeys_poisson_enclosed_mass_of_identical_cubic_kernel() {
        let h=1.8e6;
        let intervals=12_000;
        let dq=2.0/intervals as f64;
        let mut cumulative=0.0;
        let probes=[0.2,0.6,0.95,1.2,1.8,2.0];
        for j in 0..intervals {
            let q=(j as f64+0.5)*dq;
            let r=q*h;
            cumulative+=4.0*std::f64::consts::PI*
                sph3::cubic_spline_3d(r,h)*r*r*h*dq;
            for &target in &probes {
                if ((j+1) as f64*dq-target).abs()<dq*0.51 {
                    let enclosed=enclosed_mass_fraction(target*h,h);
                    near(enclosed,cumulative,8e-5);
                }
            }
        }
        near(cumulative,1.0,1e-7);
        near(enclosed_mass_fraction(3.0*h,h),1.0,1e-14);
    }
}
