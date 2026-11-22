//! BG-PHY-26: joint three-moment free-surface depth and normal volume-bias
//! identification for truncated isotropic SPH kernels.
//!
//! Given outward normal `n_out = -Q_i / |Q_i|` from the raw first moment
//! `Q_i = sum_j (m_j / rho_j) W_ij (x_j - x_i) / h_i`, let
//! `z_ij = n_out . (x_j - x_i) / h_i` and measure
//!
//!   M0 = sum_j (m_j / rho_j) W_ij
//!   M1 = sum_j (m_j / rho_j) W_ij z_ij
//!   M2 = sum_j (m_j / rho_j) W_ij z_ij^2.
//!
//! Under the local planar linear-volume-bias model
//! `V_est(z) = c * (1 + beta * z) * V_true(z)` truncated at `z = a` (`0 <= a <= 2`),
//!
//!   M0 = c * [I0(a) + beta * I1(a)]
//!   M1 = c * [I1(a) + beta * I2(a)]
//!   M2 = c * [I2(a) + beta * I3(a)],
//!
//! where `I_k(a) = integral_{|q|<2, q.n_out<=a} (q.n_out)^k W(|q|, 1) d^3q`.
//!
//! Eliminating `c` via `r1 = M1/M0`, `r2 = M2/M0` yields
//!   beta(a) = [I1(a) - r1 * I0(a)] / [r1 * I1(a) - I2(a)],
//!   F(a)    = [I2(a) + beta(a) * I3(a)] / [I0(a) + beta(a) * I1(a)] - r2.
//!
//! Algebraically clearing the removable factor `r1 * I1(a) - I2(a)` gives the
//! smooth rational residual
//!   r2_hat(a) = [(I2^2 - I1*I3) + r1*(I0*I3 - I1*I2)] / (I0*I2 - I1^2),
//! whose denominator `I0*I2 - I1^2 > 0` strictly on `[0, 2]` by Cauchy-Schwarz.
//! Every candidate zero is bracketed across `[0, 2]`, polished to machine
//! precision, verified directly on the original rational `beta(a)` and `F(a)`
//! equations, and checked for `c > 0`, `1 + beta*z > 0` across the occupied
//! half-space `[-2, a]` (and realized neighbor samples), uniqueness, and
//! Jacobian conditioning `|dF/da| >= min_df_da`. Ambiguous or ill-conditioned
//! solutions are marked unidentified and never impose a vacuum `H = 0` ghost.

use crate::sph3::{Fluid3, V3};
use crate::sph_adaptive::AdaptiveState;
use crate::sph_kernel::{self, DensityKernel};
use crate::sph_neighbors::CellNeighbors;
use crate::white_dwarf::DegenerateElectrons;
use std::f64::consts::PI;

/// 8-point Gauss-Legendre nodes and weights on `[-1, 1]`, exact for polynomials
/// up to degree 15 (maximum integrand degree for C4 `I_3` is `8 + 5 = 13`).
const GL8_NODES: [f64; 8] = [
    -0.960_289_856_497_536_3,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_329,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const GL8_WEIGHTS: [f64; 8] = [
    0.101_228_536_290_376_26,
    0.222_381_034_453_374_47,
    0.313_706_645_877_887_27,
    0.362_683_783_378_362,
    0.362_683_783_378_362,
    0.313_706_645_877_887_27,
    0.222_381_034_453_374_47,
    0.101_228_536_290_376_26,
];

const BRACKET_INTERVALS: usize = 256;
pub const DEFAULT_MIN_DF_DA: f64 = 1e-3;
const POSITIVITY_FLOOR: f64 = 1e-4;

/// Stage B gate: inferred last-shell depth must be identified on >= 90% of
/// outermost-shell particles and achieve mean absolute depth error <= 0.15h
/// before a vacuum H=0 ghost can be promoted as physically trustworthy.
pub const GEOMETRY_GATE_MIN_IDENTIFIED_FRACTION: f64 = 0.90;
pub const GEOMETRY_GATE_MAX_DEPTH_ABS_ERROR_H: f64 = 0.15;

#[derive(Clone, Copy, Debug)]
pub struct PlanarMoments {
    pub a: f64,
    pub i: [f64; 4],
    pub di_da: [f64; 4],
}

#[derive(Clone, Debug)]
pub struct MomentTable {
    pub kernel: DensityKernel,
    pub norm: f64,
    pub grid: Vec<PlanarMoments>,
}

fn gl8_segment<F: Fn(f64) -> f64>(lo: f64, hi: f64, f: F) -> f64 {
    if hi <= lo {
        return 0.0;
    }
    let mid = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let mut acc = 0.0;
    for (&x, &w) in GL8_NODES.iter().zip(&GL8_WEIGHTS) {
        acc += w * f(mid + half * x);
    }
    acc * half
}

fn sorted_knots_for_a(a: f64, kernel: DensityKernel) -> ([f64; 4], usize) {
    let a_clamped = a.clamp(0.0, 2.0);
    match kernel {
        DensityKernel::CubicM4 => {
            let mut k = [0.0, 1.0, 2.0, a_clamped];
            k.sort_by(|x, y| x.total_cmp(y));
            (k, 4)
        }
        DensityKernel::WendlandC4 => {
            let mut k = [0.0, a_clamped, 2.0, 2.0];
            k[..3].sort_by(|x, y| x.total_cmp(y));
            (k, 3)
        }
    }
}

/// Unnormalized spherical integral `I_0(2) = int_{|q|<2} W(q, 1) d^3q`.
pub fn kernel_volume_norm(kernel: DensityKernel) -> f64 {
    let (knots, len) = sorted_knots_for_a(2.0, kernel);
    let mut total = 0.0;
    for seg in 0..len - 1 {
        let (lo, hi) = (knots[seg], knots[seg + 1]);
        if hi - lo <= 1e-15 {
            continue;
        }
        total += gl8_segment(lo, hi, |q| {
            4.0 * PI * q * q * sph_kernel::kernel_and_dh(q, 1.0, kernel).0
        });
    }
    total
}

/// Exact 1D piecewise Gauss-Legendre quadrature for `I_0..I_3(a)` and exact
/// Leibniz derivatives `dI_k/da = 2*pi*a^k * int_a^2 q W(q, 1) dq / I_0(2)`.
pub fn planar_moments_and_deriv_with_norm(
    a: f64,
    kernel: DensityKernel,
    norm: f64,
) -> PlanarMoments {
    assert!((0.0..=2.0).contains(&a) && a.is_finite());
    assert!(norm > 0.0 && norm.is_finite());
    let (knots, len) = sorted_knots_for_a(a, kernel);
    let mut moments = [0.0_f64; 4];
    let mut slice_mass = 0.0_f64;
    for seg in 0..len - 1 {
        let (lo, hi) = (knots[seg], knots[seg + 1]);
        if hi - lo <= 1e-15 {
            continue;
        }
        let mid = 0.5 * (lo + hi);
        let inside_sphere = mid <= a;
        for (k, slot) in moments.iter_mut().enumerate() {
            let ku = k as i32;
            *slot += gl8_segment(lo, hi, |q| {
                let w = sph_kernel::kernel_and_dh(q, 1.0, kernel).0;
                if inside_sphere {
                    if k % 2 == 1 {
                        0.0
                    } else {
                        4.0 * PI * q.powi(ku + 2) * w / (k as f64 + 1.0)
                    }
                } else {
                    let sign = if (k + 1) % 2 == 0 { 1.0 } else { -1.0 };
                    let poly = a.powi(ku + 1) * q - sign * q.powi(ku + 2);
                    2.0 * PI * w * poly / (k as f64 + 1.0)
                }
            });
        }
        if !inside_sphere {
            slice_mass += gl8_segment(lo, hi, |q| {
                2.0 * PI * q * sph_kernel::kernel_and_dh(q, 1.0, kernel).0
            });
        }
    }
    for m in &mut moments {
        *m /= norm;
    }
    slice_mass /= norm;
    let di_da = if a >= 2.0 {
        [0.0; 4]
    } else {
        [
            slice_mass,
            a * slice_mass,
            a * a * slice_mass,
            a * a * a * slice_mass,
        ]
    };
    PlanarMoments {
        a,
        i: moments,
        di_da,
    }
}

pub fn planar_moments_and_deriv(a: f64, kernel: DensityKernel) -> PlanarMoments {
    let norm = kernel_volume_norm(kernel);
    planar_moments_and_deriv_with_norm(a, kernel, norm)
}

pub fn planar_moments(a: f64, kernel: DensityKernel) -> [f64; 4] {
    planar_moments_and_deriv(a, kernel).i
}

impl MomentTable {
    pub fn new(kernel: DensityKernel) -> Self {
        let norm = kernel_volume_norm(kernel);
        assert!((norm - 1.0).abs() < 1e-13, "kernel normalization drift: {norm}");
        let grid = (0..=BRACKET_INTERVALS)
            .map(|idx| {
                let a = 2.0 * idx as f64 / BRACKET_INTERVALS as f64;
                planar_moments_and_deriv_with_norm(a, kernel, norm)
            })
            .collect();
        Self { kernel, norm, grid }
    }

    pub fn eval(&self, a: f64) -> PlanarMoments {
        planar_moments_and_deriv_with_norm(a, self.kernel, self.norm)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThreeMomentStatus {
    Identified,
    IllConditionedInterior,
    NoPhysicalRoot,
    AmbiguousMultipleRoots,
    InvalidMoments,
}

impl ThreeMomentStatus {
    pub fn label(self) -> &'static str {
        match self {
            Self::Identified => "identified",
            Self::IllConditionedInterior => "ill_conditioned_interior",
            Self::NoPhysicalRoot => "no_physical_root",
            Self::AmbiguousMultipleRoots => "ambiguous_multiple_roots",
            Self::InvalidMoments => "invalid_moments",
        }
    }
    pub fn is_identified(self) -> bool {
        self == Self::Identified
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CandidateRoot {
    pub a: f64,
    pub c: f64,
    pub beta: f64,
    pub residual_f: f64,
    pub df_da: f64,
    pub depth_sensitivity: f64,
    pub min_volume_factor_halfspace: f64,
    pub min_volume_factor_samples: f64,
}

#[derive(Clone, Debug)]
pub struct ThreeMomentFit {
    pub status: ThreeMomentStatus,
    pub root: Option<CandidateRoot>,
    pub physical_root_count: usize,
    pub r1: f64,
    pub r2: f64,
}

/// Smooth removable-pole-free residual `r2_hat(a) - r2` and its exact derivative.
fn smooth_residual_and_deriv(pm: &PlanarMoments, r1: f64, r2: f64) -> (f64, f64) {
    let [i0, i1, i2, i3] = pm.i;
    let [d0, d1, d2, d3] = pm.di_da;
    let num = (i2 * i2 - i1 * i3) + r1 * (i0 * i3 - i1 * i2);
    let den = (i0 * i2 - i1 * i1).max(1e-30);
    let dnum = (2.0 * i2 * d2 - d1 * i3 - i1 * d3)
        + r1 * (d0 * i3 + i0 * d3 - d1 * i2 - i1 * d2);
    let dden = d0 * i2 + i0 * d2 - 2.0 * i1 * d1;
    let val = num / den - r2;
    let deriv = (dnum * den - num * dden) / (den * den);
    (val, deriv)
}

/// Evaluate the exact Issue #31 rational equations at trial `a` and check
/// positivity and Jacobian conditioning.
fn verify_root(
    a: f64,
    m0: f64,
    r1: f64,
    r2: f64,
    sample_z_bounds: Option<(f64, f64)>,
    table: &MomentTable,
) -> Option<CandidateRoot> {
    if !(0.0..=2.0).contains(&a) || !a.is_finite() {
        return None;
    }
    let pm = table.eval(a);
    let [i0, i1, i2, i3] = pm.i;
    let [d0, d1, d2, d3] = pm.di_da;
    let denom_beta = r1 * i1 - i2;
    if !denom_beta.is_finite() || denom_beta.abs() <= 1e-12 {
        return None;
    }
    let num_beta = i1 - r1 * i0;
    let beta = num_beta / denom_beta;
    if !beta.is_finite() {
        return None;
    }
    let denom_f = i0 + beta * i1;
    if !denom_f.is_finite() || denom_f <= 1e-12 {
        return None;
    }
    let num_f = i2 + beta * i3;
    let residual_f = num_f / denom_f - r2;
    if !residual_f.is_finite() || residual_f.abs() > 1e-9 {
        return None;
    }
    let c = m0 / denom_f;
    if !c.is_finite() || c <= 0.0 {
        return None;
    }
    // Occupied half-space inside compact support |q| <= 2 has z in [-2, a].
    let min_halfspace = (1.0 - 2.0 * beta).min(1.0 + a * beta);
    if min_halfspace <= POSITIVITY_FLOOR {
        return None;
    }
    let min_samples = if let Some((z_min, z_max)) = sample_z_bounds {
        (1.0 + beta * z_min).min(1.0 + beta * z_max)
    } else {
        min_halfspace
    };
    if min_samples <= POSITIVITY_FLOOR {
        return None;
    }
    let dbeta = ((d1 - r1 * d0) * denom_beta - num_beta * (r1 * d1 - d2))
        / (denom_beta * denom_beta);
    let dnum_f = d2 + dbeta * i3 + beta * d3;
    let dden_f = d0 + dbeta * i1 + beta * d1;
    let df_da = (dnum_f * denom_f - num_f * dden_f) / (denom_f * denom_f);
    if !df_da.is_finite() {
        return None;
    }
    let depth_sensitivity = if df_da.abs() > 1e-30 {
        1.0 / df_da.abs()
    } else {
        f64::INFINITY
    };
    Some(CandidateRoot {
        a,
        c,
        beta,
        residual_f,
        df_da,
        depth_sensitivity,
        min_volume_factor_halfspace: min_halfspace,
        min_volume_factor_samples: min_samples,
    })
}

fn polish_bracket(
    mut lo: f64,
    mut hi: f64,
    mut flo: f64,
    mut fhi: f64,
    r1: f64,
    r2: f64,
    table: &MomentTable,
) -> f64 {
    if flo == 0.0 {
        return lo;
    }
    if fhi == 0.0 {
        return hi;
    }
    let mut a = 0.5 * (lo + hi);
    for _ in 0..64 {
        if (hi - lo).abs() <= 1e-14 {
            break;
        }
        let pm = table.eval(a);
        let (fa, dfa) = smooth_residual_and_deriv(&pm, r1, r2);
        if fa.abs() <= 1e-15 {
            return a;
        }
        if (flo < 0.0 && fa > 0.0) || (flo > 0.0 && fa < 0.0) {
            hi = a;
            fhi = fa;
        } else {
            lo = a;
            flo = fa;
        }
        let newton = if dfa.abs() > 1e-14 {
            a - fa / dfa
        } else {
            f64::NAN
        };
        a = if newton.is_finite() && newton > lo + 1e-14 && newton < hi - 1e-14 {
            newton
        } else {
            let _ = fhi;
            0.5 * (lo + hi)
        };
    }
    0.5 * (lo + hi)
}

/// Jointly infer `(a, c, beta)` from `(M0, M1, M2)` with safeguarded brackets,
/// positivity checks over `[-2, a]` and `sample_z_bounds`, uniqueness check,
/// and Jacobian conditioning gate `|dF/da| >= min_df_da`.
pub fn solve_three_moments_with_table(
    m0: f64,
    m1: f64,
    m2: f64,
    sample_z_bounds: Option<(f64, f64)>,
    table: &MomentTable,
    min_df_da: f64,
) -> ThreeMomentFit {
    if !m0.is_finite() || !m1.is_finite() || !m2.is_finite() || m0 <= 0.0 || m2 <= 0.0 {
        return ThreeMomentFit {
            status: ThreeMomentStatus::InvalidMoments,
            root: None,
            physical_root_count: 0,
            r1: f64::NAN,
            r2: f64::NAN,
        };
    }
    if m1 * m1 > m0 * m2 * (1.0 + 1e-11) {
        return ThreeMomentFit {
            status: ThreeMomentStatus::InvalidMoments,
            root: None,
            physical_root_count: 0,
            r1: m1 / m0,
            r2: m2 / m0,
        };
    }
    let r1 = m1 / m0;
    let r2 = m2 / m0;
    let vals: Vec<f64> = table
        .grid
        .iter()
        .map(|pm| smooth_residual_and_deriv(pm, r1, r2).0)
        .collect();

    let mut raw_roots: Vec<f64> = Vec::new();
    if vals[0].abs() <= 1e-13 {
        raw_roots.push(0.0);
    }
    for idx in 0..BRACKET_INTERVALS {
        let (f0, f1) = (vals[idx], vals[idx + 1]);
        // Skip flat machine-epsilon collars near a -> 2 where (2-a)^8 < 1e-15.
        if f0.abs().max(f1.abs()) <= 1e-13 {
            continue;
        }
        let (a0, a1) = (table.grid[idx].a, table.grid[idx + 1].a);
        if f0 == 0.0 {
            raw_roots.push(a0);
        } else if f1 == 0.0 {
            raw_roots.push(a1);
        } else if (f0 < 0.0 && f1 > 0.0) || (f0 > 0.0 && f1 < 0.0) {
            raw_roots.push(polish_bracket(a0, a1, f0, f1, r1, r2, table));
        }
    }
    if vals[BRACKET_INTERVALS].abs() <= 1e-13 {
        // Collapse any root inside the flat machine-epsilon collar [1.98, 2.0] to 2.0.
        raw_roots.retain(|&r| r < 1.98);
        raw_roots.push(2.0);
    }

    let mut dedup: Vec<f64> = Vec::new();
    for r in raw_roots {
        if dedup.iter().all(|&prev| (prev - r).abs() > 1e-5) {
            dedup.push(r);
        }
    }

    let physical: Vec<CandidateRoot> = dedup
        .into_iter()
        .filter_map(|a| verify_root(a, m0, r1, r2, sample_z_bounds, table))
        .collect();

    match physical.len() {
        0 => ThreeMomentFit {
            status: ThreeMomentStatus::NoPhysicalRoot,
            root: None,
            physical_root_count: 0,
            r1,
            r2,
        },
        1 => {
            let root = physical[0];
            let status = if root.df_da.abs() >= min_df_da {
                ThreeMomentStatus::Identified
            } else {
                ThreeMomentStatus::IllConditionedInterior
            };
            ThreeMomentFit {
                status,
                root: Some(root),
                physical_root_count: 1,
                r1,
                r2,
            }
        }
        count => ThreeMomentFit {
            status: ThreeMomentStatus::AmbiguousMultipleRoots,
            root: None,
            physical_root_count: count,
            r1,
            r2,
        },
    }
}

pub fn solve_three_moments(
    m0: f64,
    m1: f64,
    m2: f64,
    kernel: DensityKernel,
) -> ThreeMomentFit {
    let table = MomentTable::new(kernel);
    solve_three_moments_with_table(m0, m1, m2, None, &table, DEFAULT_MIN_DF_DA)
}

/// Stage A.3 diagnostic: evaluate spherical-cap moments of radius `r_sphere_h`
/// (in units of `h`) when the query particle sits `a_true_h` inside the sphere,
/// isolating pure surface curvature from density gradients and particle disorder.
pub fn spherical_cap_moments(
    a_true_h: f64,
    r_sphere_h: f64,
    c_scale: f64,
    beta_slope: f64,
    kernel: DensityKernel,
) -> [f64; 3] {
    assert!(r_sphere_h > a_true_h && a_true_h >= 0.0);
    let r_particle = r_sphere_h - a_true_h;
    let norm = kernel_volume_norm(kernel);
    let nq = 128usize;
    let dq = 2.0 / nq as f64;
    let mut out = [0.0_f64; 3];
    for iq in 0..nq {
        let q = (iq as f64 + 0.5) * dq;
        let w = sph_kernel::kernel_and_dh(q, 1.0, kernel).0 / norm;
        if w <= 0.0 {
            continue;
        }
        // |r_particle * e_z + q * u|^2 <= r_sphere_h^2
        // => r_p^2 + 2 r_p q mu + q^2 <= R^2 => mu <= (R^2 - r_p^2 - q^2)/(2 r_p q).
        let mu_cap = if r_particle > 1e-12 {
            ((r_sphere_h * r_sphere_h - r_particle * r_particle - q * q)
                / (2.0 * r_particle * q))
                .clamp(-1.0, 1.0)
        } else if q <= r_sphere_h {
            1.0
        } else {
            -1.0
        };
        if mu_cap <= -1.0 {
            continue;
        }
        let a0 = mu_cap + 1.0;
        let a1 = 0.5 * (mu_cap * mu_cap - 1.0);
        let a2 = (mu_cap.powi(3) + 1.0) / 3.0;
        let a3 = 0.25 * (mu_cap.powi(4) - 1.0);
        let shell = 2.0 * PI * q * q * w * dq;
        let i0 = shell * a0;
        let i1 = shell * q * a1;
        let i2 = shell * q * q * a2;
        let i3 = shell * q.powi(3) * a3;
        out[0] += c_scale * (i0 + beta_slope * i1);
        out[1] += c_scale * (i1 + beta_slope * i2);
        out[2] += c_scale * (i2 + beta_slope * i3);
    }
    out
}

/// Stage A.3 diagnostic: planar half-space moments under a nonlinear quadratic
/// volume bias `V_est / V_true = c * (1 + beta * z + gamma * z^2)` (modeling steep
/// convex degenerate-surface density falloff where SPH smoothing overestimates
/// surface density and creates curvature `gamma != 0` in `m / rho_SPH`).
pub fn quadratic_volume_bias_moments(
    a_true_h: f64,
    c_scale: f64,
    beta_slope: f64,
    gamma_curv: f64,
    kernel: DensityKernel,
) -> [f64; 3] {
    let norm = kernel_volume_norm(kernel);
    let nq = 128usize;
    let dq = 2.0 / nq as f64;
    let mut ik = [0.0_f64; 5];
    for iq in 0..nq {
        let q = (iq as f64 + 0.5) * dq;
        let w = sph_kernel::kernel_and_dh(q, 1.0, kernel).0 / norm;
        let mu_max = if q <= a_true_h { 1.0 } else { (a_true_h / q).clamp(-1.0, 1.0) };
        for (k, slot) in ik.iter_mut().enumerate() {
            let ku = k as i32;
            let sign = if (k + 1) % 2 == 0 { 1.0 } else { -1.0 };
            let ang = (mu_max.powi(ku + 1) - sign) / (k as f64 + 1.0);
            *slot += 2.0 * PI * q.powi(ku + 2) * w * ang * dq;
        }
    }
    [
        c_scale * (ik[0] + beta_slope * ik[1] + gamma_curv * ik[2]),
        c_scale * (ik[1] + beta_slope * ik[2] + gamma_curv * ik[3]),
        c_scale * (ik[2] + beta_slope * ik[3] + gamma_curv * ik[4]),
    ]
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThreeMomentCandidate {
    /// Affine enthalpy blend with H=0 ghost at three-moment depth `a*h`
    /// using scale-normalized weights `w_ij / c_i` (falls back to original
    /// scalar density whenever the three-moment plane is unidentified).
    ThreeMomentGhostBlend,
    /// Affine enthalpy blend with H=0 ghost at `a*h` AND locally de-biased
    /// neighbor weights `w_ij / (c_i * (1 + beta_i * z_ij))`.
    ThreeMomentDebiasedGhostBlend,
    /// Locally de-biased affine enthalpy blend without a vacuum ghost,
    /// isolating volume-gradient correction from the H=0 constraint.
    ThreeMomentDebiasedAffineBlend,
}

impl ThreeMomentCandidate {
    pub fn label(self) -> &'static str {
        match self {
            Self::ThreeMomentGhostBlend => "three_moment_ghost_blend",
            Self::ThreeMomentDebiasedGhostBlend => "three_moment_debiased_ghost_blend",
            Self::ThreeMomentDebiasedAffineBlend => "three_moment_debiased_affine_blend",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThreeMomentGeometry {
    pub m0: f64,
    pub m1: f64,
    pub m2: f64,
    pub first_moment_vec: V3,
    pub first_moment_norm: f64,
    pub outward: V3,
    pub normal_reliable: bool,
    pub status: ThreeMomentStatus,
    pub physical_root_count: usize,
    pub inferred_depth_h: Option<f64>,
    pub inferred_c: Option<f64>,
    pub inferred_beta: Option<f64>,
    pub df_da: Option<f64>,
    pub depth_sensitivity: Option<f64>,
    pub ideal_support_at_depth: Option<f64>,
    pub ghost_weight: f64,
    pub blend_weight: f64,
}

#[derive(Clone, Debug)]
pub struct ThreeMomentField {
    pub mode: ThreeMomentCandidate,
    pub rho: Vec<f64>,
    pub unidentified_fraction: f64,
    pub fit_fallback_fraction: f64,
    pub maximum_pivot_spread: f64,
}

#[derive(Clone, Debug)]
pub struct ThreeMomentBatch {
    pub fields: Vec<ThreeMomentField>,
    pub geometry: Vec<ThreeMomentGeometry>,
}

#[derive(Clone, Copy, Debug)]
struct WeightedSample {
    q: V3,
    z: f64,
    w: f64,
    h: f64,
}

fn basis(q: V3) -> [f64; 4] {
    [1.0, q.x, q.y, q.z]
}

fn weighted_affine(samples: &[WeightedSample], ghost: Option<(V3, f64)>) -> Option<([f64; 4], f64)> {
    if samples.len() < 6 {
        return None;
    }
    let mut a = [[0.0_f64; 4]; 4];
    let mut b = [0.0_f64; 4];
    for s in samples {
        let v = basis(s.q);
        for i in 0..4 {
            b[i] += s.w * v[i] * s.h;
            for j in 0..4 {
                a[i][j] += s.w * v[i] * v[j];
            }
        }
    }
    if let Some((q, w)) = ghost {
        let v = basis(q);
        for i in 0..4 {
            for j in 0..4 {
                a[i][j] += w * v[i] * v[j];
            }
        }
    }
    let mut scale = [0.0; 4];
    for i in 0..4 {
        if a[i][i] <= 0.0 || !a[i][i].is_finite() {
            return None;
        }
        scale[i] = a[i][i].sqrt();
    }
    let mut l = [[0.0_f64; 4]; 4];
    let (mut min_pivot, mut max_pivot) = (f64::INFINITY, 0.0_f64);
    for i in 0..4 {
        for j in 0..=i {
            let mut value = a[i][j] / (scale[i] * scale[j]);
            for k in 0..j {
                value -= l[i][k] * l[j][k];
            }
            if i == j {
                if !value.is_finite() || value <= 1e-10 {
                    return None;
                }
                min_pivot = min_pivot.min(value);
                max_pivot = max_pivot.max(value);
                l[i][i] = value.sqrt();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    let spread = max_pivot / min_pivot;
    if !spread.is_finite() || spread > 1e8 {
        return None;
    }
    let mut y = [0.0; 4];
    for i in 0..4 {
        let mut value = b[i] / scale[i];
        for j in 0..i {
            value -= l[i][j] * y[j];
        }
        y[i] = value / l[i][i];
    }
    let mut beta = [0.0; 4];
    for i in (0..4).rev() {
        let mut value = y[i];
        for j in i + 1..4 {
            value -= l[j][i] * beta[j];
        }
        beta[i] = value / l[i][i];
    }
    for i in 0..4 {
        beta[i] /= scale[i];
    }
    if !beta.iter().all(|v| v.is_finite()) {
        return None;
    }
    Some((beta, spread))
}

fn smooth_gate_from_ideal_support(ideal_s: f64) -> f64 {
    let t = (2.0 * (1.0 - ideal_s)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub fn reconstruct_three_moment(
    p: &[Fluid3],
    state: &AdaptiveState,
    kernel: DensityKernel,
    eos: DegenerateElectrons,
) -> ThreeMomentBatch {
    let n = p.len();
    assert!(n >= 16 && state.valid_for(n));
    let rho = &state.densities_kg_m3;
    assert!(rho.iter().all(|x| x.is_finite() && *x > 0.0));
    let enthalpy: Vec<f64> = rho.iter().map(|&x| eos.enthalpy_m2_s2(x)).collect();
    let table = MomentTable::new(kernel);
    let max_h = state.smoothing_m.iter().copied().fold(0.0_f64, f64::max);
    let candidates = CellNeighbors::new(p, max_h).all_candidates();

    let mut geometry = Vec::with_capacity(n);
    let mut ghost_blend = Vec::with_capacity(n);
    let mut debiased_ghost_blend = Vec::with_capacity(n);
    let mut debiased_affine_blend = Vec::with_capacity(n);

    let mut unidentified = 0usize;
    let (mut fail_ghost, mut fail_debiased_ghost, mut fail_debiased_affine) = (0usize, 0usize, 0usize);
    let (mut spread_ghost, mut spread_debiased_ghost, mut spread_debiased_affine) =
        (0.0_f64, 0.0_f64, 0.0_f64);

    for i in 0..n {
        let h = state.smoothing_m[i];
        let mut raw_samples: Vec<(V3, f64, f64)> = Vec::new();
        let mut m0 = 0.0;
        let mut first_vec = V3::default();
        for &j in &candidates[i] {
            let d = p[j].pos_m - p[i].pos_m;
            let r = d.len();
            if r >= 2.0 * h {
                continue;
            }
            let w = p[j].mass_kg / rho[j] * sph_kernel::kernel_and_dh(r, h, kernel).0;
            if w <= 0.0 {
                continue;
            }
            let q = d * (1.0 / h);
            m0 += w;
            first_vec += q * w;
            raw_samples.push((q, w, enthalpy[j]));
        }
        assert!(m0 > 0.0);
        let q_norm = first_vec.len();
        let normal_reliable = q_norm > 1e-8 * m0;
        let outward = if normal_reliable {
            first_vec * (-1.0 / q_norm)
        } else {
            V3::default()
        };

        let (mut m1, mut m2) = (0.0_f64, 0.0_f64);
        let (mut z_min, mut z_max) = (0.0_f64, 0.0_f64);
        let mut samples = Vec::with_capacity(raw_samples.len());
        for &(q, w, hj) in &raw_samples {
            let z = if normal_reliable { outward.dot(q) } else { 0.0 };
            m1 += w * z;
            m2 += w * z * z;
            z_min = z_min.min(z);
            z_max = z_max.max(z);
            samples.push(WeightedSample { q, z, w, h: hj });
        }

        let fit = if normal_reliable {
            solve_three_moments_with_table(
                m0,
                m1,
                m2,
                Some((z_min, z_max)),
                &table,
                DEFAULT_MIN_DF_DA,
            )
        } else {
            ThreeMomentFit {
                status: ThreeMomentStatus::InvalidMoments,
                root: None,
                physical_root_count: 0,
                r1: 0.0,
                r2: 0.0,
            }
        };

        if !fit.status.is_identified() {
            unidentified += 1;
            geometry.push(ThreeMomentGeometry {
                m0,
                m1,
                m2,
                first_moment_vec: first_vec,
                first_moment_norm: q_norm,
                outward,
                normal_reliable,
                status: fit.status,
                physical_root_count: fit.physical_root_count,
                inferred_depth_h: fit.root.map(|r| r.a),
                inferred_c: fit.root.map(|r| r.c),
                inferred_beta: fit.root.map(|r| r.beta),
                df_da: fit.root.map(|r| r.df_da),
                depth_sensitivity: fit.root.map(|r| r.depth_sensitivity),
                ideal_support_at_depth: None,
                ghost_weight: 0.0,
                blend_weight: 0.0,
            });
            // Stage A.5 invariant: preserve fallback to original scalar density
            // whenever unidentified or poorly conditioned; never impose H=0 ghost.
            ghost_blend.push(rho[i]);
            debiased_ghost_blend.push(rho[i]);
            debiased_affine_blend.push(rho[i]);
            continue;
        }

        let root = fit.root.unwrap();
        let ideal_s = table.eval(root.a).i[0];
        let ghost_w = (1.0 - ideal_s).clamp(0.0, 0.5);
        let alpha = smooth_gate_from_ideal_support(ideal_s);

        geometry.push(ThreeMomentGeometry {
            m0,
            m1,
            m2,
            first_moment_vec: first_vec,
            first_moment_norm: q_norm,
            outward,
            normal_reliable,
            status: fit.status,
            physical_root_count: fit.physical_root_count,
            inferred_depth_h: Some(root.a),
            inferred_c: Some(root.c),
            inferred_beta: Some(root.beta),
            df_da: Some(root.df_da),
            depth_sensitivity: Some(root.depth_sensitivity),
            ideal_support_at_depth: Some(ideal_s),
            ghost_weight: ghost_w,
            blend_weight: alpha,
        });

        let scaled_samples: Vec<WeightedSample> = samples
            .iter()
            .map(|s| WeightedSample {
                w: s.w / root.c,
                ..*s
            })
            .collect();
        let debiased_samples: Vec<WeightedSample> = samples
            .iter()
            .map(|s| WeightedSample {
                w: s.w / (root.c * (1.0 + root.beta * s.z)),
                ..*s
            })
            .collect();

        let ghost_opt = if ghost_w > 0.0 {
            Some((outward * root.a, ghost_w))
        } else {
            None
        };

        let to_rho = |res: Option<([f64; 4], f64)>| {
            res.filter(|(b, _)| b[0] > 0.0)
                .map(|(b, sp)| (eos.density_from_enthalpy(b[0]), sp))
                .filter(|(val, _)| val.is_finite() && *val > 0.0)
        };

        let fit_g = to_rho(weighted_affine(&scaled_samples, ghost_opt));
        let fit_dg = to_rho(weighted_affine(&debiased_samples, ghost_opt));
        let fit_da = to_rho(weighted_affine(&debiased_samples, None));

        if fit_g.is_none() {
            fail_ghost += 1;
        }
        if fit_dg.is_none() {
            fail_debiased_ghost += 1;
        }
        if fit_da.is_none() {
            fail_debiased_affine += 1;
        }
        if let Some((_, sp)) = fit_g {
            spread_ghost = spread_ghost.max(sp);
        }
        if let Some((_, sp)) = fit_dg {
            spread_debiased_ghost = spread_debiased_ghost.max(sp);
        }
        if let Some((_, sp)) = fit_da {
            spread_debiased_affine = spread_debiased_affine.max(sp);
        }

        let rho_g = fit_g
            .map(|(v, _)| (1.0 - alpha) * rho[i] + alpha * v)
            .unwrap_or(rho[i]);
        let rho_dg = fit_dg
            .map(|(v, _)| (1.0 - alpha) * rho[i] + alpha * v)
            .unwrap_or(rho[i]);
        let rho_da = fit_da
            .map(|(v, _)| (1.0 - alpha) * rho[i] + alpha * v)
            .unwrap_or(rho[i]);

        ghost_blend.push(rho_g);
        debiased_ghost_blend.push(rho_dg);
        debiased_affine_blend.push(rho_da);
    }

    let unid_frac = unidentified as f64 / n as f64;
    ThreeMomentBatch {
        fields: vec![
            ThreeMomentField {
                mode: ThreeMomentCandidate::ThreeMomentGhostBlend,
                rho: ghost_blend,
                unidentified_fraction: unid_frac,
                fit_fallback_fraction: fail_ghost as f64 / n as f64,
                maximum_pivot_spread: spread_ghost,
            },
            ThreeMomentField {
                mode: ThreeMomentCandidate::ThreeMomentDebiasedGhostBlend,
                rho: debiased_ghost_blend,
                unidentified_fraction: unid_frac,
                fit_fallback_fraction: fail_debiased_ghost as f64 / n as f64,
                maximum_pivot_spread: spread_debiased_ghost,
            },
            ThreeMomentField {
                mode: ThreeMomentCandidate::ThreeMomentDebiasedAffineBlend,
                rho: debiased_affine_blend,
                unidentified_fraction: unid_frac,
                fit_fallback_fraction: fail_debiased_affine as f64 / n as f64,
                maximum_pivot_spread: spread_debiased_affine,
            },
        ],
        geometry,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_adaptive::{update_adaptive_densities, AdaptiveConfig};
    use crate::sph_glass;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::profile_for_mass;

    #[test]
    fn exact_1d_quadratures_match_analytic_limits_and_leibniz_derivatives() {
        let m4_0 = planar_moments(0.0, DensityKernel::CubicM4);
        let m4_2 = planar_moments(2.0, DensityKernel::CubicM4);
        assert!((m4_0[0] - 0.5).abs() < 1e-15);
        assert!((m4_0[1] - (-31.0 / 140.0)).abs() < 1e-15);
        assert!((m4_0[2] - 0.15).abs() < 1e-15);
        assert!((m4_0[3] - (-127.0 / 1008.0)).abs() < 1e-15);
        assert!((m4_2[0] - 1.0).abs() < 1e-15);
        assert!(m4_2[1].abs() < 1e-16);
        assert!((m4_2[2] - 0.3).abs() < 1e-15);
        assert!(m4_2[3].abs() < 1e-16);

        let c4_0 = planar_moments(0.0, DensityKernel::WendlandC4);
        let c4_2 = planar_moments(2.0, DensityKernel::WendlandC4);
        assert!((c4_0[0] - 0.5).abs() < 1e-15);
        assert!((c4_0[1] - (-35.0 / 192.0)).abs() < 1e-15);
        assert!((c4_0[2] - (4.0 / 39.0)).abs() < 1e-15);
        assert!((c4_0[3] - (-15.0 / 208.0)).abs() < 1e-15);
        assert!((c4_2[0] - 1.0).abs() < 1e-15);
        assert!(c4_2[1].abs() < 1e-16);
        assert!((c4_2[2] - (8.0 / 39.0)).abs() < 1e-15);
        assert!(c4_2[3].abs() < 1e-16);

        for kernel in [DensityKernel::CubicM4, DensityKernel::WendlandC4] {
            assert!((kernel_volume_norm(kernel) - 1.0).abs() < 1e-14);
            let eps = 1e-6;
            for a in [0.05, 0.35, 0.85, 1.05, 1.55] {
                let mid = planar_moments_and_deriv(a, kernel);
                let plus = planar_moments(a + eps, kernel);
                let minus = planar_moments(a - eps, kernel);
                assert!(mid.i[0] * mid.i[2] - mid.i[1] * mid.i[1] > 0.0);
                for k in 0..4 {
                    let fd = (plus[k] - minus[k]) / (2.0 * eps);
                    assert!(
                        (mid.di_da[k] - fd).abs() < 1e-9,
                        "{kernel:?} a={a} k={k}: deriv={} vs fd={fd}",
                        mid.di_da[k]
                    );
                }
            }
        }
    }

    #[test]
    fn recovers_synthetic_plane_depths_and_volume_slopes_while_q_over_s_fails() {
        for kernel in [DensityKernel::CubicM4, DensityKernel::WendlandC4] {
            let table = MomentTable::new(kernel);
            for &a_true in &[0.0, 0.1, 0.25, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0] {
                let [i0, i1, i2, i3] = table.eval(a_true).i;
                for &c_true in &[0.5, 0.85, 1.0, 1.25, 1.5] {
                    for &b_true in &[-0.4, -0.2, 0.0, 0.2, 0.4] {
                        let m0 = c_true * (i0 + b_true * i1);
                        let m1 = c_true * (i1 + b_true * i2);
                        let m2 = c_true * (i2 + b_true * i3);
                        let exact = solve_three_moments_with_table(m0, m1, m2, None, &table, 0.0);
                        assert_eq!(
                            exact.physical_root_count, 1,
                            "{kernel:?} a={a_true} c={c_true} b={b_true}"
                        );
                        let r = exact.root.unwrap();
                        assert!((r.a - a_true).abs() < 1e-8, "a mismatch: {} vs {a_true}", r.a);
                        assert!((r.c - c_true).abs() < 1e-7, "c mismatch: {} vs {c_true}", r.c);
                        assert!(
                            (r.beta - b_true).abs() < 1e-6,
                            "beta mismatch: {} vs {b_true}",
                            r.beta
                        );

                        if a_true <= 1.2 {
                            let conditioned = solve_three_moments_with_table(
                                m0,
                                m1,
                                m2,
                                None,
                                &table,
                                DEFAULT_MIN_DF_DA,
                            );
                            assert!(
                                conditioned.status.is_identified(),
                                "expected well-conditioned surface root at {kernel:?} a={a_true} b={b_true}"
                            );
                        }
                    }
                }
            }
            let a_ref = 0.5;
            let [i0, i1, i2, i3] = table.eval(a_ref).i;
            let (c_ref, b_ref) = (0.85, -0.35);
            let m0 = c_ref * (i0 + b_ref * i1);
            let m1 = c_ref * (i1 + b_ref * i2);
            let m2 = c_ref * (i2 + b_ref * i3);
            let fit = solve_three_moments_with_table(m0, m1, m2, None, &table, DEFAULT_MIN_DF_DA);
            assert!(fit.status.is_identified());
            assert!((fit.root.unwrap().a - a_ref).abs() < 1e-11);
            let biased_q_over_s = m1.abs() / m0;
            let ideal_q_over_s = i1.abs() / i0;
            assert!((biased_q_over_s - ideal_q_over_s).abs() > 0.04);
        }
    }

    #[test]
    fn spherical_curvature_and_quadratic_volume_bias_shift_planar_depth_inward() {
        for kernel in [DensityKernel::CubicM4, DensityKernel::WendlandC4] {
            let table = MomentTable::new(kernel);
            // 1. Spherical curvature: a sphere of radius R=3h curves inward below its tangent plane,
            // causing a planar estimator to underestimate radial depth a_true=0.45h, and converging
            // monotonically toward a_true as R/h -> infinity.
            let a_true = 0.45;
            let m_tight = spherical_cap_moments(a_true, 2.8, 1.0, -0.2, kernel);
            let m_wide = spherical_cap_moments(a_true, 25.0, 1.0, -0.2, kernel);
            let fit_tight = solve_three_moments_with_table(
                m_tight[0], m_tight[1], m_tight[2], None, &table, DEFAULT_MIN_DF_DA,
            );
            let fit_wide = solve_three_moments_with_table(
                m_wide[0], m_wide[1], m_wide[2], None, &table, DEFAULT_MIN_DF_DA,
            );
            assert!(fit_tight.status.is_identified() && fit_wide.status.is_identified());
            let a_tight = fit_tight.root.unwrap().a;
            let a_wide = fit_wide.root.unwrap().a;
            assert!(
                a_tight < a_wide && (a_wide - a_true).abs() < 0.02,
                "{kernel:?}: tight={a_tight}, wide={a_wide}, true={a_true}"
            );

            // 2. Nonlinear quadratic volume bias gamma != 0 confounds the 3-moment linear-bias model.
            let m_quad = quadratic_volume_bias_moments(0.45, 0.75, -0.35, -0.18, kernel);
            let fit_quad = solve_three_moments_with_table(
                m_quad[0], m_quad[1], m_quad[2], None, &table, DEFAULT_MIN_DF_DA,
            );
            if let Some(r) = fit_quad.root {
                assert!(
                    (r.a - 0.45).abs() > 0.02,
                    "{kernel:?}: quadratic volume curvature must shift linear-bias depth"
                );
            }
        }
    }

    #[test]
    fn constant_affine_and_curved_enthalpy_profiles_with_debiased_weights() {
        let table = MomentTable::new(DensityKernel::WendlandC4);
        let a_true = 0.55;
        let c_true = 1.25;
        let beta_true = -0.30;
        let step = 0.22_f64;
        let d3 = step.powi(3);
        let mut raw = Vec::new();
        let (mut m0, mut m1, mut m2) = (0.0, 0.0, 0.0);
        let mut idx = 0_u64;
        let mut xi = -1.98;
        while xi <= 1.98 {
            let mut yi = -1.98;
            while yi <= 1.98 {
                let mut zi = -1.98;
                while zi <= a_true {
                    idx = idx.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    let jx = (((idx >> 33) as f64 / (1u64 << 31) as f64) - 0.5) * 0.04;
                    idx = idx.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    let jy = (((idx >> 33) as f64 / (1u64 << 31) as f64) - 0.5) * 0.04;
                    let q = V3::new(xi + jx, yi + jy, zi);
                    let r = q.len();
                    if r < 2.0 && q.z <= a_true {
                        let w_true = d3 * sph_kernel::kernel_and_dh(r, 1.0, DensityKernel::WendlandC4).0;
                        let w_est = c_true * (1.0 + beta_true * q.z) * w_true;
                        m0 += w_est;
                        m1 += w_est * q.z;
                        m2 += w_est * q.z * q.z;
                        raw.push((q, w_est));
                    }
                    zi += step;
                }
                yi += step;
            }
            xi += step;
        }
        let fit = solve_three_moments_with_table(m0, m1, m2, Some((-1.98, a_true)), &table, DEFAULT_MIN_DF_DA);
        assert!(fit.status.is_identified());
        let root = fit.root.unwrap();
        assert!((root.a - a_true).abs() < 0.03, "irregular cloud depth: {} vs {a_true}", root.a);
        assert!((root.beta - beta_true).abs() < 0.06, "irregular cloud beta: {} vs {beta_true}", root.beta);

        let ideal_s = table.eval(root.a).i[0];
        let ghost_w = (1.0 - ideal_s).clamp(0.0, 0.5);
        let n_out = V3::new(0.0, 0.0, 1.0);

        let const_samples: Vec<WeightedSample> = raw
            .iter()
            .map(|&(q, w)| WeightedSample {
                q,
                z: q.z,
                w: w / (root.c * (1.0 + root.beta * q.z)),
                h: 12.0,
            })
            .collect();
        let (c_fit, _) = weighted_affine(&const_samples, None).unwrap();
        assert!((c_fit[0] - 12.0).abs() < 1e-11);

        let affine_samples: Vec<WeightedSample> = raw
            .iter()
            .map(|&(q, w)| WeightedSample {
                q,
                z: q.z,
                w: w / (root.c * (1.0 + root.beta * q.z)),
                h: 10.0 * (a_true - q.z),
            })
            .collect();
        let (a_fit, _) = weighted_affine(&affine_samples, Some((n_out * root.a, ghost_w))).unwrap();
        assert!((a_fit[0] - 10.0 * a_true).abs() < 0.15);

        let curved_samples: Vec<WeightedSample> = raw
            .iter()
            .map(|&(q, w)| {
                let d = a_true - q.z;
                WeightedSample {
                    q,
                    z: q.z,
                    w: w / (root.c * (1.0 + root.beta * q.z)),
                    h: 10.0 * d + 2.5 * d * d,
                }
            })
            .collect();
        let (q_fit, _) = weighted_affine(&curved_samples, Some((n_out * root.a, ghost_w))).unwrap();
        let h_true = 10.0 * a_true + 2.5 * a_true * a_true;
        assert!(q_fit[0].is_finite() && q_fit[0] > 0.0 && (q_fit[0] - h_true).abs() < 0.5 * h_true);
    }

    #[test]
    fn impossible_nonplanar_and_ill_conditioned_cases_report_unidentified() {
        for kernel in [DensityKernel::CubicM4, DensityKernel::WendlandC4] {
            let table = MomentTable::new(kernel);
            let impossible = solve_three_moments_with_table(0.7, -0.15, 0.85 * 0.7, None, &table, DEFAULT_MIN_DF_DA);
            assert_eq!(impossible.status, ThreeMomentStatus::NoPhysicalRoot);
            assert!(impossible.root.is_none());

            let bad_cs = solve_three_moments_with_table(1.0, -0.6, 0.2, None, &table, DEFAULT_MIN_DF_DA);
            assert_eq!(bad_cs.status, ThreeMomentStatus::InvalidMoments);

            let [i0, i1, i2, i3] = table.eval(1.98).i;
            let deep = solve_three_moments_with_table(
                i0 - 0.1 * i1,
                i1 - 0.1 * i2,
                i2 - 0.1 * i3,
                None,
                &table,
                DEFAULT_MIN_DF_DA,
            );
            assert_eq!(deep.status, ThreeMomentStatus::IllConditionedInterior);
            assert!(!deep.status.is_identified());
        }
    }

    #[test]
    fn white_dwarf_surface_identifiability_and_score_only_depth_audit() {
        let star = profile_for_mass(
            0.65 * SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),
            10_000.0,
        )
        .unwrap();
        let center = V3::new(1e8, 0.0, 0.0);
        let base = seed_cubic_shells_oriented(&star, center, V3::default(), 96, 17).0;
        let p = sph_glass::relax_antipodal_angles(&base, &star, center, 1.2, 20).0;
        for (kernel, eta) in [(DensityKernel::CubicM4, 1.2), (DensityKernel::WendlandC4, 1.7)] {
            let mut cfg = AdaptiveConfig::from_star_radius(star.radius_m)
                .with_cubic_gravity()
                .with_newton_density();
            if kernel == DensityKernel::WendlandC4 {
                cfg = cfg.with_wendland_c4();
            }
            cfg.eta = eta;
            let mut state = AdaptiveState::new(&p, star.radius_m, cfg);
            update_adaptive_densities(&p, &mut state, cfg);
            let batch = reconstruct_three_moment(&p, &state, kernel, star.eos);
            for f in &batch.fields {
                assert!(f.rho.iter().all(|v| v.is_finite() && *v > 0.0));
            }
            // Unidentified particles must strictly retain exact original scalar density.
            for (idx, g) in batch.geometry.iter().enumerate() {
                if !g.status.is_identified() {
                    for f in &batch.fields {
                        assert_eq!(f.rho[idx], state.densities_kg_m3[idx]);
                    }
                }
            }
        }
    }
}
