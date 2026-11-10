# BG-PHY-19: Local quadratic pressure and enthalpy reconstruction

**Date:** 10 November 2026  
**Predecessor:** [BG-PHY-18](PRESSURE_GRADIENT_CONSISTENCY.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Problem

BG-PHY-18 demonstrated that the first-moment inverse corrects
linear fields exactly, but greatly overestimates the derivative of
the **nonlinear** reference white-dwarf pressure near the surface.
At N1024/seed0/glass/C4 with 200 realized neighbors, the existing
conservative pressure-vector RMS was 0.1204, the direct numerical
difference gradient 0.1552, and the linear moment-corrected numerical
gradient 0.3285. The reference-pressure oracle's corresponding
direct/corrected errors were 0.4474 and 1.1851. The linear reference
oracle's corrected outer radial force was 2.48 times the hydrostatic
target despite a well-conditioned local matrix.

The numerical SPH density is also biased upward near the stellar
surface, with outer-quartile mean log(rho_SPH/rho_ref)=+0.242 at
N1024/C4/200, and the first-moment trace is only 0.615.

The next test asks whether a **local quadratic fit**, which accounts
for second-order pressure curvature, can reconstruct the hydrostatic
field more reliably. A second scalar field, specific enthalpy, makes
the near-vacuum gradient better conditioned.

## Thermodynamic enthalpy

For the existing cold barotropic degenerate-electron EOS,

```text
dH = dP/rho,
a_pressure = -(1/rho) grad(P) = -grad(H).
```

The independent 1D white-dwarf reference solves

```text
dH_ref/dr = -G M_ref(<r)/r^2.
```

The code already implements the exact relativistic-electron
specific enthalpy `DegenerateElectrons::enthalpy_m2_s2(rho)`.
Near the vacuum surface this has a finite radial derivative
equal to the hydrostatic gravity. Pressure itself approaches
zero steeply, so differentiating a wide-support pressure fit
and dividing by a small local density can amplify nonlocal
curvature errors.

BG-PHY-19 reconstructs both numerical and reference H fields,
alongside both numerical and reference P fields. It leaves the
existing EOS implementation, pressure force and gravity force
unchanged.

## Weighted local polynomial fits

For each particle i, use normalized offset

```text
q_ij = (x_j-x_i)/h_i.
```

The linear basis is `[qx,qy,qz]`. The quadratic basis contains
`[qx,qy,qz, qx²/2,qy²/2,qz²/2, qx*qy,qx*qz,qy*qz]`.
The scalar-value difference at particle j is approximated by

```text
f_j-f_i ≈ phi(q_ij) dot beta_i.
```

Using exactly the SPH kernel and numerical volumes
`V_j=m_j/rho_SPH,j`, minimize

```text
sum_{j!=i} V_j W_ij(h_i)
    [phi(q_ij)·beta_i -(f_j-f_i)]².
```

For all four scalar fields the design matrix is identical, so a
single preconditioned normal matrix and Cholesky factorization
solve four independent right-hand sides. Coefficients for qx,
qy,qz yield the gradient after division by h_i.

The code checks positive pivots (threshold 1e-10 in normalized
normal-matrix units), a maximum pivot-spread proxy of 1e8,
and enough supported neighbors for the polynomial order.
It excludes particles whose quadratic and linear fits cannot
both be solved; the CSV reports the usable fraction for each
radial quartile. The direct matrix inverse from BG-PHY-18
is independently available for comparison.

The eight predicted *outward* pressure accelerations are:

```text
pressure_num_linear       = -grad(P_SPH)_linear/rho_SPH
pressure_num_quadratic    = -grad(P_SPH)_quadratic/rho_SPH
pressure_ref_linear       = -grad(P_ref)_linear/rho_ref
pressure_ref_quadratic    = -grad(P_ref)_quadratic/rho_ref

enthalpy_num_linear       = -grad(H_SPH)_linear
enthalpy_num_quadratic    = -grad(H_SPH)_quadratic
enthalpy_ref_linear       = -grad(H_ref)_linear
enthalpy_ref_quadratic    = -grad(H_ref)_quadratic
```

Every estimator is compared with the independent hydrostatic
`g_ref(r_i) rhat_i`, and normalized to the same mass-weighted
RMS reference-gravity scale. Global, inner, and outer mass
quartiles each report vector errors and radial pressure/reference
ratios. Pairwise conservation is **not** claimed for any of these
directed fits.

## Experiment and tests

The matrix covers N256 and N1024 with shell orientation seeds
0/17, and N2048 seed0. Each seed is evaluated on both the
original antipodal shell arrangement and the fixed-reference
angular glass. The control is M4 eta=1.20, with C4 at the same
eta and at 100,150,200 realized neighbors calibrated by exact
h/rho Newton roots.

The new `sph_quadratic_audit` verifies exact reproduction of
known quadratic polynomial coefficients on an irregular point
cloud, exact annihilation of constants, and translation invariance
with an actual C4 star fixture. `phy19_bench` records the
original conservative SPH vector-force error and first-order
BG-PHY-18 results alongside the fitted polynomial estimates.

```bash
cargo test --all-targets

cargo run --release --bin phy19_bench -- \
  --n 1024 --seeds 0,17 --targets 100,150,200 \
  --glass-sweeps 20 --calibration-steps 12 \
  --csv bg-phy19-1024.csv

cargo run --release --bin phy19_bench -- \
  --n 64 --seeds 0 --targets 20,32 \
  --glass-sweeps 4 --calibration-steps 6 \
  --csv bg-phy19-smoke.csv
```

The [BG-PHY-19 GitHub Actions matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy19.yml)
publishes each N256/1024/2048 CSV; regular Rust and Windows CI
test the N64 smoke and compiler build.

## Decision

A quadratic oracle improvement near the surface would support
**pressure curvature** as the source of the large BG-PHY-18
first-moment overestimate. A reference-enthalpy reconstruction
close to the independent hydrostatic gravity would confirm that
enthalpy is a numerically better-conditioned field for a
future, separately designed conservative formulation.

Compare the numerical vs reference field at the same order to
quantify persistent SPH density sampling bias. Report any
rank-deficient polynomial fit cases explicitly; no denominator
or physical-gate tuning is permitted to force acceptance.

This is a diagnostic experiment only. No energy/momentum
conservation claim is transferred to the fitted directed
gradients. The original playable forces and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)
physical acceptance remain unchanged.

## References

- Bonet and Lok (1999), *Variational and momentum preservation
  aspects of Smooth Particle Hydrodynamic formulations*,
  https://doi.org/10.1016/S0045-7825(99)00051-1.
- García-Senz, Cabezón and Escartín (2012), *Improving smoothed
  particle hydrodynamics with an integral approach to calculating
  gradients*, https://doi.org/10.1051/0004-6361/201117939.
- Dehnen and Aly (2012), *Improving convergence in smoothed
  particle hydrodynamics simulations without pairing instability*,
  https://doi.org/10.1111/j.1365-2966.2012.21439.x.
