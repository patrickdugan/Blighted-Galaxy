# BG-PHY-20: Free-surface density reconstruction and enthalpy bias

**Date:** 10 November 2026  
**Predecessor:** [BG-PHY-19](QUADRATIC_ENTHALPY_SURFACE.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Research question

BG-PHY-19 established that the independent hydrostatic **enthalpy**
field is recoverable with a local polynomial on the stellar seed.
For N1024, C4 angular glass, 200 neighbors, its quadratic reference
enthalpy gradient has dimensionless force-vector error 0.0716. The
same estimator fed *numerical* SPH enthalpy has error 0.1561,
and an independent density audit measured mean outer-quartile
log(rho_SPH/rho_ref) of +0.242. At N2048, those enthalpy errors
are 0.0454 (reference) and 0.1113 (numerical).

The present experiment tests the density-sampling hypothesis
directly, while keeping the particle geometry and gradient
operator fixed. A method must reduce **both** the outer density
log-error and the numerical quadratic enthalpy-gradient error
at matched resolution, kernel and sample positions to qualify
as an improvement.

## Candidate estimators without reference-profile labels

For each particle i, the baseline density and volumes are

```text
rho_i = sum_j m_j W(|x_i-x_j|, h_i),
V_j = m_j/rho_j.
```

Six density fields are scored:

1. **Original SPH**: exactly the prior density values.
2. **Shepard normalization**:
   ```text
   S_i = sum_j V_j W_ij(h_i),
   rho_i,Shepard = (sum_j m_j W_ij)/S_i.
   ```
   This is a partition-of-unity correction. Its validity near a
   vacuum stellar surface is an empirical hypothesis.
3. **Local affine MLS**: weighted least-squares fit to
   existing rho_j with the basis [1, dx/h, dy/h, dz/h].
   The fitted intercept estimates the central particle density.
4. **Local quadratic MLS**: expand the same basis with
   [dx²/2, dy²/2, dz²/2, dx*dy, dx*dz, dy*dz] / h².
   The normal equations use diagonal preconditioning and Cholesky,
   positive pivots and a maximum pivot-spread proxy of 1e8.
5. **Surface-blended affine MLS (exploratory)**: use the observed
   missing partition of unity, S_i, as a geometry-only proxy for
   surface exposure. A locally planar vacuum interface passes
   approximately half the compact kernel's support, motivating:
   ```text
   x_i = clamp((1-S_i)/0.5, 0, 1)
   alpha_i = 3*x_i² - 2*x_i³
   rho_blended,i = (1-alpha_i)*rho_SPH,i + alpha_i*rho_MLS1,i.
   ```
   This is a post-initial-matrix hypothesis motivated by the N256
   result that affine MLS improved outer density while sometimes
   hurting the global enthalpy gradient. The mapping is fixed for
   all subsequent N and orientation seeds, and the candidate is
   compared against both pure endpoints.
6. **Spherical radial mass secant**: sort particles by radius
   about the mass centroid and estimate rho ~ Delta M/Delta V
   using a local rank window of approximately sqrt(N)/2.
   This uses known *particle masses and positions* but explicitly
   assumes spherical symmetry. It is a reference diagnostic for an
   isolated equilibrium star and is unsuitable as an arbitrary
   three-dimensional disruption density estimator.

All estimators read only the frozen particle data, baseline SPH
density and kernel. The independently integrated 1D hydrostatic
density appears **only in scoring**, never as estimator input.
The initial seed's radial quantiles already follow the reference
mass CDF by construction. This could favor the radial-mass secant
on the initial configuration. To expose that dependency, a
separate N256 pilot applies **3 seconds of ordinary M4 damping**
before freezing the particle positions and rerunning every estimator.

The simulation's particle masses and conserved pair forces are never
modified. At any particle where MLS is singular or predicts nonpositive
density, the experiment explicitly falls back to the original
rho_i and records the fallback fraction. No acceptance is awarded to
a method with poor fit coverage.

## Hold the enthalpy-gradient experiment fixed

BG-PHY-19's `sph_quadratic_audit` now accepts an optional frozen
candidate density field. This replaces the *numerical* pressure
and enthalpy values and local density used for dividing the
pressure-gradient estimators. It leaves:

- all particle coordinates and masses unchanged;
- the exact-root smoothing lengths h_i unchanged;
- each local polynomial's original volume weight V_j=m_j/rho_SPH,j
  and kernel unchanged;
- the independent reference pressure, enthalpy and gravity fields
  unchanged;
- all conservative hydrodynamic and gravity force calculations unchanged.

A unit test verifies that supplying the original density reproduces
the previous results exactly, and that a changed candidate cannot
alter the reference-enthalpy oracle results.

This identifies the effect of the candidate density samples on the
gradient field, conditional on the fixed SPH volume weights.
Any future fully self-consistent kernel normalization will require
a separate, conservative derivation.

## Metrics and population

For each candidate, evaluate

```text
L_i = ln(rho_candidate,i / rho_ref(|x_i-C|))
bias_log = sum_i m_i L_i / sum_i m_i
rms_log  = sqrt(sum_i m_i L_i² / sum_i m_i).
```

Report these scores globally and in the inner and outer mass
quartiles, together with the arithmetic mean density ratio,
the Shepard partition sum S_i, GLS/MLS fallback and condition
diagnostics, realized neighbor population, and unchanged
BG-PHY-19 enthalpy-gradient scores. Crucially the **quadratic
reference enthalpy-gradient error** is an internal control that
must remain unchanged for all candidate fields.

The primary frozen matrix spans N256, N1024, N2048;
shell-orientation seeds 0/17 at N<=1024, seed0 at N2048;
original cubic-shell and angular-glass layouts;
matched M4 eta1.20 and C4 eta1.20, then C4 calibrated to
100,150,200 realized neighbors with safeguarded Newton roots.

A separate N256 seed0 experiment repeats the entire frozen
measurement after identical 300-step, dt=0.01s M4 damped
preparation for each geometric layout. The physical acceptance
criteria in Issue #8 are unchanged, and these 3s snapshots
alone cannot pass the 3-dynamical-time undamped hold.

## Reproduction

```bash
cargo test --all-targets

cargo run --release --bin phy20_bench -- \
  --n 1024 --seeds 0,17 --targets 100,150,200 \
  --glass-sweeps 20 --calibration-steps 12 \
  --csv bg-phy20-n1024.csv

cargo run --release --bin phy20_bench -- \
  --n 256 --seeds 0 --targets 100,150,200 \
  --relax-steps 300 --dt 0.01 --csv bg-phy20-n256-damped.csv

# Light unit and numerical smoke, independent of physical acceptance.
cargo run --release --bin phy20_bench -- \
  --n 64 --seeds 0 --targets 20,32 \
  --glass-sweeps 4 --calibration-steps 6 \
  --csv bg-phy20-smoke.csv
```

The [dedicated BG-PHY-20 workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy20.yml)
publishes the frozen matrices and damped pilot as separate CSV
artifacts; ordinary CI runs N64 and Windows release compilation.

## Decision

Rank by *outer-quartile log-density RMS* and numerical
quadratic enthalpy-vector error, requiring improvements across
both seed orientations. Check that the reference enthalpy-gradient
oracle remains invariant under every candidate intervention.

The Shepard correction can worsen density estimates when
SPH kernel smoothing already mixes a high-density interior
into the low-density stellar surface. Likewise the radial
mass secant is a strong symmetry-dependent comparator on the
initial exact-CDF seed, which is why the post-damping pilot is
necessary before recommending it for initialization.

No corrected density field is permitted to enter the live
pressure force or variable-h closure until the associated
conservative momentum/energy equations and support-correction
consistency are derived and tested. Frontiere, Raskin and Owen's
conservative reproducing-kernel SPH provides a separate
architectural starting point.

## Literature

- Frontiere, Raskin and Owen (2017), *CRKSPH: A Conservative
  Reproducing Kernel Smoothed Particle Hydrodynamics Scheme*,
  J. Comput. Phys. 332, 160–209, DOI 10.1016/j.jcp.2016.12.004.
- *The Mathematics of Smoothed Particle Hydrodynamics (SPH)
  Consistency* (2021), Frontiers in Applied Mathematics and
  Statistics, DOI 10.3389/fams.2021.797455.
- Bonet and Lok (1999), *Variational and momentum preservation
  aspects of Smooth Particle Hydrodynamic formulations*,
  DOI 10.1016/S0045-7825(99)00051-1.

## Completed matrix: geometry-gated MLS survives the density and enthalpy tests

The [final BG-PHY-20 four-job matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42797977435) completed N256 (two independent seed orientations), N1024 (two orientations), and N2048 (seed 0), plus a separate N256 case frozen after the **same ordinary M4 3-second damping**. All numerical jobs succeeded and emitted CSV artifacts with density, kernel partition, fallback and gradient diagnostics.

**Fixed configuration: angular glass, Wendland C4, exact-root 200 realized neighbors, orientation seed 0.**

| Candidate, metric | N256 initial | N1024 initial | N2048 initial |
|---|---:|---:|---:|
| Original SPH outer density log-RMS | 0.5529 | 0.4135 | 0.3524 |
| Surface-gated MLS outer density log-RMS | **0.1762** | **0.1736** | **0.1553** |
| Original SPH global density log-RMS | 0.3220 | 0.2103 | 0.1776 |
| Surface-gated MLS global density log-RMS | **0.1910** | **0.0951** | **0.0808** |
| Original SPH global quadratic enthalpy-vector error | 0.3654 | 0.1561 | 0.1113 |
| Surface-gated MLS global quadratic enthalpy-vector error | **0.2930** | **0.1369** | **0.0950** |
| Original SPH outer quadratic enthalpy-vector error | 0.2611 | 0.2237 | 0.1857 |
| Surface-gated MLS outer quadratic enthalpy-vector error | **0.2333** | **0.1881** | **0.1349** |

At N1024, orientation seed 17 independently reproduces the improved pair of measurements: outer density log-RMS 0.4131→0.1602, global enthalpy-gradient error 0.1552→0.1412. Both seeds retain the same fitted hydrostatic reference and original SPH volume weights.

**Paired C4-target evaluation:** 36 unique test conditions comprise N256, N1024 and N2048, shell/glass layouts, 100/150/200 realized neighbor targets, independent shell-orientation seeds where available, and the separate N256 post-damping experiments. In *all 36*, the surface-gated MLS candidate reduces each of:
- global density log-RMS (median reduction 40.1%);
- outer density log-RMS (median reduction 54.0%);
- global quadratic numerical-enthalpy gradient error (median reduction 14.9%);
- outer quadratic numerical-enthalpy gradient error (median reduction 26.1%).

The exact **reference-enthalpy oracle** is invariant across all six candidate fields to the precision of the emitted results, and all fits completed with zero fallback. This is a diagnostic comparison using one fixed reference profile across repeated conditions, rather than 36 statistically independent stars or a formal generalization guarantee.

**After 3 seconds M4 damping (N256, seed0, angular glass, C4/200 frozen audit):**

| Measure | Original SPH | Surface-gated MLS |
|---|---:|---:|
| Outer density log-RMS | 0.2907 | **0.0928** |
| Global density log-RMS | 0.2664 | **0.2303** |
| Global quadratic enthalpy-vector error | 0.4014 | **0.3369** |
| Outer quadratic enthalpy-vector error | 0.1921 | **0.1324** |

### Failures of the other estimators

Shepard normalization **worsens all four metrics in all 36 C4 target conditions**. The baseline outer kernel partition average is 0.690 at N1024 C4/200 (0.761 at N2048), so division by that incomplete partition raises already positively biased surface densities. A constant-completeness correction alone does not solve the white-dwarf free-surface problem.

Affine MLS alone reduces outer density log-RMS in all 36 cases, while it reduces global enthalpy-vector error in only 3 of 36. Its global gradient performance commonly deteriorates because a local intercept correction extends into parts of the star with adequate support. Quadratic MLS is worse than the original SPH density globally and at the outer boundary for these cases.

The spherical radial mass secant gives global frozen enthalpy errors 0.0758 at N1024 and 0.0497 at N2048, close to the reference oracle. This comparison is especially favorable on the *initial exact reference mass-CDF seed*. After ordinary M4 damping at N256/C4/200, its outer enthalpy-gradient error grows from the SPH baseline 0.1921 to 0.3216. The spherical symmetry assumption therefore remains unsuitable for general three-dimensional stellar weapon and disruption dynamics.

### Next physical implementation gate

BG-PHY-20 provides a promising *density field* and reproducible acceptance metrics, with a differentiable smoothstep boundary exposure determined by the existing particle partition. A live implementation of rho_blended requires derivatives with respect to both particle positions and smoothing lengths, including the MLS normal-matrix coefficients, partition S_i, and h/rho closure. Simply assigning the blended density inside the previous symmetric pressure-force formula would break its Lagrangian derivation.

A dedicated conservative reproducing-kernel milestone should derive compatible pair forces and thermodynamic work terms, then verify total linear momentum, angular momentum, energy, and translation/rotation invariance on small isolated systems before the same 0.65-Msun white dwarf is re-evaluated under Issue #8's undamped physical gate. Current gameplay dynamics are unchanged.

