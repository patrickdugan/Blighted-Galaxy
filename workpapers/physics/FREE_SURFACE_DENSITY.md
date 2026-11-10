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

Five density fields are scored:

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
5. **Spherical radial mass secant**: sort particles by radius
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
