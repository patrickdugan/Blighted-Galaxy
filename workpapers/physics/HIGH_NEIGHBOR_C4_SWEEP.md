# BG-PHY-17: High-neighbor Wendland C4 regime with exact-root calibration

**Date:** 10 November 2026  
**Depends on:** [BG-PHY-16](WENDLAND_GLASS_SEEDS.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Research question

BG-PHY-16 found that the constrained angular-glass initializer lowers
initial pressure-force error, while its Wendland C4 experiments at
approximately 35–80 support neighbors still produce large errors,
especially near the stellar surface. Since the intended benefit of
Wendland kernels is freedom from the high-neighbor pairing instability,
BG-PHY-17 measures their performance at **100, 150, and 200 realized
kernel-support neighbors**, holding the independent hydrostatic
reference, EOS, stellar mass and sampling radii fixed.

Dehnen and Aly (2012) demonstrate that Wendland kernels can support
larger neighbor populations without the pairing instability that limits
traditional B-splines. Zhu, Hernquist and Li (2015) require
simultaneously growing total particle and neighbor counts, with
decreasing smoothing length, for formal SPH convergence. The present
N=1024 and 2048 tests are finite-resolution experiments; no asymptotic
convergence claim follows from one sweep.

## Why calibrate actual neighbors

The uniform-density estimate is

```text
N_neighbors ≈ (32*pi/3)*eta^3
```

for support radius 2h. The white dwarf has a steep radial density
gradient, an outer free surface, inhomogeneous particle separations,
and configured h caps, so this relationship is only a starting
estimate. Every target is therefore calibrated independently on
the actual particle coordinates:

```text
h_i^3 rho_i(h_i) = eta^3 m_i,
rho_i(h) = sum_j m_j W_ij(h_i),
Nbar(eta) = (1/N) sum_i number_of_j[ |r_i-r_j| < 2h_i ].
```

We use the safeguarded Newton h/rho solver from BG-PHY-14 to evaluate
Nbar for each eta. A bounded eta bisection evaluates [0.75, 2.85] in
12 stages and returns the nearest achieved mean, with an explicit
`calibration_bracketed` flag for targets outside the reachable range.
The count includes the particle itself. Actual inner and outer radial
quartile counts, min/max, h caps, mean h, and Omega values are reported.

The optimizer operates on the **exact same physical positions** for
every eta within a layout. Angular glass itself uses a fixed
reference eta=1.20 and 20 monotonic gradient sweeps, so the target
sweep does not re-optimize the seed geometry as eta changes.
Seed=0 and seed=17 remain independent shell-orientation replicates.

## Factorial design

| Resolution | Seed orientations | Layouts | Kernels | Target support counts |
|---:|---|---|---|---|
|1024|0, 17|shell, glass|matched M4, matched C4|100, 150, 200|
|2048|0|shell, glass|matched M4, matched C4|100, 150, 200|

At each combination, record eta=1.20 and the predecessor growth-law
eta(N) as baselines. M4 is a *frozen-state control* at matched actual
neighbor count. It may experience a pairing instability during free
dynamics at these large smoothing lengths, so this matrix deliberately
makes no high-neighbor M4 dynamics claim.

Each frozen particle arrangement is audited using both the original
24-pass underrelaxed density closure and an independent exact Newton
root. The force reports record the same hydrostatic pressure, gravity,
net vector errors, radial force ratios, tangential noise, eight
mass-quantile shells, and the exact initial sampling-radius bias as
BG-PHY-15 and BG-PHY-16.

## New numerical consistency measures

The discrete zeroth-order gradient residual is measured as a
dimensionless per-particle vector:

```text
E0_i = h_i sum_(j != i) [m_j/rho_j] grad_i W(|r_i-r_j|,h_i).
```

Its RMS is split into inner (0–25% enclosed particle mass) and outer
(75–100%) quartiles. The inner tangential RMS distinguishes angular
noise from the radial surface truncation, while the signed outer radial
E0 mean indicates boundary bias.

The isotropic trace of the discrete first-moment gradient matrix is

```text
E1_trace_i = (1/3) sum_j [m_j/rho_j]
              (r_j-r_i) dot grad_i W_ij(h_i).
```

It should approach 1 for linearly consistent interior gradients. The
benchmark reports a global and outer-quartile mean to flag loss of
first-gradient consistency, as well as the pressure-force comparison
against the independent continuum solution.

These E0/E1 observables use the exact same density kernel, smoothing
length and SPH density as the force calculation. They are diagnostic
moments, not a modified pressure force law.

## Reproduction and CI

```bash
cargo test --all-targets

cargo run --release --bin phy17_bench -- \
  --n 1024 --seeds 0,17 --targets 100,150,200 \
  --calibration-steps 12 --glass-sweeps 20 --shells 8 \
  --csv bg-phy17-n1024.csv --shell-csv bg-phy17-n1024-shells.csv

cargo run --release --bin phy17_bench -- \
  --n 2048 --seeds 0 --targets 100,150,200 \
  --calibration-steps 12 --glass-sweeps 20 --shells 8 \
  --csv bg-phy17-n2048.csv --shell-csv bg-phy17-n2048-shells.csv

# Fast smoke, including unreachable-target diagnostics when needed.
cargo run --release --bin phy17_bench -- \
  --n 64 --seeds 0 --targets 20,32 --calibration-steps 6 \
  --glass-sweeps 4 --shells 4 --csv bg-phy17-smoke.csv
```

The dedicated [BG-PHY-17 GitHub Actions matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy17.yml)
runs N1024 and N2048 separately and uploads summary and radial-shell
CSVs. Ordinary CI checks unit tests, a small N64 numerical smoke and
the Windows release binary.

## Follow-up dynamics gate

After the frozen matrix, shortlist C4 choices that satisfy all of:

- exact-root residual below 1e-8 at uncapped particles;
- calibrated mean neighbor count within 5% of the requested target;
- no Omega floor, and less than 5% particles at h caps;
- smaller initial pressure-vector error than the BG-PHY-16 C4
  growth-law baseline, with outer-shell pressure support tracked.

Rank those by pressure-force RMS and inner tangential E0 residual
across both N1024 seed orientations. Explicitly report seed variance
and the effect of halving h at increasing resolution.

Then run the best qualified high-neighbor C4 candidate and the original
M4 control through the *existing* 3.00s constant-tau damped preparation
and at least 7.88s of undamped evolution. The same independent physical
gate from Issue #8 applies: radius and R99 drift, Rmax radius bias,
radial mass-CDF RMS, momentum and energy conservation, and completed
three-dynamical-time hold. A successful force snapshot alone does not
satisfy stellar equilibrium acceptance.

If no C4 candidate qualifies, report the failures and target the
source of the remaining pressure-gradient bias before longer dynamics.

## References

- W. Dehnen and H. Aly, *Improving convergence in smoothed particle
  hydrodynamics simulations without pairing instability*, MNRAS 425
  (2012), DOI 10.1111/j.1365-2966.2012.21439.x.
- Q. Zhu, L. Hernquist and Y. Li, *Numerical Convergence in Smoothed
  Particle Hydrodynamics*, ApJ 800:6 (2015),
  DOI 10.1088/0004-637X/800/1/6.
- D. Price and J. Monaghan, *An energy-conserving formalism for adaptive
  gravitational force softening*, MNRAS 374 (2007),
  DOI 10.1111/j.1365-2966.2006.11241.x.
