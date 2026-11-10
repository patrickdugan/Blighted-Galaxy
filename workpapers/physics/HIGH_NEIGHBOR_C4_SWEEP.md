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

## Completed frozen-state measurements (10 November 2026)

The [N1024 and N2048 Actions matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42767907376)
succeeded and published separate summary and radial-shell CSV artifacts.
Every 100, 150, 200 neighbor target was reached with <0.1 mean-neighbor
absolute calibration error, at both resolutions, each geometry and each
kernel. No target case hits an upper smoothing-length cap or the Omega
floor. Exact Newton closures have <1e-9 dimensionless free-particle
h/rho equation residual, matching the prescribed tolerance.

**Wendland C4 on angular-glass seeds, exact Newton closure:**

| N | seed | target neighbors | eta | pressure vector RMS | gravity vector RMS | outer pressure/reference ratio | inner E0 RMS | outer E0 RMS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|1024|0|100|1.4503|0.3288|0.0199|0.8299|0.1021|0.5298|
|1024|0|150|1.6605|0.1662|0.0208|0.7781|0.0572|0.5723|
|1024|0|200|1.8410|0.1204|0.0232|0.7495|0.0354|0.6097|
|1024|17|100|1.4529|0.3451|0.0205|0.8348|0.1125|0.5375|
|1024|17|150|1.6611|0.1576|0.0211|0.7807|0.0568|0.5787|
|1024|17|200|1.8420|0.1133|0.0236|0.7516|0.0340|0.6158|
|2048|0|100|1.4483|0.5592|0.0156|0.9282|0.1195|0.4683|
|2048|0|150|1.6477|0.2398|0.0138|0.8559|0.0633|0.4797|
|2048|0|200|1.8154|0.1416|0.0146|0.8192|0.0394|0.5081|

Increasing target neighbors from 100 to 200 reduces C4 initial pressure
force RMS and inner E0 residual sharply in every tested arrangement.
The outer E0 residual grows at N1024, while the outer radial pressure
support decreases from ~0.83 to ~0.75. This is a clear smoothing/boundary
bias tradeoff and motivates tracking both indicators in dynamics.

The matched-neighbor **M4 glass control** still has lower frozen pressure
error at 200 neighbors: 0.1175 (N1024, seed 0) and 0.0903 (N2048, seed
0). M4's high-neighbor pairing/clumping stability is a separate physical
question, so these snapshots cannot establish superior dynamics for
either kernel.

At 200 neighbors, C4 glass pressure RMS changes 0.1204 (N1024) to
0.1416 (N2048). This finite-resolution comparison shows no monotone
decrease at fixed neighbor count. The neighbor-growth rule in BG-PHY-14
requires a larger suite of N with corresponding neighbor growth
for any convergence conclusion.

### Exact-root cross-check

At N1024, seed 0, glass layout, C4 target 200, the 24-pass closure has
maximum uncapped residual 1.78e-7; Newton closes to 4.14e-10. Both
pressure vector errors round to 0.120367. At N2048 the corresponding
max residuals are 1.99e-7 and 4.23e-10. This eliminates significant
underrelaxed-closure error as an explanation for remaining pressure bias.

### Selected physical holds

The [BG-PHY-17 full-hold workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy17-hold.yml)
runs three N1024 jobs, all using the same glass initializer and seed 0:

- Wendland C4 calibrated to 200 neighbors: best initial pressure-vector
  RMS among C4 target candidates at N1024; outer support ~0.75
- Wendland C4 calibrated to 150 neighbors: intermediate smoothing length
  with greater outer support (~0.78) and pressure RMS ~0.166
- Matched M4 control at eta=1.20: reference playable physics configuration
  with a safe lower realized neighbor population

Each physical hold keeps 3 seconds damping, 7.88 seconds undamped,
the reference radius and CDF gates, energy/momentum invariants,
and the exact recorded initial sampling offset. The new
`phy17_hold` runner emits a longitudinal 30-snapshot trace.
The workflow publishes output after each job completes.

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

## Completed full-duration holds (10 November 2026)

The [BG-PHY-17 hold workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42768741280) completed **all three N1024 seed-0 angular-glass simulations**, each with 3.00s tau=0.7s damping and 7.88s (3.005 t_dyn) undamped evolution. Every numerical job succeeded and emitted a summary CSV plus a 30-snapshot trajectory CSV.

| Measured quantity | M4, eta=1.20 | C4, target 150 | C4, target 200 |
|---|---:|---:|---:|
| Initial mean neighbors | 58.0020 | 149.9590 | 199.9297 |
| Post-damping mean neighbors | 56.6914 | 143.6973 | 190.3965 |
| Initial frozen pressure-vector RMS | 0.27472 | 0.16621 | 0.12037 |
| Post-damping Rmax (km) | 7252.654 | 7204.123 | 7170.033 |
| Post-damping radius bias vs 1D | −13.6816% | −14.2592% | −14.6649% |
| Contraction vs exact sampled seed | −10.2422% | −10.8428% | −11.2647% |
| Post-damping pressure-vector RMS | 0.07144 | 0.08131 | 0.08976 |
| Post-damping gravity-vector RMS | 0.06221 | 0.07606 | 0.08581 |
| Post-damping net-force RMS | 0.03492 | 0.02854 | 0.02451 |
| Worst undamped Rmax drift | 8.0972% | 10.0256% | 11.7982% |
| Worst undamped R99 drift | 3.4451% | 5.4304% | 7.1010% |
| Worst mass-quantile RMS | 0.021684 | 0.024423 | 0.027824 |
| Maximum relative energy drift | 1.285e−8 | 7.440e−9 | 4.853e−9 |
| Maximum relative momentum drift | 6.72e−17 | 5.30e−17 | 4.93e−17 |
| Full duration / physical gate | yes / **fail** | yes / **fail** | yes / **fail** |

All cases preserve mass, energy and momentum within the numerical assertions. The original physical acceptance remains a post-damping absolute radius error below 5%, Rmax and R99 excursion below 5%, and radial mass-CDF RMS below 5% for at least 3 dynamical times. The C4 200-neighbor case has the lowest *initial* pressure-vector error, but the worst physical-radius bias and undamped outer-radius excursion of these three.

### Interpretation and next hypothesis

With increasing C4 support, inner E0 noise and initial pressure-vector RMS decline, whereas the free-surface first-moment trace and pressure-gradient support fall. The larger smoothing length broadens the kernel beyond the free stellar surface. The comparison demonstrates a **finite-neighbor vs free-boundary gradient-bias tradeoff** in this discretization.

The near-zero post-damping net force expresses a rebalanced, contracted configuration rather than agreement with the independent reference star. The 200-neighbor C4 star has a smaller post-damping net-force RMS than M4, yet is more compact and less stable through the hold.

The next discriminating experiment should measure an explicitly **diagnostic** renormalized kernel-gradient and linear-consistency correction on frozen particles, especially the outer two mass shells. A corrected gradient evaluated independently of dynamics can expose the contribution of moment-matrix truncation versus sampled density error. Only a separately derived conservative force/energy formulation should be considered for dynamical adoption. See Bonet & Lok (1999), DOI 10.1016/S0045-7825(99)00051-1 and García-Senz, Cabezón & Escartín (2012), DOI 10.1051/0004-6361/201117939.

