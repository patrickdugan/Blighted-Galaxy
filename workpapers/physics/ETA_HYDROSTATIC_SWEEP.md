# BG-PHY-13: Smoothing-length sensitivity and hydrostatic gravity convergence

**Date:** 10 November 2026  
**Status:** Controlled parameter-sweep implementation. Physical acceptance pending measured matched-duration results.  
**Issue:** [#8, cold white-dwarf equilibrium convergence](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).  
**Baseline:** [BG-PHY-12](RADIAL_FORCE_RELAX_SCHEDULE.md).  
**Files:** `src/bin/phy13_bench.rs`, `.github/workflows/bg-phy13.yml`, README and workpaper index.

The reconstructed *Blighted Galaxy* stellar physics now uses a cold-electron equation of state, a Poisson-matched cubic-M4 gravitational kernel, adaptive hydrodynamic smoothing, and a shell-stratified equal-mass particle initializer. BG-PHY-12 showed that the SPH star can reach a near-unit global pressure-to-gravity virial ratio while its radius remains 17–21% below the independent hydrostatic reference. This experiment targets a specific remaining uncertainty: the effective force resolution set by **the adaptive smoothing-length parameter**.

## 1. Source-grounded parameter

Price and Monaghan (2007), *An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes*, [MNRAS 374:1347](https://academic.oup.com/mnras/article/374/4/1347/963634), set the smoothing length as

    h_i = eta * (m_i / rho_i)^(1/3).

Here eta is dimensionless and controls the number of neighbors inside the kernel's compact support at radius 2h. In a spatially uniform 3D equal-mass fluid, the expected count is approximately

    N_support ~ (4*pi/3) * (2*h)^3 * (rho/m)
              = (32*pi/3) * eta^3.

For example, eta=1.2 corresponds to approximately 58 interacting particles in a uniform distribution, near the paper's stated ~60-neighbor reference. The real star is stratified, and its discretized density estimate includes the self-contribution. Therefore the benchmark reports both **actual kernel-support neighbors** and total candidate counts from the spatial cell index. Candidate counts can be much larger and are not a physically equivalent neighbor-number measure.

Changing eta changes the adopted smoothing lengths and thereby the pressure and softened gravitational force *resolution*. The cold electron EOS, total stellar mass, hydrostatic target, particle seed, gravitational law and step size remain fixed.

The 2007 work also notes that iterative h/rho convergence affects accuracy. Our existing finite 24-pass underrelaxed density-h solver, hard h caps and Omega floor are additional approximation sources. They will be monitored rather than silently assumed converged.

Pakmor et al. (2012), *Stellar GADGET*, [MNRAS 424:2222](https://academic.oup.com/mnras/article/424/3/2222/977734), demonstrate the sensitivity of white-dwarf particle stars to the initial radial distribution and hydrodynamic relaxation. The authors' stable models use much higher particle counts and more sophisticated shell tessellation. Our low-N sweep is a numerical diagnostic, with no claim of comparable astrophysical resolution.

## 2. Controlled experiment design

The independent reference remains a 0.65-solar-mass, mu_e=2, zero-temperature white dwarf with R_ref=8,402.207 km, as integrated in the BG-PHY-06 Newtonian hydrostatic solver. The dynamical time is

    t_dyn = sqrt(R_ref^3/(G*M)) ~= 2.62223 s.

For each resolution, the runner creates **one canonical** nearly cubic equal-mass shell seed and clones it across eta variants. Each simulation uses the identical cold electron EOS, cubic-M4 Poisson gravity, physical timestep, and **constant damping tau=0.7 s**. This retains the more stable schedule selected in BG-PHY-12 as the reference while varying only eta.

The high-resolution workflow compares these planned settings:

| Particle count | eta values | Relaxation | Undamped hold |
|---:|---|---:|---:|
| 256 | 0.85, 1.05, 1.20, 1.45 | 3.00 s | 7.88 s |
| 512 | 0.85, 1.20, 1.45 | 3.00 s | 7.88 s |
| 1024 | 0.85, 1.20, 1.45 | 3.00 s | 7.88 s |

The kernel-support estimate in a uniform sample is:

| eta | Expected neighbors |
|---:|---:|
| 0.85 | 20.6 |
| 1.05 | 38.8 |
| 1.20 | 57.9 |
| 1.45 | 102.1 |

The benchmarks measure the actual neighborhood count, which can differ from this approximation, especially near the stellar surface.

A small N=64, 0.2-s relaxation and 0.3-s undamped run is included in ordinary Rust CI solely as a compilation/finite-numerics smoke test. Its physical acceptance flag must remain false because it covers less than three dynamical times.

## 3. Radial hydrostatic force diagnostics

After damping, the benchmark calls the BG-PHY-12 force decomposition. The spherical reference field is the independently integrated enclosed mass:

    g_ref(r) = G * M_hydro_enclosed(r) / r^2.

The SPH inward gravitational acceleration and outward pressure acceleration are mass-weighted in eight radial equal-enclosed-mass bins. The new scalar gravity-reference discrepancy is

    E_g = sqrt[
       sum_bins w_b (g_sph,b - g_ref,b)^2 /
       sum_bins w_b g_ref,b^2
    ],

where w_b is the bin's mass fraction. This metric distinguishes a pressure/gravity equilibrium *inside the numerical model* from agreement with the stellar hydrostatic structure. The benchmark also retains the BG-PHY-12 pressure-to-gravity force-virial and RMS net/gravity acceleration ratio.

Additional diagnostics include initial shell field ratios at the innermost and outermost mass fractions, actual neighbor support counts, Omega minimum, smoothing cap fraction, post-damping Rmax and R99, normalized mass-CDF RMS, and full energy/momentum invariants.

## 4. Undamped acceptance gate

Immediately after damping, each variant is evolved under the same conservative SPH dynamics, with no artificial drag, for **at least three hydrostatic dynamical times**.

Each hold is sampled throughout rather than only at its end. The runner stores the greatest absolute fractional change in outermost radius and in 99%-enclosed-mass radius, and the largest normalized radial mass-CDF discrepancy.

Physical acceptance requires all of:

- after-damping outer radius within 5% of independent R_ref;
- maximum absolute Rmax and R99 excursions below 5% during a hold of at least 3*t_dyn;
- mass-CDF RMS below 5% (additional scientific criterion);
- numerical mass and momentum consistency, with bounded relative energy drift.

The benchmark asserts the numerical invariants, while the *physical* criteria are reported as explicit true/false fields. A simulator passing CI does not mean the white dwarf has passed physical equilibrium.

It also reports an **exploratory within-run score**, the sum of the absolute percentage radius error, maximum Rmax and R99 drift percentages, and 100 times the maximum normalized CDF mismatch. The minimum-score eta is a hypothesis for **independent** confirmation at another resolution and longer hold, not an accepted tuned parameter.

## 5. Run and artifact contract

Quick smoke:

    cargo run --release --bin phy13_bench --       --n 64 --etas 0.85,1.20 --relax-steps 20 --hold-steps 30       --csv bg-phy13-smoke.csv --shell-csv bg-phy13-smoke-shells.csv

Long fixed-duration validation:

    cargo run --release --bin phy13_bench --       --n 512 --etas 0.85,1.20,1.45 --relax-steps 300 --hold-steps 788       --csv bg-phy13-512.csv --shell-csv bg-phy13-shells-512.csv

The dedicated GitHub Actions **BG-PHY-13 eta and resolution sweep** workflow runs at three resolutions in parallel, publishing one summary CSV and one radial shell CSV per resolution. The Windows release build includes the headless `phy13_bench.exe` executable.

The summary CSV contains all eta settings and measured physical gates. The shell CSV exposes eight mass-fraction pressure/gravity comparisons per eta after damping and at the end of the undamped hold.

## 6. Decision and open questions

The BG-PHY-12 constant damping baseline is retained in the playable physics modes. The sweep changes an experimental runner parameter and does not retune the cold electron EOS. Any eta that improves the radius while damaging force balance, density CDF, or undamped stability will remain a research result.

If radius and field errors change strongly with eta, the next physical tasks are a more accurate h/rho root solver (the Price–Monaghan per-particle Newton procedure), neighbor-count convergence in outer layers, and higher-fidelity shell tilings. If eta has little effect on the residual bias, structural initialization or missing physical corrections deserve priority. A validated stellar explosion or high-fidelity white-dwarf shredding scene remains downstream of [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).
