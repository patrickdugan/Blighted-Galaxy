# BG-PHY-27: Radially Graded-Mass White-Dwarf Equilibrium

## 1. Motivation and Hypothesis

BG-PHY-24 through BG-PHY-26 established that the ~15%–17% contraction of equal-mass SPH cold white dwarfs (`0.65 M_sun`, carbon–oxygen `mu_e = 2`, `R_ref = 8,402.21 km`) is driven by the free-surface density/pressure-gradient deficit in the outermost shells, and that post-hoc free-boundary reconstructions (one-sided MLS ghost points in BG-PHY-25 and three-moment planar depth inversions in BG-PHY-26) suffer from spherical curvature bias and inter-shell variational coupling when inserted into a conservative Hamiltonian.

**BG-PHY-27** tests a purely Lagrangian, strictly conservative alternative that introduces **no** boundary ghost particles, **no** EOS modifications, and **no** confining potential: **radially graded particle masses** (`m_i` decreasing smoothly toward the stellar surface) under the existing exact unequal-mass Price–Monaghan variational SPH pressure + Poisson-matched self-gravity Hamiltonian (`sph_adaptive` + `sph_kernel`).

Because equal-mass sampling concentrates particle number density in proportion to `rho(r)`, the tenuous outer envelope of a degenerate white dwarf is severely under-resolved (`R_max` starts `6.68%` below `R_ref` at `N=256` and `5.06%` below at `N=512`). Grading particle number density as

$$n_\alpha(r) \propto \rho(r)^{1 - \alpha}, \qquad \alpha \in [0, 1],$$

places progressively finer radial and angular resolution in the low-density envelope (`alpha = 0` recovers equal mass; `alpha = 1` gives constant number density per unit volume and steep outer mass contrast).

---

## 2. Method (`src/sph_graded_seed.rs`)

1. **Graded radial CDF:** Using the independent 1D hydrostatic profile `white_dwarf::solve_hydrostatic_for_mass(0.65 M_sun)`, compute the cumulative number-density integral
   $$C_\alpha(r) = \frac{\int_0^r \rho(r')^{1-\alpha} r'^2\,dr'}{\int_0^{R_{\text{ref}}} \rho(r')^{1-\alpha} r'^2\,dr'}.$$
2. **Antipodal shell pairing with exact continuum enclosed mass:** Partition `N` particles (even `N`) into nearly cubic radial shells (`shell_count` chosen to match radial thickness to angular tile width, preserving exact `m_pair` symmetry across antipodal pairs `x_2k = -x_{2k+1}`). Each shell's enclosed mass `Delta M_k = M(r_{k+1}) - M(r_k)` is integrated directly from the independent 1D hydrostatic mass profile `M(r)` and divided equally among the `n_k` particles in shell `k`, guaranteeing:
   - Exact total mass `sum_i m_i = M_star` (`0.65 M_sun`),
   - Exact zero center-of-mass (`sum_i m_i x_i = 0`) and zero bulk momentum (`sum_i m_i v_i = 0`),
   - Smooth monotonic mass taper from core to surface.
3. **Unequal-mass antipodal angular glass (`relax_antipodal_angles_graded`):** Relaxes angular positions strictly on each shell's sphere (preserving exact `r_i`, `m_i`, and antipodal pairing) using local neighbor repulsion weighted by local shell spacing.
4. **Unmodified conservative dynamics:** Each seed is evolved with the exact unequal-mass variational SPH equations (`sph_adaptive` with Poisson-matched M4 at `eta = 1.20` and Poisson-matched Wendland C4 at `eta = 1.70`), verified independently in `sph_graded_seed::tests::unequal_mass_variational_forces_match_energy_gradient` against central finite differences of total energy `E = U_cold + Phi_grav` to `< 5e-5` relative error with linear momentum conserved to `< 1e-13`.
5. **Protocol and unchanged Issue #8 gates (`src/bin/phy27_bench.rs`):**
   - `3.00 s` constant damping (`300` steps at `dt = 0.01 s`), followed by
   - `7.88 s` (`3.005 t_dyn`, `788` steps at `dt = 0.01 s`) conservative undamped leapfrog hold.
   - **Issue #8 physical acceptance gates (strictly unchanged):**
     - **R (radius bias):** `|post_Rmax - R_ref| / R_ref <= 5%`
     - **H (undamped hold stability):** max excursion of both `R_max` and `R_99` over `>= 3 t_dyn` `<= 5%`
     - **C (mass-CDF accuracy):** max radial mass-CDF RMS against the 1D hydrostatic profile `<= 0.05`

---

## 3. Measured Results

Committed machine-readable datasets:
- [`data/bg-phy27/bg-phy27-n256-seed0.csv`](data/bg-phy27/bg-phy27-n256-seed0.csv)
- [`data/bg-phy27/bg-phy27-n256-seed17.csv`](data/bg-phy27/bg-phy27-n256-seed17.csv)
- [`data/bg-phy27/bg-phy27-n512-seed0.csv`](data/bg-phy27/bg-phy27-n512-seed0.csv)

### 3.1 N = 256 (`seed = 0` and `seed = 17`)

| N | Seed | Kernel | `alpha` | `m_max/m_min` | `initR` bias | Last-shell `P/g` | `postR` bias | Hold exc `Rmax / R99` | `maxCDF` | `dE/E` | Gate `(R, H, C)` |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 256 | 0 | M4 (`1.20`) | 0.00 | 1.0 | -6.68% | 0.731 | -17.00% | 11.00% / 9.53% | 0.0436 | `1.56e-8` | `F, F, T` |
| 256 | 0 | C4 (`1.70`) | 0.00 | 1.0 | -6.68% | 0.681 | -17.26% | 15.30% / 13.28% | 0.0528 | `1.22e-8` | `F, F, F` |
| 256 | 0 | M4 (`1.20`) | 0.25 | 4.4 | -4.26% | 0.735 | -15.06% | 10.47% / 7.88% | 0.0380 | `1.11e-8` | `F, F, T` |
| 256 | 0 | C4 (`1.70`) | 0.25 | 4.4 | -4.26% | 0.702 | -15.04% | 15.72% / 12.58% | 0.0446 | `7.99e-9` | `F, F, T` |
| 256 | 0 | M4 (`1.20`) | 0.50 | 31.9 | -2.28% | 0.758 | -13.76% | 8.44% / 7.10% | 0.0302 | `3.69e-8` | `F, F, T` |
| 256 | 0 | C4 (`1.70`) | 0.50 | 31.9 | -2.28% | 0.729 | -14.44% | 9.39% / 9.41% | 0.0324 | `2.43e-8` | `F, F, T` |
| 256 | 0 | M4 (`1.20`) | 0.75 | 545.7 | -0.90% | 0.779 | -9.44% | 36.93% / 5.71% | 0.0252 | `2.70e-7` | `F, F, T` |
| 256 | 0 | C4 (`1.70`) | 0.75 | 545.7 | -0.90% | 0.753 | -11.84% | 8.97% / 8.20% | 0.0315 | `4.68e-8` | `F, F, T` |
| 256 | 0 | M4 (`1.20`) | 1.00 | 47,144 | -0.19% | 0.817 | -5.23% | 147.88% / 8.24% | 0.0264 | `4.55e-6` | `F, F, T` |
| 256 | 0 | C4 (`1.70`) | 1.00 | 47,144 | -0.19% | 0.797 | -6.02% | 20.75% / 9.45% | 0.0320 | `2.72e-7` | `F, F, T` |
| 256 | 17 | M4 (`1.20`) | 0.00 | 1.0 | -6.68% | 0.707 | -17.59% | 10.73% / 8.78% | 0.0421 | `1.41e-8` | `F, F, T` |
| 256 | 17 | C4 (`1.70`) | 0.00 | 1.0 | -6.68% | 0.668 | -17.34% | 15.48% / 12.79% | 0.0495 | `1.20e-8` | `F, F, T` |
| 256 | 17 | M4 (`1.20`) | 0.25 | 4.4 | -4.26% | 0.740 | -14.43% | 11.19% / 7.80% | 0.0369 | `1.12e-8` | `F, F, T` |
| 256 | 17 | C4 (`1.70`) | 0.25 | 4.4 | -4.26% | 0.704 | -14.90% | 15.24% / 11.60% | 0.0432 | `6.85e-9` | `F, F, T` |
| 256 | 17 | M4 (`1.20`) | 0.50 | 31.9 | -2.28% | 0.761 | -11.25% | 10.03% / 5.84% | 0.0290 | `2.14e-8` | `F, F, T` |
| 256 | 17 | C4 (`1.70`) | 0.50 | 31.9 | -2.28% | 0.730 | -12.26% | 7.04% / 9.28% | 0.0311 | `6.02e-9` | `F, F, T` |
| 256 | 17 | M4 (`1.20`) | 0.75 | 545.7 | -0.90% | 0.788 | -13.34% | **4.04% / 4.32%** | 0.0295 | `2.16e-7` | `F, T, T` |
| 256 | 17 | C4 (`1.70`) | 0.75 | 545.7 | -0.90% | 0.758 | -13.44% | 7.00% / 7.15% | 0.0277 | `1.51e-7` | `F, F, T` |
| 256 | 17 | M4 (`1.20`) | 1.00 | 47,144 | -0.19% | 0.822 | -7.28% | 7.04% / 6.93% | 0.0314 | `5.46e-7` | `F, F, T` |
| 256 | 17 | C4 (`1.70`) | 1.00 | 47,144 | -0.19% | 0.801 | -7.58% | 33.02% / 9.71% | 0.0308 | `2.60e-7` | `F, F, T` |

### 3.2 N = 512 (`seed = 0`)

| N | Seed | Kernel | `alpha` | `m_max/m_min` | `initR` bias | Last-shell `P/g` | `postR` bias | Hold exc `Rmax / R99` | `maxCDF` | `dE/E` | Gate `(R, H, C)` |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 512 | 0 | M4 (`1.20`) | 0.00 | 1.0 | -5.06% | 0.765 | -15.47% | 9.35% / 5.62% | 0.0292 | `1.18e-8` | `F, F, T` |
| 512 | 0 | C4 (`1.70`) | 0.00 | 1.0 | -5.06% | 0.720 | -16.04% | 12.40% / 8.81% | 0.0354 | `6.78e-9` | `F, F, T` |
| 512 | 0 | M4 (`1.20`) | 0.50 | 43.8 | -1.53% | 0.823 | -11.81% | 7.11% / 4.02% | 0.0204 | `2.21e-8` | `F, F, T` |
| 512 | 0 | C4 (`1.70`) | 0.50 | 43.8 | -1.53% | 0.780 | -13.95% | 6.00% / 6.30% | 0.0221 | `1.16e-8` | `F, F, T` |
| 512 | 0 | M4 (`1.20`) | 0.75 | 995.0 | -0.54% | 0.842 | **-9.82%** | **4.82% / 4.26%** | **0.0178** | `1.25e-7` | `F, T, T` |
| 512 | 0 | C4 (`1.70`) | 0.75 | 995.0 | -0.54% | 0.803 | -10.21% | 5.94% / 5.44% | 0.0192 | `1.02e-7` | `F, F, T` |

---

## 4. Physical Interpretation and Verdict

1. **Direct confirmation of the outer-resolution bottleneck:** Grading particle mass toward the surface (`alpha = 0 -> 0.75`) raises initial last-shell pressure support `P/g` from `0.765` to `0.842` at `N=512`, shrinks the initial sampling radius shortfall by ~10x (`-5.06% -> -0.54%`), cuts post-damping radial contraction by ~37% (`-15.47% -> -9.82%`), and cuts radial mass-CDF RMS by ~39% (`0.0292 -> 0.0178`), all with `0%` smoothing-length cap hits, `dE/E <= 1.25e-7`, and `|dP| ~ 1e-16`.
2. **First 3-dynamical-time undamped hold excursion pass (`<= 5%`):** At `N=512, seed=0, alpha=0.75, M4`, both `R_max` (`4.82%`) and `R_99` (`4.26%`) stay within the `5%` excursion gate over `3.005 t_dyn` (`N=256, seed=17, alpha=0.75, M4` also passes at `4.04% / 4.32%`).
3. **High mass-ratio surface splash (`alpha >= 0.75` at low `N`, `alpha = 1.00` generally):** When `m_max/m_min` exceeds `~500` at coarse angular resolution (`N=256`), an isolated light outer particle can suffer seed-dependent ejection (`36.93%` at `N=256, seed=0, alpha=0.75` M4; `147.88%` at `alpha=1.00` M4; `33.02%` at `seed=17, alpha=1.00` C4) even while `R_99` remains tight (`5.7%–9.7%`). Moderate grading (`alpha = 0.50`, `m_max/m_min ~ 32–44`) is free of splash across all seeds and kernels, while `N=512` stabilizes `alpha = 0.75` under both M4 and C4.
4. **Issue #8 Verdict (HONEST):** **Issue #8 remains OPEN.** Although `N=512, alpha=0.75` M4 passes both the **H** (`<= 5%` excursion over `3.005 t_dyn`) and **C** (`<= 0.05` mass-CDF RMS) gates, its post-damping outer radius (`7,577 km`, `-9.82%` vs `8,402 km`) still exceeds the `<= 5%` **R** radius-bias gate. The dedicated workflow `.github/workflows/bg-phy27.yml` reproduces the full `N=256` (`seed=0,17`) and `N=512` (`seed=0`) matrix in CI, while higher-resolution (`N >= 1024`) graded-mass runs remain the next step to close the remaining radius gap without altering playable game physics.

## 5. Reproduction

~~~sh
cargo test --lib sph_graded_seed
cargo run --release --bin phy27_bench -- --n 256 --seeds 0,17 --alphas 0.0,0.25,0.50,0.75,1.0 --kernel both --relax-steps 300 --hold-steps 788 --dt 0.01 --csv bg-phy27-n256.csv
cargo run --release --bin phy27_bench -- --n 512 --seeds 0 --alphas 0.0,0.50,0.75 --kernel both --relax-steps 300 --hold-steps 788 --dt 0.01 --csv bg-phy27-n512.csv
~~~
