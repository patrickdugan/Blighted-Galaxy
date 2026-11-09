# BG-PHY-10: Poisson-matched M4 gravity and approximately cubic stellar shells

**Status:** Experimental clean-room Rust reconstruction, 10 November 2026.  
**Predecessor:** [BG-PHY-09](GRAD_H_GRAVITY_SHELL_SEED.md).  
**Convergence tracker:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).  
**Code:** `src/sph_gravity_kernel.rs`, `src/sph_adaptive.rs`, `src/sph_sampling.rs`, `src/bin/phy10_bench.rs`, `src/bin/white_dwarf_lab.rs`.

## 1. Sources and design goal

- Price & Monaghan (2007), *An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes*, MNRAS 374, 1347–1358, [DOI](https://doi.org/10.1111/j.1365-2966.2006.11241.x), especially Appendix A (equations A1–A3). Their equations provide the potential, force and smoothing-length derivative of the standard cubic SPH mass kernel.
- Pakmor et al. (2012), *Stellar GADGET*, MNRAS 424, 2222–2231, [Section 5](https://academic.oup.com/mnras/article/424/3/2222/977734): approximate equal-volume stellar particle placement in radial shells, with angular tiles and a subsequent isolated relaxation.
- [Phantom SPH method review](https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/phantom-a-smoothed-particle-hydrodynamics-and-magnetohydrodynamics-code-for-astrophysics/182F569916933E20686B8032111A6182), Section 2.12.2: adaptive gravitational softening, symmetrization, grad-h correction.

The physical target remains an ideal Newtonian hydrostatic cold white dwarf of **0.65 solar masses**, radius **8,402.21 km**, with a fully relativistic zero-temperature electron EOS at mean electron molecular weight mu_e=2. The stellar-fluid samples must converge to this independently calculated continuum reference *and* remain stable when damping is disabled. A smaller numerical force residual or a closer short-duration radius alone cannot establish hydrostatic equilibrium.

The original fiction and workpapers are preserved in the archival folders. The SPH algorithms in this workpaper are newly reconstructed scientific prototypes.

## 2. M4 density kernel and the matching gravitational potential

The three-dimensional cubic-spline density kernel with smoothing length h and q=r/h is

    W(r,h) = 1/(pi*h^3) * f(q)
    f(q) = 1 - (3/2)q² + (3/4)q³          (0 <= q < 1)
         = (1/4)(2-q)³                    (1 <= q < 2)
         = 0                               (q >= 2).

BG-PHY-09 used a differentiable Plummer potential with h_i-dependent gravitational softening. BG-PHY-10 constructs the gravitational pair potential from **Poisson's equation for the same cubic-spline mass kernel**, giving Newtonian gravity exactly outside its compact support.

For the unit-mass potential phi(r,h), with sign convention phi = -1/r at r >= 2h:

    phi(r,h) = (1/h) * F(q)

    F(q) = (2/3)q² - (3/10)q⁴ + (1/10)q⁵ - 7/5    (0 <= q < 1)
         = (4/3)q² - q³ + (3/10)q⁴ - (1/30)q⁵
           - 8/5 + 1/(15q)                       (1 <= q < 2)
         = -1/q                                   (q >= 2).

The nonnegative radial derivative of the potential is

    dphi/dr = (1/h²) * F_force(q)

    F_force(q) = (4/3)q - (6/5)q³ + (1/2)q⁴           (0 <= q < 1)
               = (8/3)q - 3q² + (6/5)q³
                 - (1/6)q⁴ - 1/(15q²)              (1 <= q < 2)
               = 1/q²                              (q >= 2).

The partial derivative with respect to h at fixed physical r is

    dphi/dh = (1/h²) * F_h(q)

    F_h(q) = -2q² + (3/2)q⁴ - (3/5)q⁵ + 7/5         (0 <= q < 1)
           = -4q² + 4q³ - (3/2)q⁴
              + (1/5)q⁵ + 8/5                      (1 <= q < 2)
           = 0                                      (q >= 2).

These functions are implemented independently in `sph_gravity_kernel.rs`. Tests check potential/force continuity near q=1,2, analytic radial and h derivatives against central finite differences, and the Poisson enclosed-mass identity

    r² * dphi/dr = 4pi * integral_0^r [W(s,h) s² ds].

For r>=2h the right-hand side is exactly 1, and the gravitational force returns to Newton's inverse-square law. This property is absent from Plummer softening, which remains modified at every finite separation.

## 3. Symmetric variable-h gravitational force

Each SPH parcel has mass m_i and smoothing length h_i, defined by

    h_i = eta*(m_i/rho_i)^(1/3)
    rho_i = sum_j m_j W(r_ij,h_i).

The symmetric pair potential energy is

    U_ij = (G m_i m_j/2) [phi(r_ij,h_i)+phi(r_ij,h_j)].

The attractive inverse-square-style contribution to acceleration is

    a_i,pair = G m_j (r_j-r_i)/r_ij *
               0.5 [phi_r(r_ij,h_i)+phi_r(r_ij,h_j)].

Because h varies with density, differentiating the *same* energy also adds the gravitational grad-h term, with

    Omega_i = 1 + (h_i/(3rho_i)) sum_j m_j dW(r_ij,h_i)/dh_i
    zeta_i  = -(h_i/(3rho_i)) sum_(j != i) m_j dphi(r_ij,h_i)/dh_i.

The additional acceleration is

    a_i,grad-h = -(G/2) sum_j m_j [
        (zeta_i/Omega_i) grad_i W(r_ij,h_i) +
        (zeta_j/Omega_j) grad_i W(r_ij,h_j)].

Total pair accelerations use equal-and-opposite updates. The gravitational energy diagnostic now uses the matching kernel potential for each pair, so changes in h are included in the same modeled energy expression.

**Practical caveats:** h and rho are solved through a finite, damped fixed-point iteration. The h bounds and Omega floor can break exact Lagrangian consistency. The black-hole central body retains a fixed Plummer source field, which is isolated from the new stellar-self-gravity experiments. The code uses global O(N²) gravitational pairs and does not yet implement an adaptive neighbor tree or multipole approximation. Even if the conservation tests pass, a variable-h kernel needs numerical gradient and time-step convergence auditing.

The selectable gravity choices are `FixedPlummer`, `MatchedAdaptivePlummer` (BG-PHY-09) and `MatchedCubicSpline` (BG-PHY-10). The old models remain available for A/B testing.

## 4. Approximately cubical stellar shell initialization

Pakmor et al. divide a hydrostatic stellar profile into mass-equivalent angular cells and radial layers whose thickness is comparable to an angular cell edge. BG-PHY-10 implements a deterministic approximation, keeping the exact mass-CDF radial quantiles and antipodal pairs.

For a candidate shell containing k *pairs* (2k equal-mass particles), let its radial mass-CDF boundaries be r_lower and r_upper, and let r_mid be its mass-median radius. The approximate equal-area angular tile edge is

    ell_tan ≈ r_mid sqrt(4pi/(2k)) = r_mid sqrt(2pi/k).

The desired shell thickness is `Delta r = r_upper-r_lower ≈ ell_tan`. A greedy search chooses k by minimizing the absolute logarithm of the shell-thickness/tile-edge ratio, with a small penalty for tiny shells. Each shell receives a separately rotated Fibonacci arrangement of k angular directions and their antipodes. Every equal-mass pair inherits its own radial CDF quantile. Therefore the initial star conserves total mass, center of mass and bulk velocity, while reducing the tight radius/angle correlation of earlier low-N samplers.

**This is not a full HEALPix implementation or an exact equal-volume Voronoi tessellation.** Finite k, the greedy shell split and the antipodal construction introduce geometric noise. The new initializer reports shell count, aspect ratios and conservation error, and the benchmark tests whether those choices improve the *physical hydrostatic density* after SPH relaxation.

## 5. Numerical benchmark contract

Run:

    cargo test --all-targets
    cargo run --release --bin phy10_bench -- --csv bg-phy10-comparison.csv
    cargo run --release --bin white_dwarf_lab

At each N in 64, 128, 256, 512 and 1024, the benchmark uses the same relaxation duration for every available variant at that N. The variants include legacy sampling, the previous radially grouped shells, the approximately cubic shells with adaptive Plummer gravity, and the approximately cubic shells with M4 Poisson-matched gravity. For computing cost, the absolute damping duration depends on N; a complete cross-resolution convergence claim will require repeated runs with common physical duration and controlled timestep.

The output CSV records particle count, initializer, gravity law, damping interval, radial shell shape, residual-force ratio, logarithmic density-profile mismatch, relaxed radius and its bias, minimum Omega, smoothing cap incidence, and timing. The Poisson gravity model also receives undamped three-dynamical-time holds at selected particle counts. CI publishes the resulting CSV without suppressing failing physics assertions.

### Acceptance criteria from Issue #8

- Relaxed radius error <5% from the independent hydrostatic reference.
- Residual density-profile error decreases under increased resolution.
- Isolated undamped evolution over at least three stellar dynamical times holds radius within 5% of its initial relaxed value and reports energy/momentum drift.
- Tests reproduce the physical kernel's Poisson equation, analytic force derivative, and total-energy-gradient consistency.

## 6. First CI results and physical interpretation

The first complete Linux CI benchmark passes the M4 potential/force derivative, Poisson enclosed mass, energy-gradient and symmetric-momentum checks. It runs the reference star at N=64,128,256,512,1024, using 3.0 s damping at N<=128, 1.2 s at N=256 and 512, and 0.6 s at N=1024. The stellar dynamical time for the continuum reference is **2.62223 s**, so the largest-N cases are only short-time probes.

Here, **legacy** means the original antipodal Fibonacci initializer; **near-cubic** means the new radial/angle cell-aspect initializer. The first two gravity variants use variable Plummer softening, and the final gravity variant uses the new Poisson-matched M4 law.

| N | Damping | Legacy + variable Plummer radius | Near-cubic + variable Plummer radius | Near-cubic + cubic-M4 radius | Cubic-M4 density log-RMS |
|---:|---:|---:|---:|---:|---:|
| 64 | 3.0 s | 6,830 km | 7,076 km | 6,705 km | 0.3033 |
| 128 | 3.0 s | 6,923 km | 7,186 km | 6,885 km | 0.2584 |
| 256 | 1.2 s | 7,620 km | 7,595 km | 7,549 km | 0.2216 |
| 512 | 1.2 s | 7,879 km | 7,756 km | 7,735 km | 0.1826 |
| 1024 | 0.6 s | 8,114 km | 8,001 km | 7,999 km | 0.1667 |

The cubic-M4 model produces a smaller logarithmic density error in these runs, especially at N=1024 where its 0.1667 compares with 0.5670 for legacy+Plummer. That represents about **71% lower log-RMS discrepancy** for the short-time density sample. However, at the same N and duration the legacy+Plummer radius is closer to the hydrostatic reference (−3.43% versus cubic-M4's −4.80%), and the cubic-M4 force-residual ratio is still **0.528**. The 1024-particle star was damped for only **0.23 dynamical times**. Its radius matching the 5% tolerance at this moment does not mean its particles are in force-balanced hydrostatic equilibrium.

The new M4 kernel is distinctly compact: it recovers the exact Newtonian pair force outside two smoothing lengths. Plummer softening remains gravity-modified at arbitrary separation. Its potential consistency comes from the Poisson integral, while its physical usefulness must be checked against the live star's resolved density and motion.

### Undamped stability probes

After a three-dynamical-time free evolution (7.87 seconds) with the cubic-M4 / near-cubic configuration:

| N | Radius change | Relative total-energy drift | Relative momentum drift |
|---:|---:|---:|---:|
| 128 | −10.97% | 3.15×10^−8 | 1.14×10^−16 |
| 256 | −14.39% | 3.87×10^−8 | 4.33×10^−17 |

These stars contract substantially despite excellent integration energy and momentum diagnostics. The samples are not yet stable fluid equilibria. **Issue #8 stays open**. To pass it, test at comparable physical damping durations, refine particle counts and local smoothing behavior, and verify a post-relaxation undamped radius variation within 5% over three dynamical times.

CI outputs `bg-phy10-comparison.csv` under the `bg-phy10-reference` artifact. The interactive laboratory exposes **L** for the new nearly cubic particle layout and M4 gravitational kernel, and keeps the earlier H/J/K modes for reference.

## 7. Scientific status and next steps

The new kernel law is mathematically consistent with the cubic SPH density kernel. The remaining astrophysical credibility requires empirical convergence. Early particle counts have an unrealistically small radius after damping and some variants expand markedly during an undamped hold.

After evaluating this milestone, prioritize a real neighbor index for pressure/density updates and a suitable tree/FMM approximation for long-range gravity, so N>=4096 becomes practical. Use matched physical durations across resolutions, adaptive error-controlled timesteps, and more faithful equal-mass shell tessellation. A stable isolated star should then be inserted into a white dwarf/black-hole tidal encounter with GR capture and appropriate equation-of-state limitations kept explicit.
