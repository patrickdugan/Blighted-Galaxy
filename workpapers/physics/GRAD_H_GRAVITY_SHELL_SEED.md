# BG-PHY-09: Variational adaptive gravity and shell-balanced white dwarfs

**Date:** 10 November 2026  
**Status:** Research implementation with regression and headless A/B benchmarks.  
**Continuity:** [BG-PHY-08](ADAPTIVE_SPH_MOVING_OBSERVER.md) and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).  
**Affected modules:** `src/sph_adaptive.rs`, `src/sph_sampling.rs`, `src/bin/phy09_bench.rs`, and `src/bin/white_dwarf_lab.rs`.

The surviving Blighted Galaxy manuscripts use tidally disrupted compact objects and stellar-remnant weaponry. BG-PHY-09 continues reconstruction of the **numerical initial conditions** that these scenes will eventually need. This code was developed as a new prototype using published SPH methods; it is not recovered source from the lost engine.

## Mathematical references

- D. J. Price and J. J. Monaghan (2007), *An energy-conserving formalism for adaptive gravitational force softening in SPH and N-body codes*, MNRAS 374:1347–1358, [DOI 10.1111/j.1365-2966.2006.11241.x](https://doi.org/10.1111/j.1365-2966.2006.11241.x). This derives gravitational grad-h corrections when softening length depends on local SPH density.
- Pakmor, Edelmann, Röpke and Hillebrandt (2012), *Stellar GADGET*, MNRAS 424:2222–2231, [arXiv:1205.5806](https://arxiv.org/abs/1205.5806). Section 5 describes mapping a hydrostatic stellar density profile into regular equal-mass, angularly tiled shells, followed by relaxation.
- [BG-PHY-06](KERR_CAMERA_WHITE_DWARF.md) documents the Newtonian reference star and the fully relativistic zero-temperature electron EOS.

## 1. Differentiable adaptive Plummer potential

The BG-PHY-08 solver used an adaptive hydrodynamic smoothing length `h_i` with symmetric grad-h pressure, while **gravity retained a fixed Plummer softening** `epsilon`. This new module permits a separate optional gravitational model called `MatchedAdaptivePlummer`, whose pair potential and pair force depend on the same dynamically evolved `h_i`.

For two fluid parcels i,j, mass `m_i,m_j` and physical separation `r_ij`, define a negative unit potential

    phi(r,h) = -1/sqrt(r²+h²)

and pair gravitational energy

    U_ij = G m_i m_j [phi(r_ij,h_i)+phi(r_ij,h_j)] / 2.

At fixed h, the pair acceleration is attractive and symmetric, computed with the average of two Plummer gradients. To account for h depending on the local SPH density, also differentiate the energy through

    h_i = eta (m_i/rho_i)^(1/3),
    rho_i = sum_j m_j W(r_ij,h_i).

The correction to the density derivative uses

    Omega_i = 1 + h_i/(3rho_i) sum_j m_j dW_ij(h_i)/dh_i

and the gravitational adaptive-softening coefficient is

    zeta_i = -(h_i/(3rho_i)) sum_(j != i) m_j dphi(r_ij,h_i)/dh_i

with

    dphi/dh = h/(r²+h²)^(3/2).

Using the gradient of the cubic-spline hydrodynamic kernel `W` from BG-PHY-08, the extra acceleration on i is

    a_i,grad-h-grav =
       -(G/2) sum_j m_j [
          (zeta_i/Omega_i) grad_i W(r_ij,h_i)
         +(zeta_j/Omega_j) grad_i W(r_ij,h_j)
       ].

The implementation uses explicit equal-and-opposite pair updates, including the gravitational correction, and computes total potential energy from the same symmetric Plummer law.

**Important limits:** Plummer softening is a differentiable approximation to gravity, but its density distribution is **not the cubic-spline hydrodynamic kernel**. This is a variationally differentiated adaptive Plummer model, not a cubic-kernel Poisson-matched gravity solver. A finite h/rho fixed-point iteration, hard smoothing caps and a protective Omega floor can break ideal variational conservation. The central black hole retains a fixed Plummer field in this test, with no h-dependent correction from it. The gravitational model can still be used as an isolated self-gravitating star comparison.

The tests verify the analytic `dphi/dh` against finite differences, force symmetry, gravitational potential and energy-gradient consistency, short-duration energy drift, and finite state under force integration.

## 2. Equal-mass shell-balanced particle initialization

The original BG-PHY-06 initializer used antipodal Fibonacci directions, with radius determined by the cumulative mass fraction. It preserved mass and bulk momentum but correlated angular placement and radius, especially when only a few hundred fluid particles represented a star.

`seed_shell_balanced` keeps the precise radial **mass quantiles** from the original hydrostatic profile, while dividing the star into physical radial bands. Each band independently receives a rotated set of approximately uniform Fibonacci angular directions. Every generated sample has a counterpart at the opposite position, with identical mass and bulk velocity. Hence total mass, center of mass, and bulk momentum are preserved to floating-point precision. The deterministic seed also gives a reproducible sequence for physics comparisons.

This is a simple shell-stratified low-discrepancy approximation inspired by the logic of Pakmor et al. It does not implement their HEALPix tessellation or enforce exactly cubic equal-volume cells. A more faithful hydrostatic shell discretization would choose shell boundaries according to enclosed mass and enforce angular cell edge lengths comparable to shell widths, followed by multiple dynamical times of relaxation.

The code records the shell count, particle counts in each band, radial extrema and conserved initial bulk quantities.

## 3. Experimental comparisons

Hydrostatic reference from BG-PHY-06: cold carbon/oxygen `mu_e=2`, `M=0.65 M_sun`, `R=8,402.21 km`. The BG-PHY-09 benchmark tests three combinations using identical damping durations at each particle count:

- **Legacy/fixed:** previous antipodal seed, density-adaptive pressure h_i, fixed Plummer gravitational epsilon.
- **Shell/fixed:** shell-balanced seed, same adaptive pressure and fixed gravitational softening.
- **Shell/matched:** shell-balanced seed, adaptive pressure and symmetric h_i-dependent Plummer gravity with the gravitational zeta/Omega correction.

Each experiment logs radius, relative radius bias, RMS force ratio, logarithmic density mismatch and minimum Omega. The shell/matched cases at low N also undergo an undamped hold over three hydrostatic dynamical times.

### Results from the initial Linux CI benchmark

| N | Damping interval | Legacy/fixed radius | Shell/fixed radius | Shell/matched radius |
|---:|---:|---:|---:|---:|
| 64 | 3.0 s | 6,465 km | 6,687 km | 7,084 km |
| 128 | 3.0 s | 6,687 km | 6,941 km | 7,149 km |
| 256 | 1.2 s | 7,591 km | 7,635 km | 7,686 km |
| 512 | 1.2 s | 7,885 km | 7,868 km | 7,882 km |

At N=64–256, the `Shell/matched` model gives the closest radius for these damping intervals. Its remaining shortfall is **15.7%, 14.9%, and 8.5%**, respectively. At N=512, the previous legacy/fixed model has a marginally smaller radius error (−6.15%) than shell/matched (−6.20%), and shell/fixed is −6.36%.

The shell layout improves the logarithmic density mismatch at N=512. For `Legacy/fixed`, `Shell/fixed` and `Shell/matched`, the final log-RMS density errors are approximately **0.435, 0.248 and 0.246**. Their final-to-initial RMS force ratios are **0.215, 0.479 and 0.495**. These tests establish a measurable particle-placement tradeoff at short damping duration, while the longer undamped hold remains the decisive equilibrium gate.

Under undamped evolution for three dynamical times (7.87 physical seconds), the shell/matched star's radius expands by **28.3% at N=64**, **24.7% at N=128**, and **31.5% at N=256**. This is a more severe oscillation than the previous BG-PHY-08 setup. The experiments underscore that smoother or weaker gravity can improve a damped radius without producing a stable star.

**Conclusion for Issue #8:** Keep it open. The acceptance gate requires a radius error below 5% **and** an undamped evolution within 5% over at least three dynamical times. More particles and a closer match between hydrostatic density and SPH-resolved density are necessary, together with suitable timestep/convergence controls.

## 4. Numerical reproducibility

    cargo test --all-targets
    cargo run --release --bin phy09_bench -- --csv bg-phy09-comparison.csv
    cargo run --release --bin white_dwarf_lab

In the interactive stellar lab, **H** selects the older fixed-h relaxation, **J** selects adaptive pressure smoothing with fixed gravity, and **K** selects shell-balanced particles together with the new matched adaptive gravity. The post-relaxation dynamics use the corresponding gravitational model and their own energy diagnostic. The old methods remain available as references.

The new CI job publishes `bg-phy09-comparison.csv` with the complete A/B table and the undamped hold diagnostics. Tests and Windows compilation must pass before any changes reach `main`.

## 5. Further work

The next experiment should implement the actual Poisson-consistent SPH gravitational kernel paired with the hydrodynamic cubic spline and its matched grad-h correction, then audit finite-difference force/energy gradients at multiple smoothing lengths. The current Plummer model is a controlled intermediate.

The **shell sampling** should be extended to equal-volume tessellations with shell thickness tied to angular cell size, as in Pakmor et al. Large N and a proper neighbor index are required before claiming a converged white-dwarf SPH density profile. A gradual damping schedule and an isolated undamped hold over several sound-crossing times must be verified on the same particle realization.

No simulation of physically realistic tidal destruction should use these new stars as calibrated inputs until their free-evolution gates pass.
