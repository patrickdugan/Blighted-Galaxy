# BG-PHY-24 | Radial homology energy and virial diagnosis

**Date:** 10 November 2026  
**Issue:** [#27](https://github.com/patrickdugan/Blighted-Galaxy/issues/27)  
**Depends on:** [BG-PHY-23](GRAVITY_COUPLED_HAMILTONIAN.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Why this experiment

BG-PHY-23 coupled the analytic full cold-electron EOS internal-energy derivative of the surface-gated affine MLS density with Poisson-matched variable-h gravity, and checked the complete Hamiltonian force against re-solved scalar energy finite differences.

The N256 matched kernel/density factorial then established that the boundary-gated density reduces undamped Rmax oscillation from 10.95% to 5.34% under M4 at eta=1.20, and from 15.10% to 5.43% under Wendland C4 at eta=1.70. It shifts the damped equilibrium farther inside the independent 8402.207 km continuum reference, however: M4 radius bias moves from -17.00% to -20.13%, and C4 from -17.26% to -22.49%.

BG-PHY-24 measures the energy and radial force that cause those displacements. The model is already an exact barotropic Hamiltonian, so a bounded **radial homology scan** can identify the discrete energy's preferred scale without introducing a force-tuning parameter.

## Homology and virial identity

For the mass centroid C and fixed initial physical displacements d_i = x_i-C, define a family of geometries:

```text
x_i(s) = C + s d_i,   s>0.

V(s) = U_Fermi[rho_hat(x(s))] + Phi_G[x(s)]

rho_hat = rho_SPH or the BG-PHY-20 surface-gated affine MLS density

Phi_G = G/2 sum_(i<j) m_i m_j
        [phi(r_ij,h_i(x)) + phi(r_ij,h_j(x))].
```

Every evaluation recomputes the exact safeguarded h/rho root, updated kernel volumes, affine MLS moment inverse and support partition. Neither the continuum radius nor any reference density enters these estimators.

The complete force from BG-PHY-23 satisfies F_i=-grad_i V. Consequently

```text
dV/ds = -sum_i d_i · F_i(s)

dV_P/ds = -sum_i d_i · F_pressure,i(s)
dV_G/ds = -sum_i d_i · F_gravity,i(s).
```

The use of **unscaled** d_i matters: x_i(s)-C=s d_i, so the virial sum computed with the physical positions at scale s is s times the homology derivative. At s=1 they coincide.

The analytic derivative is compared against the independent scalar evaluation:

```text
[dV/ds]_FD =
 [V(s+epsilon)-V(s-epsilon)]/(2 epsilon).

epsilon=1e-4.
```

The numerical discrepancy is normalized by
`abs(dV_P/ds)+abs(dV_G/ds)`, which remains well-defined even when the pressure and gravity terms nearly cancel at a stationary radius. The measured error must be below 1e-5 on every smooth active-set evaluation.

For each finite-difference pair, the experiment compares the **identity of every particle** on an h cap and records the MLS fallback fraction. If a cap or fallback branch changes, the point is flagged as a potentially nonsmooth derivative test rather than being counted as a passing central difference. Extrapolation across such a branch requires one-sided checks.

## Finding the preferred radial scale

Sample the bounded interval 0.75 ≤ s ≤ 1.30 in increments of 0.05. A candidate **minimum** requires a negative-to-positive sign change of the analytic derivative. The solver refines each bracket with fourteen safeguarded bisection iterations and evaluates the finite-difference energy curvature at the stationary point.

A positive-curvature stationary point is retained. If multiple stationary minima occur, choose the lowest-energy qualifying candidate. If no such bracket exists, the benchmark reports **unbracketed** and records both endpoint energies. It never replaces this with an invented equilibrium radius.

The implied radius bias is

```text
Rmax,stationary / Rref - 1
  = s_star * Rmax,starting / Rref - 1.
```

This is a one-parameter homology displacement of a given discrete particle geometry, not the unique relaxed solution of its full N-particle Hamiltonian. Angular rearrangements and nonhomologous motions may change the actual equilibrium.

## Eight-shell decomposition

At each frozen seed or post-damped configuration, sort particles by radius about the mass centroid, preserving equal mass quantile bins. For shell k, measure

```text
pressure_support_k =
  [sum_(i in k) m_i (a_P,i dot rhat_i)] /
  [sum_(i in k) m_i g_ref(r_i)]

gravity_inward_k =
  [-sum_(i in k) m_i (a_G,i dot rhat_i)] /
  [sum_(i in k) m_i g_ref(r_i)]

virial_P,k = sum_(i in k) (x_i-C) dot F_P,i

virial_G,k = sum_(i in k) (x_i-C) dot F_G,i

g_ref(r)=G M_1D(<r)/r².
```

The benchmark also reports RMS tangential pressure and net acceleration normalized by mass-weighted reference gravity, shell mean radius, original SPH density log-bias against the hydrostatic reference where available, mean kernel partition S, fractions on h caps, minimum Omega and MLS fallback fraction.

The **outer two shells**, containing 25% of the stellar mass, are summarized separately. The eight shell virial contributions must add back to the complete radial derivative at s=1 to machine precision.

## Controls and provenance

The factorial holds *all* of the following fixed:

- 0.65 solar mass, mu_e=2 cold Fermi EOS and independent Newtonian hydrostatic 1D reference;
- Matched M4 at eta=1.20 and matched Wendland C4 at eta=1.70, each tested with original and partition-gated affine MLS density;
- Initial shell and geometry-relaxed angular-glass particle positions at N256;
- Independent shell-orientation seeds 0 and 17;
- N256 3.00-second constant tau=0.7s damping as in BG-PHY-23 for the post-damping configuration, with no EOS tuning.

The **initial** stage compares all four Hamiltonian choices on identical frozen positions for each seed/layout. The **post** stage allows each Hamiltonian to create its own damped geometry and then repeats the energy scan and shell-force audit. An additional fixed-position cross-audit recomputes original and blended density forces on each post-damped geometry, distinguishing the direct density effect from shifts in the geometry.

The original stellar hold gate remains 3 dynamical times undamped, post-damped Rmax bias within 5%, Rmax and R99 excursion within 5%, mass-profile RMS within 5%, and established mass/energy/momentum assertions. A homology minimum alone does not satisfy those gates.

## Reproduction

```sh
cargo test --lib sph_homology -- --nocapture

# Reference-profile seeds: two layouts, two kernels, two densities.
cargo run --release --bin phy24_bench -- \
  --n 256 --seeds 0,17 --layout both --stage initial \
  --kernel both --density both \
  --s-lower 0.75 --s-upper 1.30 --s-step 0.05 --fd-eps 0.0001 \
  --csv bg-phy24-summary.csv \
  --profile-csv bg-phy24-homology.csv \
  --shell-csv bg-phy24-shells.csv

# Original damped preparation for the selected model; cross-evaluate
# both density methods at the same post-damping coordinates.
cargo run --release --bin phy24_bench -- \
  --n 256 --seeds 0 --layout glass --stage post \
  --kernel both --density both --damp-steps 300 \
  --dt 0.01 --damp-tau 0.7 \
  --csv bg-phy24-post-summary.csv \
  --profile-csv bg-phy24-post-homology.csv \
  --shell-csv bg-phy24-post-shells.csv \
  --cross-csv bg-phy24-post-fixed-position.csv
```

The [BG-PHY-24 workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy24.yml) publishes separate N256 initial and post-damping artifacts by shell-orientation seed. Regular CI includes an N64 smoke and Windows release build.

## Decision rule

First verify full analytic-versus-scalar radial energy agreement. Compare the initial particle homology minima to the independent reference radius, then test whether the post-damping configurations sit near their own homology energy minima. This distinguishes equilibrium contraction intrinsic to the discrete Hamiltonian from transient bulk oscillation.

Use the matched frozen-position controls to identify which part of a change in pressure support arises from the density functional and which comes from the radius shift. Identify the shell(s) with the greatest deficit of outward pressure relative to independently integrated gravity. The next conservative boundary intervention must improve the physical equilibrium radius while retaining the demonstrated dynamical stabilization and Hamiltonian conservation.

## Related

- [BG-PHY-23 gravity-coupled Hamiltonian results](GRAVITY_COUPLED_HAMILTONIAN.md).
- [BG-PHY-22 analytic surface-density variation](ANALYTIC_SURFACE_VARIATION.md).
- [BG-PHY-20 density and enthalpy reconstruction](FREE_SURFACE_DENSITY.md).
