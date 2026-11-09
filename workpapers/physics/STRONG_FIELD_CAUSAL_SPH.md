# BG-PHY-04: Geodesics, causal records, and self-gravitating breakup

**Date:** 9 November 2026  
**Status:** New Rust modules with executable tests and an interactive lab.  
**Previous stage:** [BG-PHY-03](SI_GAMEPLAY_BRIDGE.md).  
**Primary source design:** [Game Design Treatment](../design/game-design-treatment.txt), [*Patricia Tries*](../../stories/manuscripts/patricia-tries.txt) and [*White Dwarf Super Shotgun*](../../stories/manuscripts/white-dwarf-super-shotgun.txt).

This is the next physics milestone in the reconstruction. The original lost Rust engine remains unrecovered; all code in this milestone is new. Physical calculations are distinguished from speculative setting technology.

## 1. Exact test-particle Schwarzschild geodesics

**Reference:** Oregon State University, *Geometry of General Relativity*, [Geodesic Equation](https://sites.science.oregonstate.edu/physics/coursewikis/GGR/book/ggr/schwgeo.html) and [Circular Orbits](https://sites.science.oregonstate.edu/physics/coursewikis/GGR/book/ggr/ocircle.html).

The equation of motion describes a freely falling negligible-mass particle in the equatorial plane around an isolated, **nonspinning**, uncharged black hole. It is exact in this background spacetime outside the event horizon. It does not include the test particle's gravitational backreaction, spin coupling, electromagnetic propulsion, Kerr frame dragging or multi-black-hole interactions.

Set gravitational radius `r_g=GM/c²`, dimensionless areal radius `x=r/r_g`, dimensionless proper time `s=cτ/r_g`, specific angular momentum `ℓ=L/(c r_g)`, and conserved specific energy `ε=E/(mc²)`. The equations integrated by RK4 are:

```text
dx/ds       = u
du/ds       = -1/x² + ℓ²/x³ - 3ℓ²/x⁴
dφ/ds       = ℓ/x²
d(ct/r_g)/ds= ε / (1 - 2/x)

Constraint: u² + (1 - 2/x)(1 + ℓ²/x²) = ε²
```

For circular timelike geodesics, `x>3` and:

```text
ℓ = x / sqrt(x-3)
ε = (1-2/x) / sqrt(1-3/x)
Ω_coordinate = sqrt(GM/r³)
```

The **innermost stable circular orbit (ISCO)** is at `6 r_g`. The photon sphere is at `3 r_g`; photon propagation is not integrated by the present timelike solver. For a slightly eccentric, weak-field orbit with semilatus rectum `p` in `r_g`, the leading periapsis advance is `Δφ≈6π/p`. The executable checks a resolved relativistic periapsis return against this order-of-magnitude prediction.

**Numerical method:** Fixed capped steps of `Δs≤0.04` with a fourth-order Runge–Kutta integrator and measured effective-energy residual. The integration stops at `x=2.01` because Schwarzschild coordinate time diverges at the horizon. The resulting `HorizonCoordinateCutoff` marks a chart boundary, not physical destruction. A future horizon-regular chart (Kerr–Schild or Eddington–Finkelstein) is required for smooth crossing.

**Test gates:** Schwarzschild photon/ISCO radii, circular orbital constants, one-orbit radial closure and constraint error, eccentric relativistic advance, and radial infall to a finite coordinate cutoff.

## 2. Historical light cones

**Module:** `src/causal.rs`.  
**Connected gameplay:** `src/si_bridge.rs`, the default relativistic joust.

Every observed actor has a bounded, time-ordered record of `(t, position, velocity)`. Samples must be finite, strictly increasing in coordinate time, with individually subluminal velocities **and subluminal displacement between successive samples**. There is no teleporting across a step.

A recorded trajectory is piecewise-linear over each known interval. The sensor receives a sample when:

```text
|source_position(t_emit) - observer_position(t_observe)|
    = c × (t_observe - t_emit)
```

For subluminal trajectories this equation has at most one causal root; bisection finds it inside the recorded interval. If the needed emission predates the oldest sample or lies beyond the newest available sample, the sensor returns an **unavailable** state. It does not fabricate an earlier event by extrapolating current position or velocity.

**Integration:** The main encounter records the player and enemy at 120 Hz in histories capped at 8,192 samples per participant (approximately 68 seconds before eviction). Visual ghost markers and hostile firing directions now come from recorded causal histories. The enemy can only take a targeted shot after historical observations have reached it; when its view is unavailable it defers firing. Worldline position is linearly interpolated; a dense causal event queue or polynomial ephemeris is deferred.

**Test gates:** Stationary 1-light-second delay, analytic receding source at 0.5c, accelerating piecewise trajectory, missing past and future, sample capacity eviction, strictly ordered time, rejection of teleports even when endpoint velocities are causal, and bounded retarded shots.

**Limitations:** Flat-spacetime light propagation is used for the encounter. Gravitational null geodesics, Shapiro delay, curved-spacetime ray tracing, Doppler spectra and lensing are separate tasks. The visual game displays an optical marker representing recorded events; it does not render physically computed lensing.

## 3. Pairwise self-gravitating tidal debris

**Module:** `src/debris.rs`.  
**Connected gameplay:** `SiTidalEncounter` after the tidal diagnostic exceeds 1.

A white-dwarf breakup initializes 96 symmetric tracer parcels whose summed mass, center of mass and bulk momentum match the progenitor to numerical precision. Each particle thereafter feels gravity from the central black hole, gravity from **all other debris particles**, and an approximate fluid-pressure force. Every gravitational and pressure interaction is implemented as an equal-and-opposite pair. The black hole responds to the gravitational recoil of the debris.

This is **2D smoothed-particle hydrodynamics (SPH-like)** rather than a 3D degenerate-star solver. The physical state of a parcel uses SI metres, seconds and kilograms. Surface density `σ` is measured in kg/m² and approximated through the compact-support cubic-spline kernel of smoothing length `h`:

```text
q = r/h
W_2D(q,h) = 10/(7πh²) ×
  1 - 1.5q² + 0.75q³,  0 ≤ q < 1
  0.25(2-q)³,         1 ≤ q < 2
  0,                    q ≥ 2

σ_i = Σ_j m_j W_ij
P_i = c_s² σ_i
a_i,pressure = -Σ_j m_j (P_i/σ_i² + P_j/σ_j²) ∇_i W_ij
```

Pairwise Plummer-softened gravity is added to each acceleration, including the central body's recoil. The integrator applies velocity Verlet to the entire coupled system, with the following prototype parameters:

- Smoothing length `h = 3×10^6 m`
- Gravitational softening `ε = 10^6 m`
- Isothermal effective sound speed `c_s = 5×10^5 m/s`
- Integration step in the default scene `Δt = 0.1 s`

**Momentum and mass:** Mass is explicitly accounted for across all parcels. The force law is pairwise symmetric, so total linear momentum is conserved to floating-point roundoff when no external force is applied. An isothermal pressure closure lacks the thermal energy equation; therefore this implementation makes **no total-energy conservation claim** for the modeled fluid.

**Test gates:** 2D kernel normalization, positive densities, negative radial kernel derivative, net pair force near zero, force direction, zero-step identity, speed bounded over short runs, and numerical momentum/mass drift across the full post-breakup encounter.

**Limitations:** The `q≥1` tidal trigger remains a simplified gravitational strength criterion. Real white-dwarf disruption involves degenerate-matter equations of state, pressure gradients, nuclear burning and radiation; a 2D isothermal model cannot reproduce those outcomes. Debris close to the horizon also requires GR and a capture/accretion treatment. Softening can obscure close-range dynamics, so convergence studies must vary particle count, smoothing length and timestep.

**Method reference:** University of Illinois [SPH educational notes](https://courses.physics.illinois.edu/cs418/fa2026/text/sph.html), with symmetric pressure interactions and compact kernels; see also [Abel (2011), *rpSPH*](https://academic.oup.com/mnras/article/413/1/271/1061267) for discussion of symmetric traditional SPH forces.

## 4. Executables and reproduction

```sh
cargo test --all-targets
cargo run --release --bin strongfield_bench
cargo run --release --bin strongfield_lab
cargo run --release --bin blighted-galaxy
```

`strongfield_lab` shows stable equatorial Schwarzschild orbit, eccentric periapsis precession, radial plunge with chart cutoff, and the 96-particle self-gravitating debris scene. HUD exposes energy constraints, coordinate/proper time, tidal stress, total traced mass and relative momentum drift.

The default game retains SI-calibrated craft motion, with enemy targeting changed to use the causal worldline history. Full visual playtesting on hardware is still required: successful automated build and unit tests establish compilation and algebraic/integrator acceptance gates.

## 5. Remaining experiments

1. Add a horizon-regular geodesic chart and investigate Kerr equatorial geodesics and frame dragging.
2. Replace piecewise-linear source histories with continuous ephemeris plus a causal event queue; implement curved-spacetime photon travel for optical combat.
3. Introduce a 3D polytropic/degenerate EOS, self-consistent pressure and thermal energy, capture treatment, convergence tests and GPU particle rendering.
4. Explore energy-momentum accounting for Wicker's relativistic projectile, the Chamber of Guf and white-dwarf superweapon set pieces as explicit speculative engineering workpapers.
