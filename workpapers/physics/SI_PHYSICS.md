# SI Physics Workpaper — BG-PHY-02

**Date:** 9 November 2026  
**Status:** Implemented reference core and unit tests in `src/physics.rs`. The existing playable scene still uses its explicitly labeled normalized game units. This is an engineering research layer that can be tested before scaling the renderer.

## 1. Source material versus mathematical extensions

The historical Game Design Treatment asks for gravity near massive objects, time-lagged shots, Hawking-class black-hole bombs, stellar shearing and precise transhuman duels. The fiction makes use of white dwarfs, neutron stars, stellar war economies and relativistic weapons. The equations below are modern physical reference models used to reason about those fictional systems; no numerical values in this workpaper are statements from the authorial story texts.

Sources for numerical constants and approximations:

- NIST 2022 CODATA recommended constants: https://physics.nist.gov/cuu/Constants/
- NIST extended constants listing: https://physics.nist.gov/cuu/pdf/all.pdf
- University of Zurich Hawking-radiation course notes, sec. 6.3: https://www.physik.uzh.ch/groups/jetzer/notes/PHY529_ART_III.pdf
- Naive lifetime model is affected by greybody factors and species channels. Example discussion: https://academic.oup.com/mnras/article/532/1/451/7656972

## 2. Units and constants

| Symbol | Value | Unit |
|---|---:|---|
| c | 299,792,458 (exact) | m/s |
| G | 6.67430 × 10^-11 | m³ kg^-1 s^-2 |
| ħ | 1.054571817 × 10^-34 | J s |
| k_B | 1.380649 × 10^-23 (exact) | J/K |
| Astronomical unit | 149,597,870,700 (exact) | m |
| Julian year | 31,557,600 | s |
| Solar mass | 1.98847 × 10^30 (representative) | kg |

The game must specify a **coordinate frame**, plus fixed metric conversion per rendering scene. Unit tests should fail on accidental mixing of SI, astronomical units, and game units.

## 3. Gravity and orbital numerical method

For bodies outside strongly relativistic regimes, acceleration on i from j:

`a_i = G m_j (r_j-r_i) / (|r_j-r_i|²+eps²)^(3/2)`.

Pairwise accumulation updates both bodies so action and reaction agree to roundoff. The semi-symplectic **velocity-Verlet** update uses a beginning and ending acceleration and is reversible at fixed dt:

`r_next = r + v dt + (1/2) a_old dt²`  
`v_next = v + (1/2)(a_old+a_new) dt`

Potential energy is `-G m_i m_j/sqrt(r²+eps²)`, matching the chosen force law. The module also records total momentum and mechanical energy. Two-body trajectories should close after their Keplerian period `P=2πsqrt(a³/(G(M+m)))` within a specified discretization error.

This module uses a **fixed maximum substep**. It remains inaccurate for extremely close stellar encounters and all strong-field black-hole orbits. Add adaptive timesteps or a dedicated close-encounter solver when tests require them.

## 4. Relativistic flight, light cones and kinetic energy

For a fixed subluminal speed v:

`beta=v/c; gamma=1/sqrt(1-beta²)`  
`tau=Delta_t/gamma`  
`E_kin=(gamma-1)mc²`

The code computes gamma-1 with `beta²/[sqrt(1-beta²)(1+sqrt(1-beta²))]` to preserve accuracy for small v. Proper time refers to a *constant-speed segment*, not a ship accelerating through changing fields.

For collinear velocities `u` and `v`, relativistic addition is `(u+v)/(1+uv/c²)`. It remains subluminal when both inputs are subluminal.

For a uniformly moving light source, solve `|x_source(t_emit)-x_observer|=c(t_observe-t_emit)` to display the source at its **retarded emission position**. This supplies the correct foundation for delayed weapon detection and optical HUDs, though gravitational light bending is separate.

**Example:** A 2 kg impactor at 0.9993c has kinetic energy vastly above a Newtonian `mv²/2` estimate. `physics_bench` prints the exact special-relativistic model value in joules and megatons TNT equivalent. No material equation of state or atmospheric blast calculation is included.

## 5. Black-hole quantities

For an isolated Schwarzschild black hole:

`R_S = 2GM/c²`  
`T_H = ħc³/(8πGMk_B)`

The basic constant-species emission approximation is:

`alpha = ħc⁴/(15360πG²)`  
`dM/dt = -alpha/M²`  
`P = alpha c²/M²`  
`t_life = M_initial³/(3alpha) = 5120πG² M³/(ħc⁴)`

At 1 solar mass, the Schwarzschild radius is approximately **2.953 km**, the temperature approximately **6.17×10^-8 K**. The naive lifetime is enormous; small holes evaporate quickly because of the cubic mass dependence.

**Precision and physics caveats:** Greybody factors, spin, charge and opening particle species change luminosity. Accretion offsets evaporation. Near the endpoint, semiclassical theory loses reliability. Rendering a bright terminal event is a **fictional game design choice** with an energy-budget upper bound based on the remaining mass; the code doesn't assert an experimentally established explosive endpoint.

**Analytic mass integration:** `M(t)=M0[1-t/t_life]^(1/3)` for `0≤t≤t_life`; returns zero after assumed endpoint as a simulation convention.

## 6. Shearing and lensing

For a target of radius R at distance d from a disruptor:

`a_tidal,near = GM[(d-R)^-2 - d^-2]`  
`a_self = Gm/R²`  
`q=a_tidal,near/a_self`

The ratio is a **diagnostic** when bodies are non-overlapping. High q suggests vulnerability to tidal breakup, but rigid material strength, spin, density gradient and impact mechanics determine real failure.

Weak-field deflection by impact parameter b is `alpha_lens≈4GM/(c²b)`. The module refuses to return this approximation if `b<20R_S`, where strong-field general relativity must be considered.

## 7. Regression gate

- 2022 CODATA constants and SI scale checks
- Solar Schwarzschild radius and Hawking temperature
- Lifetime `~M³`, Hawking power `~1/M²`, energy conservation in the analytic evaporation prescription
- 0.8c gamma and proper time; 0.5c kinetic energy; relativistic addition and speed checks
- One light-year propagation and retarded moving-emitter position
- Monotonic near-side tide and weak-lensing behavior
- Symmetric N-body acceleration, one Sun/Earth orbit, mechanical energy/momentum drift, and zero-step invariance

Run `cargo test --all-targets` or `cargo run --release --bin physics_bench`.

## 8. Integration plan

1. Retain **two numerical scales** while deriving the coordinate transform. World x/y in AU (map), local encounter meters (physical), camera pixel coordinates (render).
2. Introduce strongly typed wrappers for `Meters`, `Seconds`, `Kilograms` and `Joules` at the adapter boundary.
3. Attach a retarded optical state to every simulated emitting body and show it separately from inferred present position.
4. Advance real bodies with SI time stepping; use an encounter-window scheduler for near-contact objects.
5. Introduce density profiles and smoothed-particle hydrodynamics only after showing that the 2D tidal proxy is insufficient.
6. Replace CPU sprite scattering with a GPU particle buffer and a deterministic CPU comparison harness.

This is an explicit stage gate: **a passing SI bench does not validate the current arcade scene as fully physical**.
