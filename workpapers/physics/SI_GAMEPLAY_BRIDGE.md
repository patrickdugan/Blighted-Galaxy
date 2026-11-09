# BG-PHY-03: Physical-to-playable simulation bridge

**Status:** Code-backed prototype, 9 November 2026  
**Primary modules:** src/si_bridge.rs, src/main.rs  
**Related documents:** [GDD](../../design/GDD.md), [SI Physics](SI_PHYSICS.md), [Original Design Treatment](../design/game-design-treatment.txt)

## Source premise

The source design treatment describes relativistic jousting, momentum-dominated flight, close transhuman combat, and tidal manipulation with massive objects. *You Never Go Full Blight* records a 2 kg projectile at 0.9993c and roughly 4.6 EJ kinetic energy; *Patricia Tries* dramatizes light-delay combat and the interruption of a Hawking-class payload by feeding it mass. These are fictional story events. The current prototype tests their physical quantities without asserting that the original lost Rust engine has been recovered.

## What runs now

The default game executable launches **two distinct physically scaled gameplay encounters**. The earlier normalized prototype is preserved under bin/arcade_legacy.rs and can be run with cargo run --release --bin arcade_legacy.

### Mode 1: SI relativistic RKM interceptor

- Both craft occupy a region near 1 AU from a solar-mass star, with world positions in metres and coordinate times in seconds. Coordinate units stay f64 until converted to camera-relative pixels.
- The player directs high-energy thrust while its specific relativistic momentum u=gamma*v is integrated. Velocity is reconstructed as v=u/sqrt(1+|u|²/c²), so repeated boosts remain subluminal. The model uses a specified externally powered thrust and Newtonian solar gravity. No fuel mass budget or full relativistic equations of motion in curved spacetime are implemented.
- Each shot has rest mass 2 kg and muzzle velocity 0.72c **relative to its firing craft**. The code performs the 2D special-relativistic velocity transform for the firing craft's lab velocity.
- Projectile tracks use swept-circle collision against moving targets, at 120 Hz. The craft bodies are visually exaggerated; their fictional defended cross-sections are millions of metres wide and are shown as rings. A real 2 kg pellet is tiny in comparison.
- Damage is calculated from the projectile's kinetic energy in the target rest frame, with gamma_rel=gamma_p gamma_t(1-v_p·v_t/c²). The fraction coupled to fictional shields is explicitly game tuned (90%). The model has no atmosphere, matter penetration, beam radiation-transfer, shield constitutive law or blast-damage calibration.
- A ghost marker draws the target's estimated **light-cone apparent position** using a constant-velocity retarded emission solver. This is a local extrapolation, not a complete historical trajectory from an event queue.
- A head-up display reports range, light time, speed in units of c, Lorentz gamma, accumulated craft proper time, impact energy and target shielding.

### Mode 2: Black-hole / white-dwarf tidal disintegration

- A pair-symmetric velocity-Verlet integrator advances an approximately 10 solar-mass nonspinning black hole and 0.65 solar-mass white dwarf, initially separated by about 146 million metres.
- The calibrated near-side differential gravitational acceleration is divided by the white dwarf's Newtonian surface gravity. When q >= 1, the code fires a **proxy** tidal disruption event.
- The white dwarf splits into 96 deterministic mass-labelled tracers. Symmetric initialization preserves their combined mass, centre of mass and initial momentum within floating-point tolerance. An illustrative dispersal impulse of ±25 km/s radial and ±9 km/s tangential is imposed on matched pairs. Stellar hydrostatic pressure, radiation, nuclear burning, realistic partial stripping, and tracer self-gravity are future work.
- After fragmentation, the black hole is moved ballistically at its existing velocity while debris is propagated as ballistic Newtonian tracers in its field. This is an approximation. Conservation claims are limited to the **instant of tracer initialization**; the continued tidal simulation does not model momentum exchange with the black hole.
- Compact-object horizon drawn on the screen is exaggerated; all body marker sizes are visual abstractions.
- The simulation uses a physical dt <=0.10 s, with user-controlled time rate and camera scale.

## Units and adapter contract

The SI layer has explicit wrappers for Metres, Seconds, Kilograms and Joules at selected API and display boundaries, while vector state uses fields suffixed _m and _m_s. WorldCamera subtracts its f64 world-origin in metres before casting to f32 pixels. With an encounter centred at 1 AU, small local offsets remain numerically resolvable.

The arcade engine is maintained separately. In particular, do not accidentally pass positions from the old 1 unit ~= 1 pixel world to the SI camera.

## Regression gates

Run:

    cargo test --all-targets
    cargo test si_bridge -- --nocapture
    cargo check --all-targets
    cargo run --bin physics_bench

The bridge-specific tests cover AU conversion, screen coordinate round-trip at 1 AU, both parallel and transverse Lorentz velocity addition, high-thrust speed limits, proper time, kinetic energy relative to a target, swept collision at speed, shield energy consumption, paired tracer mass and momentum, and first tidal break threshold.

## Reproducible manual validation

1. Launch the Windows artifact blighted-galaxy.exe or cargo run --release.
2. In Mode 1, aim at the ship marked in orange and hold Space. Confirm projectiles follow the cursor, the shield ring loses energy after impacts, and the HUD reports the impact energy in joules.
3. Aim toward the yellow light-cone marker, then vary thrust to see the distinction between observed and simulated present position.
4. Press 2. Observe the white dwarf infall, track the q value, and verify the deterministic breakup at q >=1 into mass-carrying tracers.
5. Increase or reduce the time rate using E/Q and zoom using the wheel or X/Z. Reset with R. Return to Mode 1 with 1.
6. Confirm that the legacy binary arcade_legacy.exe still runs the initial normalized prototype.

The game still needs an actual graphical playtest on a Windows GPU. CI validates builds, numeric regressions and benchmark output, while the interactive effects require display inspection.

## Next work

- Replace the local retarded-position extrapolation with an ephemeris/causal event ledger and beam-time interception.
- Add a **physical unit bridge for close-range nanoborg boxing**, including structure scale, matter handling, strength-limited collisions and energy-based regeneration.
- Replace tracer breakup with 2D SPH or a validated simplified tidal-stream model, including self-gravity and work done by tidal forces.
- Derive strong-field trajectories from Schwarzschild/Kerr geodesics near an engineered compact object and restrict the Newtonian solver to an explicit weak-field validity region.
- Implement GPU particle buffers and scene-driven astrophysical light transport. Bench on the target RTX 5080 before raising swarm counts.
