# Rust engine recovery: provenance and boundaries

The historical Windows directory was recorded as C:\projects\Blighted Galaxy in April 2026. The original engine was subsequently lost. This repository was empty when reconstruction began on 9 November 2026.

## Historical recollection

The author remembers a Rust engine with sophisticated collision boxing, particle rendering, gravitational interactions between stellar remnants, tidal shearing, transhuman duels, swarm scatter, and Hawking-class detonations. That description is preserved as a functional target, not an assertion that the original code was recovered.

## Documentary basis

The surviving Game Design Treatment specifies three interlocking scales: relativistic jousting, bio-nanoborg close combat, and interstellar war economics. It calls for subluminal flight, light-time delays, gravity-assist movement, shearing during black-hole bomb encounters, and transformed bodies fighting at close range. The stories Patricia Tries, You Never Go Full Blight, Blowing Up The Galaxy, and White Dwarf Super Shotgun provide corresponding encounter scripts. Refer to the text transcriptions in recovery/.

## First prototype contract

- One coordinate plane with normalized gameplay units and a nominal light-speed ceiling.
- Plummer-softened N-body Newtonian gravity, advanced with velocity Verlet.
- Moving continuous swept-circle collision tests for high-speed projectiles; AABB broad-phase for body boxing and punches.
- A first-order tidal to self-gravity ratio that can trigger visible stellar shear.
- An analytic toy evaporation law, dm/dt = -k/m²; the final detonation produces visible debris.
- Seeded repeatable particle/swarm scatter.
- A playable 2D sandbox renderer and an automated unit-test suite.

## Deliberately deferred

A Schwarzschild/Kerr metric, relativistic geodesics, radiation transport, real stellar interior equations of state, stable mass/energy scaling, collision broad-phase acceleration structures, camera-relative relativistic effects, network play, and game content beyond the first duel encounter.

The first milestone is visually plausible with deterministically testable mechanics, and it is explicitly NOT a physically calibrated astrophysical simulator. Scientific invariants should become acceptance gates before expanding gameplay to sector strategy.

## Recovery hygiene

Do not rewrite the six historical stories in place. Record corrections and continuity reconciliations separately. Text transcriptions retain parser line wraps and pagination; original DOCX styling and at least one embedded illustration await a binary archive.
