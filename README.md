# Blighted Galaxy

A recovery project for the lost Rust-based Blighted Galaxy game, preserving its original fiction and design documents alongside a new playable prototype.

The earlier local engine was lost in 2026. This code is a clean-room reconstruction. The surviving 94 pages of source material, approximately 34,000 whitespace-delimited words, are archived as text transcriptions in [recovery/](recovery/README.md).

## Strategic campaign (BG-STR-03)

The galactic map and SI physics encounters now share a deterministic strategic campaign. The Jenkins Swarm autonomously sends 0.10c invasion packets to reachable industrial systems. The Drummers build industry and fabricate contained kugelblitz black holes, then spend them on **0.76c Halo strikes**. Weena can recruit distant systems, prepare industrial Guf relays, and fire a synchronized three-node energy array. Every strike and capture occurs after finite travel time.

**Controls:** click a named system; `I` builds industry, `K` fabricates a kugelblitz, `Y` launches a Halo strike against a hostile system, `V` recruits for Weena, `U` builds a Guf relay, `P` fires the Guf array, `J` issues a Blight attack order, `N` toggles autonomous Blight expansion, and `T` shows/hides strategic rings. Use `Space` to run the historical clock. `1` opens the SI RKM interception and `2` opens the white-dwarf tidal encounter. When your ship has traveled to the selected star system, the encounter is a **campaign mission**, and `3`/`Tab` records victory, defeat, or withdrawal. Distant selections still open sandbox encounters, without changing strategic ownership. A ship travels at 0.76c after clicking a destination and pressing `Enter`.

**System-specific missions:** The destination's industry, defenses, and stored kugelblitz assets now determine RKM enemy shield strength, opponent fire cadence, encounter distance and available friendly shield support. White-dwarf interdiction selects a stellar-remnant mass and dwarf conditions from the same system. During a campaign tidal mission, `A` and `D` spend limited radial tug burns toward/away from the black hole; a 300-second objective deadline gives that mission a defeat state. Each scored campaign encounter triggers a 250-galactic-year rearm period at its system. A map-side briefing displays the target's current battle conditions. `R` restarts the active mission with its original launch configuration.

Battles affect local shields, hostile control, and industry. The infrared-black Blight texture responds to ownership changes, with small reclaimed or defended enclaves around industrial sites, while the original 0.10c Orion envelope remains an absolute boundary. Strategic resource units and kugelblitz engineering are speculative game rules. See [strategic campaign design and test contracts](workpapers/strategy/STRATEGIC_CAMPAIGN.md).

## Playable galactic overworld (procedurally rendered)

The game opens on a **procedural Milky Way navigation atlas**. Its starfield, luminous spiral arms, dust, barred Galactic Center, infrared-darkened Jenkins Swarm, red occupied-star speckles and narrow expansion tendrils are generated from fixed numeric seeds in Macroquad. The image is rendered at runtime from world coordinates, with no baked NASA map, screenshot texture or copyrighted game assets.

The radial navigation stencils match the astronomical chart in the design reference. **Galactic longitude** is measured from the north/up axis clockwise in 30-degree sectors, with 15,000-light-year rings around the Galactic Center. A separate **Orion front** stencil has 15-degree sectors and 2,000-light-year radial bands. The maximal causal front advances at 0.10c and reaches **exactly 10,000 light-years after 100,000 years**; sector-specific Jenkins tendrils stop inside that radius as post-human resistance delays replication by thousands of years. Substantial portions of the galaxy remain visibly unoccupied.

The campaign timeline represents a parameterized *100,000 years since the Orion outbreak* scenario. The recovered [macro-history](workpapers/lore/macro-history.txt) describes a longer and partly undated human/post-human chronology; no historical episode dates are silently asserted. Schematic navigable sites include Amicus, Olympus, the Great Library, the Outer Arm exodus, Bel-Air, the Nebular Wastes, the Perseus front and the Canis Major departure. Each derives its name and narrative role from the source manuscripts. Their approximate map coordinates are gameplay staging, not recovered canon.

~~~sh
cargo run --release
# Previous direct-to-intercept startup:
cargo run --release -- --encounter
~~~

**Overworld controls:** left-click a site or sector to plot a route; **Enter** commits a 0.76c relativistic leg. The navigation clock tracks galactic-frame and ship proper years separately. **Space** pause, **Q/E** adjust years per second, **,/.** skip ±1,000 years, **Home/End** jump to 0/100,000 years. **WASD/arrows** or right-mouse drag pan; wheel or **Z/X** zoom; **F** focus on the Orion wavefront, **C** Galactic Center, **O** Orion; **G** Galactic Center grid, **L** local radial stencil, **B** occultation overlay, **H** help. **1** enters the existing relativistic intercept, **2** the tidal encounter, and **3/Tab** returns to the overworld. **Esc** quits.

See the [technical and art-direction workpaper](workpapers/overworld/PROCEDURAL_ORION_NAV.md) for deterministic generation, sector math, travel timing and source/canon distinctions.

## Main game: SI-driven encounter bridge

The encounter bridge contains a physically scaled **relativistic intercept** and an opt-in **white-dwarf tidal encounter**. These scenes are accessible from the galactic overworld with keys **1/2**, or can be launched directly with `--encounter`. The original game remains preserved as an opt-in legacy binary.

~~~sh
cargo run --release
cargo run --release --bin arcade_legacy
~~~

SI game controls: **1** relativistic intercept, **2** white-dwarf tidal encounter, **mouse** aim, **M** toggle keyboard aim, **A/D** rotate with keyboard aim, **W/S** thrust, **Space** release 2 kg kinetic projectiles, **Q/E** reduce/increase simulation time rate, **Z/X** zoom, **R** restart, **Esc** quit.

The player and enemy craft move with SI coordinate positions in metres and special-relativistic momentum-to-velocity conversion. Every shot transforms a 0.72c muzzle velocity out of the firing craft's rest frame, carries a computable kinetic energy in joules and uses swept continuous collision. Defensive envelopes and energy coupling remain **fictional game mechanics**. The second encounter calculates a white dwarf's tidal-acceleration/self-gravity ratio near a 10-solar-mass black hole. When the ratio passes one, deterministic debris tracers inherit the white dwarf's mass and momentum. The breakup and debris are an explicitly identified proxy in place of full stellar hydrodynamics.

See [SI gameplay workpaper](workpapers/physics/SI_GAMEPLAY_BRIDGE.md) and [physics research](workpapers/physics/SI_PHYSICS.md).

## Legacy arcade sandbox

The original 2026 clean-room normalized-physics sandbox remains available:

~~~sh
cargo run --release --bin arcade_legacy
~~~

The legacy scene contains a gravitationally active star, a black hole, a white dwarf, a neutron star, a nanite swarm and a dueling post-human opponent. It uses GPU-backed 2D drawing via macroquad.

Controls: W/S thrust, A/D rotate, Space fire, E contact strike, B release a short-lived Hawking payload, G scatter a swarm, mouse wheel zoom, R restart, Escape exit.

For a Windows build without installing Rust, open the latest successful [GitHub Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/ci.yml) and download the artifact named blighted-galaxy-windows-x64. Unzip it and run blighted-galaxy.exe for the new SI encounter game, or arcade_legacy.exe for the original normalized sandbox. The artifact is generated by the Windows release job and retained for 30 days.

## Canon and game design

- [Game Design Document](design/GDD.md): narrative/technical product spec and milestones.
- [Stories](stories/README.md): four historical fiction manuscripts, with recovered texts retained.
- [Workpapers](workpapers/README.md): game-design treatment, macro-history and calibrated physics workpaper.

## Calibrated physics lab

For a genuine SI-unit orbital demo, run:

~~~sh
cargo run --release --bin si_orbits
~~~

The laboratory has a Sun/Earth regression binary, a white-dwarf flyby and an eight-solar-mass black-hole binary. Time controls and diagnostics are shown on screen. Body markers are visually exaggerated, while trajectories are evolved in metres, seconds and kilograms by the separate SI reference module.

~~~sh
cargo run --release --bin physics_bench
cargo test --all-targets
~~~

The benchmark prints reference values for Hawking radiation, a near-light-speed 2 kg kinetic projectile, time dilation and orbital energy conservation. See [SI Physics Workpaper](workpapers/physics/SI_PHYSICS.md) for approximations and scope.

## Strong-field geodesics, causal observations and self-gravitating debris

The BG-PHY-04 physics milestone adds a **Schwarzschild geodesic simulator** for test particles orbiting nonspinning black holes, a recorded causal worldline solver for delayed observations, and a symmetric pair-force 2D SPH-like debris simulation with central-body recoil.

~~~sh
cargo run --release --bin strongfield_lab
cargo run --release --bin strongfield_bench
~~~

The graphical lab has four modes: stable circular orbit, precessing eccentric orbit, radial infall toward a Schwarzschild coordinate-chart cutoff, and self-gravitating white-dwarf disruption. Switch between them with 1–4; Space pauses, Q/E changes simulation speed, mouse wheel or Z/X zooms, and R resets. It displays proper/coordinate time, geodesic energy residual, and fluid mass/momentum diagnostics.

The **default combat executable** now obtains observed positions from bounded recorded histories. It withholds targeted AI fire until causal data is available, rather than extrapolating nonexistent observations.

Scientific limitations: Schwarzschild geometry describes a nonrotating isolated central mass with negligible test-particle backreaction. The 2D isothermal fluid is a toy proxy for degenerate white-dwarf matter and does not model nuclear reactions or radiation. See [BG-PHY-04 workpaper](workpapers/physics/STRONG_FIELD_CAUSAL_SPH.md) for the equations, assumptions, and acceptance tests.

## Kerr photons and 3D stellar fluid (BG-PHY-05)

The next standalone lab implements **equatorial Kerr null photon trajectories**, spin-dependent frame dragging, a **3D polytropic SPH-like fluid** after the white-dwarf's tidal breakup, and GPU-backed batched particle sprites.

~~~sh
cargo run --release --bin kerr_sph3_lab
cargo run --release --bin kerr_sph3_bench
~~~

Controls: **1** Schwarzschild photon rays, **2** spinning Kerr (+0.85), **3** reversed spin (-0.85), **4** 3D stellar debris. Q/E changes time rate; Z/X or wheel zooms, A/D rotates the debris projection, R resets.

The benchmark reports Kerr orbital reference values, radial constraint residuals, 3D fluid energy and momentum drift, and a 10,000-particle batch test. The GPU renders indexed triangles; SPH and mesh assembly still run on the CPU. The physically grounded modules remain separate from the fictional energy shields, near-horizon weapons and other source-world engineering.

Detailed assumptions: [BG-PHY-05 physics workpaper](workpapers/physics/KERR_PHOTON_SPH3.md).

## Non-equatorial Kerr photon camera and cold white dwarfs (BG-PHY-06)

The new Kerr camera integrates past-directed null trajectories with the Carter separation constant from a local zero-angular-momentum observer. It resolves black-hole silhouettes, curved photon paths, synthetic background stars, and a thin accretion disc at arbitrary viewing inclination.

The new white-dwarf module solves Newtonian hydrostatic structure with the fully relativistic zero-temperature electron equation of state. A 0.65-solar-mass carbon–oxygen reference star has a calculated radius around 8,402 km. Its radial density profile can seed a 3D SPH model using the same cold electron pressure.

~~~sh
cargo run --release --bin kerr_camera_lab
cargo run --release --bin white_dwarf_lab
cargo run --release --bin kerr_camera_bench -- --output kerr-camera.png
cargo run --release --bin white_dwarf_bench -- --output wd-mass-radius.csv
~~~

The camera lab uses keys 1–4 to select spin (0, +0.45, +0.85, -0.85), Up/Down for inclination and Z/X for field of view. The stellar lab uses keys 1–3 for mass, Space to pause, Q/E for time rate and A/D to rotate the 3D projection.

The curved Kerr rays and hydrostatic white-dwarf profiles use calibrated equations. Synthetic disc colors, idealized 3D particle sampling, isothermal assumptions in older scenes and fictional spacecraft hardware remain distinguishable from physical predictions.

Full derivations and test scope: [BG-PHY-06 workpaper](workpapers/physics/KERR_CAMERA_WHITE_DWARF.md).

## Curved Kerr optical acquisition and white-dwarf relaxation (BG-PHY-07)

A new **moving-target optical laboratory** sends past-directed Kerr photon rays through the recorded worldline of a posthuman combatant on a circular test-body orbit. It computes the emission time and the source position that an observer receives after the curved photon flight. The source's simulated present position is displayed separately.

~~~sh
cargo run --release --bin kerr_optical_lab
cargo run --release --bin optical_relax_bench
~~~

Kerr optical lab controls: **1/2/3** choose black-hole spin (0, +0.75, -0.75), **F** toggles automatic optical acquisition, **A/D** steer the manual ray, **Q/E** change time rate, and **Z/X** or mouse wheel zoom. The sensor traces **incoming** photons and targets the recorded causal past. It does not simulate an outgoing laser or full 3D moving-object intersection.

The white-dwarf lab also supports **H** to reinitialize a hydrostatic star and apply 15 simulated seconds of isolated damped SPH relaxation before entering the black-hole encounter. The new diagnostics track force reduction and density mismatch. Resolution tests reveal that the relaxed low-particle-count star contracts excessively, despite large reductions in force residual. Treat the relaxed star as a numerical experiment until the mass–radius convergence gate passes.

More details and measurements: [BG-PHY-07 workpaper](workpapers/physics/KERR_CAUSAL_RELAX.md).

## BG-PHY-08: orbiting Kerr observers and adaptive white-dwarf SPH

The new **moving-observer Kerr lab** carries an equatorial camera along a circular timelike Kerr orbit, computes its speed relative to the local zero-angular-momentum frame, applies photon aberration, and solves light-cone intersections with a **recorded, potentially noncircular source history**. Its target envelope uses a local metric-based estimate of the emitter's rest space.

~~~sh
cargo run --release --bin kerr_moving_lab
cargo run --release --bin phy08_bench -- --csv bg-phy08-adaptive.csv
~~~

Moving-observer controls: **1–3** black-hole spin, **W** switching between a ZAMO and an orbiting observer, **F** automatic photon acquisition, **A/D** manual ray direction, **Q/E** clock rate, **Z/X** zoom.

The white-dwarf lab adds **J** for isolated damping using an adaptive smoothing-length and grad-h-corrected SPH pressure force; **H** retains the original fixed-h comparison. Following relaxation, each mode uses the corresponding force method for live evolution. The new solver increases stellar radius at the measured particle counts. A 64-particle undamped hold still expands by over 12% in three dynamical times, leaving [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) open.

Read the complete [BG-PHY-08 physics workpaper](workpapers/physics/ADAPTIVE_SPH_MOVING_OBSERVER.md) for equations, CI data, convergence targets and the metric-worldtube approximation.

## BG-PHY-09: matched variable gravity and shell-balanced white dwarfs

The latest research mode includes a density-dependent **Plummer gravitational softening** with the corresponding grad-h energy correction, plus a deterministic radial-shell-balanced particle initializer inspired by stellar-SPH initial-condition methods. The code preserves each previous model for measurement.

~~~sh
cargo run --release --bin phy09_bench -- --csv bg-phy09-comparison.csv
cargo run --release --bin white_dwarf_lab
~~~

The white-dwarf lab controls now include **H** for the old fixed-h solver, **J** for adaptive hydrodynamic smoothing with fixed gravitational softening, and **K** for the new shell-balanced seed with variationally differentiated adaptive Plummer gravity. Each mode uses the corresponding gravitational model for subsequent dynamics. The new bench compares radius, density-profile error, force residual and an undamped hold.

The short-time stellar radius is closer to the independent hydrostatic reference in some low-resolution cases. The **long undamped hold still exhibits large expansion**, and at N=512 the earlier solver has a slightly smaller radius bias. [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) remains open.

See the [BG-PHY-09 workpaper](workpapers/physics/GRAD_H_GRAVITY_SHELL_SEED.md) for the matched energy/potential equations, baseline comparisons and convergence limits.

## BG-PHY-10: Poisson-matched M4 gravity and nearly cubic white-dwarf shells

The new gravity model derives its softened potential from **Poisson's equation for the same cubic-spline density kernel** used by the SPH solver. It becomes exactly Newtonian beyond 2h and includes a density-adaptive gravitational grad-h term, with the matching energy diagnostic. A new hydrostatic particle initializer chooses radial shell widths comparable to their angular tile widths.

~~~sh
cargo run --release --bin phy10_bench -- --csv bg-phy10-comparison.csv
cargo run --release --bin white_dwarf_lab
~~~

The stellar lab's **L** key selects nearly cubic shells with Poisson-matched M4 gravity. Earlier H/J/K modes remain available. The headless benchmark compares initializations and gravity models from 64 to 1,024 particles, and tests selected configurations over three dynamical times without damping.

At 1,024 particles, the 0.65-solar-mass M4 star has a radius of 7,999 km versus the 8,402 km hydrostatic reference and a density log-RMS error of 0.1667. This short run lasts only 0.6 s, approximately 0.23 stellar dynamical times. The 128- and 256-particle undamped stars still contract more than 10%; **Issue #8 remains open** until stable long-duration equilibrium is demonstrated.

See the [BG-PHY-10 physics workpaper](workpapers/physics/M4_GRAVITY_CUBIC_SHELLS.md) for the exact kernel potential and force formulas, source references, numerical tables and acceptance gates.

## BG-PHY-11: White-dwarf radial profile and undamped stability

A deterministic **SPH spatial cell neighbor index** now skips known-zero cubic-spline density and grad-h interactions while preserving the prior physical equations. Neighbor candidates are sorted into stable particle-index order and the radius is updated as smoothing lengths evolve. Long-range gravity remains all-pairs in this milestone.

The new radial mass-profile evaluator reports R50, R90, R95, R99, maximum radius, and RMS deviation from the independently integrated hydrostatic mass CDF. The headless benchmark tracks these throughout a multi-dynamical-time undamped hold.

~~~sh
cargo run --release --bin phy11_bench -- --n 128 --relax-steps 150 --hold-steps 262 --csv bg-phy11-quick.csv
cargo run --release --bin phy11_bench -- --n 1024 --relax-steps 300 --hold-steps 788 --csv bg-phy11-highres.csv
~~~

A separate [BG-PHY-11 high-resolution GitHub Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy11.yml) produces a 1,024-particle numerical artifact. Its physical acceptance verdict is recorded independently of software build success. The star continues to require a 5%-radius and 5%-undamped-stability pass before it can anchor calibrated stellar destruction mechanics.

Full methods: [BG-PHY-11 physics workpaper](workpapers/physics/NEIGHBOR_MASS_PROFILE_HOLD.md). Convergence tracker: [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

## BG-PHY-12: Pressure/gravity balance and tapered relaxation

The cold white-dwarf research engine now calculates **outward pressure and inward self-gravity separately** in enclosed-mass shells. Its diagnostics show the local hydrostatic gravity reference, radial and tangential residual acceleration, an SPH force-virial ratio, and the stellar mass-radius profile throughout undamped evolution.

Two relaxation schedules can be compared with identical particle seeds and physical time: the previous constant velocity damping and a new strong-to-weak geometric taper that turns damping off near the end of the preparation stage.

~~~sh
cargo run --release --bin phy12_bench -- --n 128 --relax-steps 300 --hold-steps 262 --csv bg-phy12-timeseries.csv --shell-csv bg-phy12-shells.csv
~~~

A dedicated [BG-PHY-12 GitHub Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy12.yml) tests the 512-particle reference over 6 seconds of relaxation and 7.88 seconds of undamped evolution for both schedules. The resulting CSV files keep numerical conservation checks and the physical hydrostatic-equilibrium verdict separate.

See the [BG-PHY-12 workpaper](workpapers/physics/RADIAL_FORCE_RELAX_SCHEDULE.md) and open [stellar-equilibrium Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

## BG-PHY-13: smoothing-length and stellar resolution experiment

The white-dwarf physics laboratory now has a separate benchmark that varies the SPH smoothing parameter **eta**, while preserving the same cold-electron equation of state, M4 Poisson-matched gravitational potential, shell-balanced particle seed, damping duration and undamped hold. It evaluates actual kernel-support neighbor counts, local SPH gravity against the independent hydrostatic field, radial density errors and conservation.

~~~sh
cargo run --release --bin phy13_bench -- --n 256 --etas 0.85,1.05,1.20,1.45 --relax-steps 300 --hold-steps 788 --csv bg-phy13.csv --shell-csv bg-phy13-shells.csv
~~~

The [BG-PHY-13 GitHub Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy13.yml) performs matched-duration tests at N=256, 512 and 1,024, followed by approximately three stellar dynamical times without damping. Results and physical acceptance flags are exported per resolution. A smaller N=64 smoke test is included in ordinary CI.

**Measured outcome:** N=256–1024 all completed their three-dynamical-time holds and all tested eta settings failed at least one stellar equilibrium gate. At N1024, eta=0.85 gives a post-damping radius 8.78% below hydrostatic and a 5.77% maximum outer-radius excursion, despite 4.07% RMS error against the independent gravity field. The playable stellar modes therefore retain the previous constant-damping baseline. Complete methodology, data tables and source references are in the [BG-PHY-13 workpaper](workpapers/physics/ETA_HYDROSTATIC_SWEEP.md), with [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) tracking unresolved stellar equilibrium.

## Simulation modules

- Reproducible seeded nanoswarm scatter.
- Newtonian N-body dynamics with softened gravity and velocity Verlet integration.
- Swept-circle projectile collision, with AABB contact-boxing.
- Simplified tidal-disruption criterion and visible shear debris.
- Analytic normalized evaporation law with end-of-life particle burst.
- Fixed 120 Hz simulation step with independently rendered frames.

The legacy arcade units remain normalized. The default physical encounter has SI quantities, while strong-field general relativity, radiation transport, fluid disruption and exotic defensive materials remain unimplemented.

## Rights and provenance

The source fiction and historical design materials remain the author's work, with rights reserved pending a chosen license. The reconstructed Rust code has no distribution license yet and is marked non-publishable in Cargo. See [RIGHTS.md](RIGHTS.md).

## Evidence

The [Game Design Treatment](recovery/design/game-design-treatment.txt) calls for relativistic jousting, transhuman duels and interstellar war economics. The [macro-history](recovery/design/macro-history.txt) outlines 22 episodes. Four recovered fiction manuscripts supply encounters and characters. All six transcriptions preserve the text extracted from the user-supplied Word documents.

See [engine recovery notes](docs/ENGINE_RECOVERY.md) for the distinction between documented game design, remembered features of the lost engine, and current reconstruction.

## Tests

~~~sh
cargo test --all-targets
cargo check --all-targets
~~~

GitHub Actions runs the same checks on pushes and pull requests. The core math tests include high-speed projectile tunneling, tidal trend, Hawking step consistency, deterministic swarm generation and velocity bounds.
