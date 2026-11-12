# BG-STR-05 | Reproducible campaigns and finite Guf optical geometry

## Player controls

F5 writes a versioned JSON campaign to blighted-galaxy-save.json in the process's current directory. F9 restores that slot after validation. The optional CLI flag --save-file campaign.json selects another path. On successful load the map reconstructs navigation and strategy, clears stale convoy-source selection, and regenerates territorial Blight pixels.

Shift+P sends aiming orders toward the selected star to every powered Weena Guf emitter. The order takes light-travel time to reach each mount, then turns at a limited angular rate. P fires only after three nodes are charged and aligned. Space, Q/E and the timeline controls advance the galactic history while waiting.

Tactical encounters operate in local SI seconds. F5/F9 save the overworld campaign, not an in-progress local battle or display pan/zoom.

## Causal beam geometry

Each emitter i has a SlewMount storing initial bearing, desired bearing, slew start and completion year, and assigned target. A light-speed command originating at Weena's controller reaches the emitter at

\[
t_{\mathrm{receive},i}=t_{\mathrm{order}}+\|x_i-x_{\mathrm{controller}}\|/c,\qquad c=1\ {\rm ly/yr}.
\]

The signed shortest rotation is

\[
\Delta\theta_i=\mathrm{wrap}_{[-\pi,\pi)}
\left(\operatorname{atan2}(y_t-y_i,x_t-x_i)-\theta_i\right).
\]

With game-scale maximum angular rate \(\omega_{\max}=2.5\times10^{-4}\,\mathrm{rad/yr}\),

\[
t_{\mathrm{aligned},i}=t_{\mathrm{receive},i}+|\Delta\theta_i|/\omega_{\max}.
\]

An array shot only launches after at least three emitters align. The existing distributed pulse logic synchronizes arrivals:

\[
t_{\mathrm{impact}}=t_{\mathrm{fire}}+\max_{i\in E}\|x_t-x_i\|/c,
\quad t_{\mathrm{launch},i}=t_{\mathrm{impact}}-\|x_t-x_i\|/c.
\]

Each packet locks its emitter's azimuth at firing time. At arrival, the beam strikes when

\[
|\mathrm{wrap}(\theta_{\mathrm{target}}-\theta_{\mathrm{beam}})|
\le 0.0008+\arctan\!\left(\frac{8\ \mathrm{ly}}{\|x_t-x_i\|}\right).
\]

The eight-light-year target hitbox represents infrastructure dispersed through the game's strategic local sector. The angular half-width, mount speed, omniscient target information and Guf ward capacity are explicitly fictional scenario parameters. They should not be interpreted as a real optical design or solar-system size. The campaign uses static source/target coordinates; regression tests move a target after launch and establish a miss.

Successful defensive Guf escort missions still charge a finite ward-energy buffer. Beams consume that buffer on impact before damaging a system. Later work should add angular ward coverage, live stellar ephemerides and source-side information constraints.

## Event epochs

A 250-year economy substep stops early at an approaching packet or convoy arrival, resolving consequences at the scheduled galactic coordinate year. Autonomous Blight planning stops at exact 2,500-year boundaries. The battle simulation has its own seconds-scale clock and only contributes a result on debrief.

## Save format v2

Source: src/campaign_save.rs. Format identifier: BLIGHTED-GALAXY-CAMPAIGN. Serde JSON with float_roundtrip preserves exact double-precision values for deterministic continuation.

The snapshot contains the galactic year, ship location, planned destination, in-flight relativistic leg, ship proper time, pause state and time rate. Industrial systems preserve control, economy, defenses, stored kugelblitz drives, Guf wards, array nodes, and tactical cooldowns. Convoys preserve escrowed resources, route indices, launch/arrival time and convoy sequence. Attack packets preserve type, direction, power, lock-on azimuth and scheduled trajectory. Each emitter mount preserves target and pending slew. The campaign also saves autonomy, global array charge, territorial revision and recent event logs. Unused tactical timestamps serialize as null.

System names and coordinates remain canonical Rust constants. Each save carries a fingerprint of the canonical map, as well as the schema version. Loading validates that fingerprint, both clocks, finite numbers, resource limits, source/target indices, packet light cones, relative ship travel, Guf mount aiming data, convoy duration, sequence IDs and log bounds. A failed read or invalid file leaves the current campaign intact.

F5 writes a temporary json.tmp file; an existing save is first moved to json.bak. A failed final move attempts to restore the previous file. The backup remains available on success. This is desktop safe-save protection, rather than a fully crash-consistent transactional database.

Schema v2 replaces the experimental v1 format because the relay targeting state and locked beam bearings were introduced after v1.

## Acceptance criteria and commands

- Exact JSON floating-point round trip, including a paused or active relativistic voyage.
- Deterministic continuation of pending Guf slews and in-flight convoys after loading.
- A durable previous-save backup and schema/resource tampering rejection.
- Signed shortest-angle rotation, finite command latency and failure of an off-axis beam to strike a displaced target.
- Exact light-cone arrival epochs, original Orion 0.10c causal bound, and Windows game compilation.
- Existing full atlas, Orion detail and tactical screenshot regressions preserved.

Run: cargo test --all-targets; cargo check --bin blighted-galaxy; cargo run --release -- --save-file campaign.json.

Remaining work includes tactical mission mid-battle persistence, moving targets, angular defense-screen footprints, fog of war and a constrained optimization solver for emitter selection and firing costs.
