# BG-STR-01 | The industrial war over a darkening galaxy

This layer extends PR #32's galactic navigation and its SI encounter modes. It uses the same `GalacticNavModel.year`, schematic waypoints, Orion origin, and 0.10c geometric limit.

## Play

Start the game, select a labeled star system, and use:

| Key | Action |
| --- | --- |
| T | Toggle strategic site rings, industrial icons, and in-flight attack links |
| I | Build industrial capacity using material and energy |
| K | Manufacture a kugelblitz/contained Halo Drive seed at a Drummer industrial site |
| V | Send a delayed Weena recruitment/industrial takeover order to the selected Drummer or independent system |
| U | Prepare a Guf relay node at a Weena-controlled industrial system |
| P | Fire the Guf array at the selected system, requiring three powered relays and a full charge |
| J | Send Blight replication toward the selected system from the nearest captured source |
| Space, Q/E | Run/pause strategic calendar and adjust speed |
| 1/2 | Enter the existing SI intercept and tidal-breakup encounters; 3/Tab returns to the overworld |

The strategic screen uses the existing waypoint geometry. An additional Orion seed is the initial Blight industrial source. A cyan auxiliary dot denotes stored kugelblitz holes; gold rings denote Guf relays. Lines indicate packets currently moving at finite speeds. The right-hand system dossier reports industry, energy, materials, stored holes, array status and packet count.

## Simulation contract

- Strategic economics are fixed game units. Production accrues in deterministic 250-year substeps, with bounded inventories.
- Kugelblitz fabrication requires a sufficiently industrialized Drummer system and stored power/material. This is explicitly a fictionalized strategic cost proxy; actual creation and confinement of kugelblitz black holes is speculative.
- Weena can recruit industrial systems, with coordination taking distance in light-years in years plus a scenario-defined 2,000-year political/technical delay. Recruitment converts industrial ownership only on arrival.
- The Chamber of Guf requires three Weena-controlled, charged industrial nodes. Each node pays its firing cost. The attack reaches the target no earlier than the longest individual emitter-to-target light path. This is the gameplay simplification of Episode XIX's galactic attack cone and Episode XX's delayed responses.
- Blight packets propagate at 0.10c plus 1,500 years of resistance, with an additional global Orion causality check. Existing resistance notches and the generated black texture remain driven by PR #32's fixed history model. The strategy rings/packets display changed ownership; rewriting the actual darkened star raster from the live strategic occupation map is a separate graphics pass.
- Player time jumps rebuild a deterministic strategic initial state at the target year. This is a scenario restart, not an undoable history of prior orders. Ship coordinate/proper-time calculations remain independently handled by the navigation model.

## History mapping

Episode XVII's hidden planetary accelerator establishes the Remnant's stored kugelblitz drives. Episode XVIII's O-star Dyson construction establishes the industrial background for hostile black-hole weapons. Episode XIX grounds Weena's supergiant/12-solar-mass remnant and multi-system Chamber of Guf. Episode XX grounds long light-time delays, counterfire from individual swarms, and the need for anticipation or evasive motion.

Site alignments, project prices, system ownership and action rates are balance assumptions. There is no assertion that the texts supply exact star coordinates.

## Next implementation milestones

1. Derive tactical mission parameters and rewards from the selected strategic encounter, and return battle results to the originating star system.
2. Render the true current strategic occupation mask into the Blight texture, applying the causal boundary and resistant frontline per pixel.
3. Add industrial convoys, dynamical star relocation, multi-beam interception, and save/load for deterministic replay.
4. Replace three-node threshold with a constrained synchronization and energy-routing solver for the Chamber, including beam width and finite steering.
5. Build scripted Episode XVII through XXII campaigns and an opposing planning AI.

## Verification

Unit tests cover scarce resource production, array node requirements, emitter-to-target light travel, Weena recruitment latency and Orion causal bounds. The dedicated GitHub workflow runs `cargo test --all-targets` and `cargo check --bin blighted-galaxy` on the PR.
