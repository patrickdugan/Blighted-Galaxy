# BG-STR-03 | Strategic industry and physical encounter missions

This layer extends PR #32's galactic navigation and its SI encounter modes. It uses the same `GalacticNavModel.year`, schematic waypoints, Orion origin, and 0.10c geometric limit.

## Play

Start the game, select a labeled star system, and use:

| Key | Action |
| --- | --- |
| T | Toggle strategic site rings, industrial icons, and in-flight attack links |
| I | Build industrial capacity using material and energy |
| K | Manufacture a contained kugelblitz/Halo Drive seed at a Drummer industrial site |
| V | Send a delayed Weena recruitment/industrial takeover order to the selected Drummer or independent system |
| U | Prepare a Guf relay node at a Weena-controlled industrial system |
| P | Fire the Guf array at the selected system, requiring three powered relays and a full charge |
| Y | Launch a Drummer Halo strike at a hostile system, consuming one stored kugelblitz and traveling at 0.76c |
| N | Enable or disable autonomous Blight expansion |
| J | Send Blight replication toward the selected system from the nearest captured source |
| Space, Q/E | Run/pause strategic calendar and adjust speed |
| 1/2 | Enter RKM intercept or white-dwarf interdiction; campaign scoring requires the ship to have reached the selected industrial system |
| A/D (tidal mission) | Spend one limited radial tug burn toward/away from the black hole |
| R (encounter) | Restart the current mission with exactly the original system snapshot |
| 3/Tab (encounter) | Debrief and apply a single victory/defeat result to the strategic site |

The strategic screen uses the existing waypoint geometry. An additional Orion seed is the initial Blight industrial source. A cyan auxiliary dot denotes stored kugelblitz holes; gold rings denote Guf relays. Lines indicate packets currently moving at finite speeds. The right-hand system dossier reports industry, energy, materials, stored holes, array status and packet count.

## Simulation contract

- Strategic economics are fixed game units. Production accrues in deterministic 250-year substeps, with bounded inventories.
- Kugelblitz fabrication requires a sufficiently industrialized Drummer system and stored power/material. This is explicitly a fictionalized strategic cost proxy; actual creation and confinement of kugelblitz black holes is speculative.
- Weena can recruit industrial systems, with coordination taking distance in light-years in years plus a scenario-defined 2,000-year political/technical delay. Recruitment converts industrial ownership only on arrival.
- The Chamber of Guf requires three Weena-controlled, charged industrial nodes. Each node pays its firing cost. Separate beam packets launch with timed offsets and converge after the longest individual emitter-to-target light path. This is a playable approximation of Episode XIX's galactic attack cone and Episode XX's delayed responses.
- Blight packets propagate at 0.10c plus 1,500 years of resistance, with an additional global Orion causality check. A deterministic autonomous planner seeds nearby reachable targets every 2,500 years, with at most two invasion packets in flight. The player can pause autonomy with N.
- The darkened galaxy texture now incorporates local territorial control around industrial stars. Blight-held worlds darken within the Orion causal limit; Drummer-held enclaves become cyan and Weena-controlled sites become defended regions. The static historical front still provides the broader backdrop.
- Player time jumps backward rebuild the initial strategic state at the target year, preserving the autonomy setting. Forward jumps advance the ongoing war and its packet queue. The ship's coordinate/proper-time calculations remain independently handled by the navigation model.

## BG-STR-03: tactical encounter configuration

Selecting an industrial system now presents a live mission briefing along the bottom of the galactic atlas. A campaign encounter requires the player's 0.76c ship to be within 500 light-years of that selected named site. Missions entered elsewhere remain SI sandbox encounters. Completing a scored battle locks that system for **250 galactic years** of rearm time, preventing repeated instantaneous rewards.

At launch, `src/mission.rs` snapshots the selected system's owner, defenses, industry, stored kugelblitz count, and current strategic year into a stable `MissionProfile`. Later economic changes leave the active battle's starting conditions unchanged.

- **RKM interception.** Defensive strength scales the opponent flagship's shield joules and engagement distance. Industrial capacity controls the opponent's starting approach speed, cross-section, firing interval and a genuine squadron of up to four independently shielded escorts. Each escort travels through the existing SI kinematics and fires its own rounds using causal player-history observations. Player rounds use earliest swept-circle intersection across the flagship and all surviving escorts. A mission victory requires disabling the full enemy fleet. Local Drummer kugelblitz stock powers a larger friendly shield envelope. The enemy's causal targeting still reads its recorded worldline history. All starting samples are reset to the newly configured positions and velocities. The scene retains SI joules, metres, seconds and relativistic projectiles; shields remain fictional gameplay systems.
- **White-dwarf interdiction.** Site conditions select the stellar black-hole mass (Weena's site uses the episode XIX 12 M_sun remnant), approximate white-dwarf mass/radius, and initial separation. The stellar remnant is distinct from the abstract micro-kugelblitz inventory. After mission launch, each A or D key press applies a bounded 80 km/s radial tug impulse to the white dwarf, toward or away from the black hole. Drummer industrial Halo stock gives extra available tug burns. The unchanged two-body numerical solver determines the resulting orbit and disruption. Successful disruption within 300 SI seconds wins the mission; reaching the limit intact loses. The tug is a clearly identified external spacecraft intervention.
- **Debrief.** Destroying a hostile shield or achieving a tidal breakup damages local strategic defenses and may transfer ownership to the Drummers. A win at a friendly site reinforces its defenses and energy stores; defeat reduces defenses and material. Withdrawal has no strategic resource effect. Replaying with R retains the launch snapshot and cannot itself claim rewards.
- **Rearm control.** `StrategicState.last_tactical_year` enforces one scored engagement per star system per 250 galactic years, and the map reports whether a system is ready.

### Physics and design tests

New tests verify deterministic changes to shield strength, escort fleet size, enemy cadence and stellar remnant masses across distinct systems; matching initial SI worldline samples and multi-ship first-hit collisions; finite, causal projectile speeds; fixed mission profiles after industrial changes; limited radial tug budgets; and 250-year tactical cooldown. Existing SI physics tests remain in `src/si_bridge.rs`.

## History mapping

Episode XVII's hidden planetary accelerator establishes the Remnant's stored kugelblitz drives. Episode XVIII's O-star Dyson construction establishes the industrial background for hostile black-hole weapons. Episode XIX grounds Weena's supergiant/12-solar-mass remnant and multi-system Chamber of Guf. Episode XX grounds long light-time delays, counterfire from individual swarms, and the need for anticipation or evasive motion.

Site alignments, project prices, system ownership and action rates are balance assumptions. There is no assertion that the texts supply exact star coordinates.

## Next implementation milestones

1. Add actual defending squadrons, supply convoys, and orbital infrastructure to the single-opponent SI encounter instead of proxying local fleet strength through one defender's shields and cadence.
2. Build mission-triggered tactical setups for the Chamber of Guf and inter-system Halo battles, with packet interception and finite beam steering.
3. Introduce persistent campaign saves, deterministic replay, and scripted Episodes XVII through XXII scenarios.
4. Evolve the visual occupation mask from independent industrial bubbles into a connected front driven by autonomous sector captures.

## Verification

Unit tests cover scarce resource production, array node requirements, emitter-to-target light travel, Weena recruitment latency, Orion causal bounds, distinct SI missions, and tactical rearm. The dedicated GitHub workflow runs `cargo test --all-targets` and `cargo check --bin blighted-galaxy` on the PR.

### Actual executable visual checks

The existing BG-OVR workflow now captures additional genuine GPU frames of the four-escort Orion fleet and the Weena 12-solar-mass tidal scenario, using the fixed command-line QA fixtures `--demo-mission 9` and `--demo-mission 6 --demo-tidal`. These scenes are read-only previews, so they cannot change campaign ownership without the ship's travel and mission debrief. The screenshot artifacts share the existing 30-day GitHub retention policy.
