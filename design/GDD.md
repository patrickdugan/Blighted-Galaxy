# Blighted Galaxy — Game Design Document

**Status:** Recovery draft v0.2 · 9 November 2026  
**Provenance:** Authored as a reconstruction from six surviving manuscripts/workpapers, and the author's account of the lost Rust prototype. The source texts are preserved unchanged in the archive. This document contains design choices for a new game in addition to ideas explicitly present in the originals.

## 1. High concept

*Blighted Galaxy* is a physics-led, transhuman space-warfare game played across incompatible scales of time and energy. Players survive and maneuver at relativistic speed, collide in morphing post-human bodies, manage an industrial war across stars, and influence a distributed civilization whose decisions travel at the speed of light. Its power fantasy is the capacity to redirect stellar remnants or wield a black-hole engine. Its dramatic center is whether people can preserve a recognizable self while acquiring that power.

The intended first release is a **2D simulation-driven campaign**, consistent with the source treatment's Blighted Galaxy 1. A later version may explore 3D star maps and flight without undermining the 2D strategy layer.

**Genre:** Astrophysical action sandbox × character action/RPG × interstellar strategy × narrative encounter.  
**Player perspective:** One persistent transhuman or post-human captain, with command over local infrastructure and allied units.  
**Primary platform:** PC, initially Windows/Linux desktop, keyboard-and-mouse, single-player.  
**Current prototype:** A clean-room Rust/Macroquad sandbox, with simplified gravitational interactions and a 2D duel. This is a foundational test bed, not evidence that the lost engine source was recovered.

### Core sources

- [Game Design Treatment](../workpapers/design/game-design-treatment.txt): three gameplay modes, power tiers, tactical parameters, meme combat, war-economy resources, initial physics prototype.
- [Macro-history](../workpapers/lore/macro-history.txt): 22 outlined episodes and their galactic chronology.
- [You Never Go Full Blight](../stories/manuscripts/you-never-go-full-blight.txt): evacuation, nanite attack, Weena's schism.
- [Patricia Tries](../stories/manuscripts/patricia-tries.txt): the twins, first contact, relativistic duel and Hawking bomb interception.
- [Blowing Up The Galaxy](../stories/manuscripts/blowing-up-the-galaxy.txt): Chamber of Guf, strategic gravity, post-human political settlement.
- [White Dwarf Super Shotgun](../stories/manuscripts/white-dwarf-super-shotgun.txt): Sauronoids, defensive beam networks, stellar-remnant siege.

## 2. Nonnegotiable play principles

**Causality is terrain.** No faster-than-light travel or communication. Travel-time, sensor delay, and light-cone separation affect maneuvers, campaigns, and intelligence. The player can see events whose consequences were decided centuries ago.

**Momentum is a commitment.** At high relative velocity, beam leads, turning arcs, gravitational assists and deliberate deceleration matter. Close combat becomes possible after velocity matching. A high-speed encounter is tactically distinct from a slow duel.

**Mass and wattage define power.** Star classes, matter inventory, propulsion expenditure, heat rejection, available computation, and structural resilience form the power budget. All offensive effects trace back to a scarce physical resource, even when exaggerated for play.

**Identity is mutable.** Characters may refactor their nanoborg bodies, change clades, fork, regenerate, merge, or ossify into dangerous beliefs. Social and psychological vulnerabilities provide a second axis of conflict.

**The enemy learns over historical time.** The Blight is an evolving network of replicated Patricia Tries. Contact, violence, reconciliation, and propagation delay alter regional behaviors without a unified instantaneous omniscience.

## 3. Play structure

### 3.1 Tactical mission loop

1. **Observe:** Read delayed sensor feeds, inferred trajectories and threat likelihoods.
2. **Commit:** Choose approach vector, thrust budget, desired speed and reserve capacity.
3. **Intercept:** Lead light-speed weapons or slower projectiles; exploit massive bodies and incoming swarms.
4. **Resolve:** Match velocity to close, transform into combat form, box/strike/grapple/regenerate or attempt memetic contact.
5. **Recover:** Salvage useful mass, refuel, rescue, record the encounter in local history, propagate intelligence into the campaign map.

The source treatment frames this as relativistic jousting, transhuman duels and interstellar war economy. Those three modes should share persistent mass/energy state.

### 3.2 Sector loop

Survey candidate stars and remnants → choose production hubs → allocate star-derived power and nanite mass → develop engines/structures → fortify or strike → advance the campaign clock long enough for information and fleets to arrive → respond to the altered map.

Stellar systems can occupy different strategic states because the relevant news has yet to reach them. Orders are issued under uncertainty. A lost battle can keep producing consequences for centuries.

### 3.3 Time scales and controls

| Scale | Typical experience | Player control |
|---|---|---|
| Nanoborg contact | Seconds to minutes | Real-time boxing and transformation |
| Local system duel | Minutes to days | Real-time or compressed |
| Interstellar intercept | Years to millennia | Program a trajectory, then accelerate the clock |
| Sector campaign | Centuries to tens of thousands of years | Pausable strategy, logistics, delayed orders |
| Galactic Patricia Trie | Millennia of backscatter | Influence and infer, rather than instant global commands |

Long-duration simulation should have deterministic, low-detail analytic propagation between high-detail encounter windows. Time acceleration must preserve causal order.

## 4. Player, encounters and combat

### 4.1 Player growth

The player's body begins as a **Tier 2 outer-arm bionanoborg** in the source's taxonomy. A progression arc can approach **Tier 3 baseline post-human** architecture and eventually black-hole-powered Tier 5 traits. Retain political and personal consequences, and treat radical transitions as campaign choices.

Proposed attributes:

| Attribute | Narrative basis | Game expression |
|---|---|---|
| Structure / viable mass | Nanoborg tissue, regenerating organs | Hit points, limb and hull damage |
| Energy throughput (W) | Star, fusion and Hawking power | DPS, repair rate, thrust and cooling |
| Control bandwidth | Sensor/data processing and computation | Aim precision, multi-limb control, response delay |
| Coherence | Ossification/perversion/fragmentation | Memetic defense, transformation risk |
| Matter inventory | Carbon/nanite stock | Manufacture, armor, swarm size |
| Identity & memory | Forks, backups and continuity of self | Recovery cost, dialogue, allegiance |

### 4.2 Combat model

At long range, high-speed projectile trajectories and beam interception govern outcomes. Every shot has a muzzle position, velocity, emission time, detection delay and a collision solution. Continuous collision is required because discrete-frame penetration fails at high speed.

Near contact, change the rhythm: local gravity and active thrust change approach geometry; boxing/grappling hitboxes, defensive shapes and regeneration decide the exchange. Fast, cheap broad-phase tests filter contacts, with swept volumes or conservative advancement during short high-speed contacts. Animation and particle rendering are driven by actual collision events.

At sufficiently close range, configurable tissue forms can attack with nanite erosion, beam emitters, high-energy strikes or target-specific sabotage. The original treatment explicitly suggests energy-limited regeneration and choosing morphology against opponents.

**Memetic combat** is a separate mechanic driven by trust, experience, a clade's pathologies and communication access. It can delay, mislead, persuade or alter another actor. It should have clear player-readable causal feedback.

### 4.3 Proposed playable roles

- **Alfonso:** experimental high-mobility nanoborg fighter and field operator.
- **Tatyana:** first-contact specialist whose research can change network behavior.
- **Weena:** post-human bomber fleet progenitor, then planetary-scale weapon architect.
- **A Blight node:** a responsive local system with a partial, delayed Patricia Trie.
- **Sauronoid defenders:** high-mobility, energy-constrained observers in later-war scenarios.

These are campaign/source-derived roles, not a finalized character roster.

## 5. Astrophysics and game feel

### 5.1 Physical contracts

The **SI physics module** is the reference layer. It uses kg, m, s, J, W, and K with documented equations and benchmark assertions. Gameplay scale conversion happens in an explicit adapter. The existing normalized 2D scene remains an arcade sandbox until the adapter and visual calibration are validated.

- **N-body gravity:** pairwise Newtonian gravity for bodies outside compact-object relativistic regimes. Symmetric forces, momentum invariant, velocity-Verlet time stepping, bounded close-approach substeps.
- **Relativistic kinematics:** Lorentz gamma and proper-time ratios for constant-velocity legs, kinetic energy at relativistic speed, and 1D Einstein velocity composition where applicable. An accelerated ship needs a proper worldline integrator later.
- **Event horizons and evaporation:** Schwarzschild radius and Hawking temperature, nominal emission power and \(M^3\)-lifetime approximation for an isolated, uncharged, nonrotating hole. A real particle spectrum and final endpoint are open scientific questions. See the [SI physics workpaper](../workpapers/physics/SI_PHYSICS.md).
- **Tidal shearing:** radial differential gravitational acceleration against the object's self-gravity, displayed as a diagnostic and optional mission event. True deformation requires a material model.
- **Swept collision:** continuous paths for RKMs, body encounters and projectiles, subject to a configurable time-resolution cap.
- **Swarm scatter:** deterministic seeded particles initially; spatial hashing/compute-particle GPU rendering after profiling.

### 5.2 Visual grammar

Energy-rich encounters should reveal their *cause*: lensing around a compact object, Doppler-colored emissions, accumulating orbital distortions, jet or beam momentum, particle fronts and tidal debris. The first slice uses an inexpensive 2D sprite/particle pipeline. Dedicated shaders, scattering and camera exposure follow profiling.

Use high-contrast legibility at extreme zoom: a single source body can change from a small marker to the entire screen. Particle budget adapts to active effects and zoom. Light-travel delays need clear visualization on HUD and mission timeline.

### 5.3 Exaggeration policy

The fiction includes white-dwarf shotguns, industrial kugelblitz generators, star-scale mirror guns and singularity-powered bodies. Treat them as **fictional engineering hypotheses**. The game should document where it models textbook equations and where its missions assume breakthroughs, extraordinary infrastructure or altered outcomes. Avoid describing the current toy evaporation flash as a validated model of a real Hawking endpoint.

## 6. Stellar war economy

Source-defined classes: stars, supergiants, carbon/generic nanite mass, white dwarfs, neutron stars, Hawking-class black holes, massive/rotating black holes, population and wattage.

Proposed object operations:

- Harvest local matter into factories, habitats and swarm fabric.
- Enclose stars to acquire industrial power and computing budget.
- Manufacture propulsion and weapons using captured mass and energy.
- Corral stellar remnants, change trajectories and employ gravity assists.
- Build and defend regional communications and interception arrays.
- Replicate combatant minds when the story allows, with explicit material and coherence costs.

A technology tree should unlock controlled physical processes and forms of matter handling rather than a generic list of arbitrary damage upgrades.

## 7. The Patricia Trie world simulation

Each star has a node with its own local history. Events generate signed records. Records propagate along light-speed links with latency determined by distance. A node updates its policies only when a message actually arrives. Distinct branches can remain inconsistent over millennia.

Proposed minimal simulation state:

- node_id, time, location, energy budget, available manufacturing mass;
- local trie root and branch provenance;
- records of contacts, attacks and cooperative interactions;
- policy weights for defense/assimilation/coexistence;
- outgoing message queue with physical delivery time.

This structure makes the episode XXII branch conflict playable as information warfare, while preserving the stories' uncertain psychology of the Blight. The distributed consensus is a *narrative rule* until a concrete protocol and benchmark implementation exist.

## 8. Scenario spine and timeline

The surviving macro-history identifies 22 planned episodes. Some early episodes supply titles only; do not infer additional events for those slots.

| Episode | Source title | Gameplay opportunity / known story |
|---|---|---|
| I | Blockchain All The Things | Origin-era background, title only |
| II | The Two Hundredth Amendment | Constitutional/strategic backstory, title only |
| III | Revenge of The Blight | Background, title only |
| IV | Weapons Of Mass Destruction | Background, title only |
| V | The Federation Strikes Back | Background, title only |
| VI | Return of The Blight | Background, title only |
| VII | The Post-Human Awakens | Background, title only |
| VIII | The Last Human | Amicus rescues last baseline humans |
| IX | The Harrowing of Amicus | Ancient group-mind challenges the Blight |
| X | A New Ambition | Emigration toward the Canis Major region |
| XI | Dim-sky Rising | Perseus Arm observes the dimming sky |
| XII | The Fall of Olympus | Titans prepare defenses |
| XIII | Perseus' Last Stand | Titans fuse under encirclement |
| XIV | Roast of The Great Library | Library civilization, war and transcendence |
| XV | Gotta Gotterdammerung | Relativistic destroyer fleet campaign |
| XVI | Interstellar Underground | Remnant evacuation and wild-hole propulsion |
| XVII | You Never Go Full Blight | Wicker's evacuation, nanite attack, Weena's defection |
| XVIII | Patricia Tries | Alfonso/Tatyana first contact and bomber interception |
| XIX | Blowing Up The Galaxy | Chamber of Guf and reconciliation |
| XX | Death Star Hop-Scotch | Coordinated stellar movement against beam attacks |
| XXI | Double-barreled White Dwarf Shotguns | Stellar-remnant siege |
| XXII | Amicus Brief | Light-speed conflict among diverged trie branches |

**First campaign slice:** Episodes XVII–XIX provide the strongest continuous sequence with surviving full prose. Episode XXI is a separate technical/visual stress test and later-campaign chapter.

## 9. Vertical slice definition

The first publicly playable benchmark mission, **The Cluster Bomber Intercept**, borrows its dramatic purpose from *Patricia Tries*, with new encounter scripting:

- Player enters a star system at high relative velocity, with a hostile bomber closing on a predicted target.
- The UI displays the light-time lag and intercept window.
- Player can match trajectory, fire a mass driver and prevent detonation by adding mass to the bomber's black-hole payload (an encounter-specific fictional control system).
- If close enough, transition into transhuman contact boxing against a morphing opponent.
- The explosion path uses an explicitly calibrated energy proxy; visible stellar shearing is triggered by a configurable toy-damage threshold until a numerical stellar model exists.
- The saved scenario records mass and energy totals, player choices, elapsed coordinate time and subjective time.

**Acceptance tests:** consistent fixed-step simulation independent of render FPS, deterministic replay from a seed, stable two-body orbit within tolerance, conservation tests where appropriate, high-speed tunneling regression, correct light-delay values, and a performance budget on the target RTX 5080.

## 10. Production sequence

**M0 recovered corpus:** preserve the source transcriptions, build searchable story/workpaper indices and GDD. Complete upon merge of this change.

**M1 SI physics bench:** calibrated constants, N-body acceleration and symplectic evolution, relativistic energy/time, nominal Hawking temperature/power/lifetime, diagnostics and executable benchmarks. Included in the same PR; still independent of the arcade renderer.

**BG-PHY-06 photon camera and cold-star bridge:** non-equatorial Kerr null geodesics with a finite-radius observer and Carter constant, cold relativistic electron EOS, a spherical hydrostatic reference profile, and stratified 3D SPH sampling. See [BG-PHY-06](../workpapers/physics/KERR_CAMERA_WHITE_DWARF.md). The live joust's optical history and camera ray tracer remain separate integrations.

**M2 calibrated encounter:** choose world-unit mapping and camera model, integrate SI time/space into a standalone or alternate playable scene, add trajectories and retarded optical states.

**M3 boxing and morphogenesis:** continuous collision for deformable entities, form state, energy-limited regeneration and spatial partitioning.

**M4 particle rendering and astrophysical set pieces:** GPU particles, tidal debris, beams, black-hole accretion/glow and measurable effect budgets.

**M5 campaign:** systems, 22-episode branching scaffolding, delayed order queues and localized Patricia Trie propagation.

### BG-PHY-07 integration gate

The Kerr camera now supports an **equatorial moving-source optical acquisition**: an observer receives light propagated along a stationary Kerr null geodesic, and the ray intersects a recorded circular-orbit combatant history. The 2D source-envelope collision remains chart-coordinate approximated, and the default solar-system joust continues using flat-spacetime causality until an observer/worldtube representation with consistent Kerr coordinates is ready.

The three-dimensional cold electron EOS now has a **damped isolated-SPH relaxation** experiment. The numerical residual force drops after extended damping, while the low-resolution particle star shrinks below its hydrostatic radius. This exposes a physical-fidelity blocker for any white-dwarf shotgun or near-hole collision level. Advance to adaptive smoothing and higher-resolution, converged stellar initial conditions before turning this into campaign damage balancing.

See [BG-PHY-07](../workpapers/physics/KERR_CAUSAL_RELAX.md).

### BG-PHY-08 fidelity and perception milestone

The encounter optics now permit an **orbiting observer** in a Kerr background. The camera accounts for its local velocity relative to the ZAMO tetrad, rotates the photon geodesic to the observer's reception azimuth, and intersects the ray with a recorded source history. The source can undergo a smooth noncircular perturbation. The target envelope is measured using a local rest-frame metric estimate, retaining finite-chart and two-dimensional constraints. The main solar-system duel still uses the established flat-spacetime causal engine.

The white-dwarf branch adds **adaptive 3D SPH** with per-particle smoothing lengths and a grad-h pressure correction. Fixed-h and adaptive-h models remain independently selectable, and their measured radius, force residual, density-profile mismatch and undamped free-evolution stability inform further engineering. The first comparison improves the predicted radius while retaining a substantial discrepancy from the hydrostatic reference. The white-dwarf equilibrium issue stays open until a resolution-converged radius within 5% and a stable multi-dynamical-time undamped hold are measured.

Refer to [BG-PHY-08 workpaper](../workpapers/physics/ADAPTIVE_SPH_MOVING_OBSERVER.md) and the open convergence issue.

### BG-PHY-09 gravity consistency and star sampling

The scientific physics bench now has a separate adaptive gravitational softening option. It differentiates the same h-dependent Plummer pair potential used for energy bookkeeping, adding a gravitational grad-h correction to the original variable-h pressure method. A deterministic shell-balanced equal-mass initializer supplies alternative 3D white-dwarf initial conditions from the recovered hydrostatic reference.

The radius and density benchmarks support continued development of this numerical branch. The matched potential improves short-duration radii at small N, while the undamped post-relaxation star expands substantially and the 512-particle radius remains outside the ±5% acceptance gate. The primary astrophysical story encounters therefore continue to use a clearly labeled fluid proxy until the resolution, stability and energy-momentum gates converge.

Full equations and measured results: [BG-PHY-09](../workpapers/physics/GRAD_H_GRAVITY_SHELL_SEED.md).

### BG-PHY-10: kernel-matched self-gravity and stable stellar debris

The physics research engine implements a Poisson-matched cubic-spline gravitational potential, including the h-dependent force correction, and a new approximately cubic equal-mass shell layout. The M4 gravity is exactly Newtonian beyond its compact support. The numerical benchmark reaches 1,024 particles with a cold-degenerate reference white dwarf.

After 0.6 simulated seconds, the 1,024-particle M4 model falls within 5% of the Newtonian hydrostatic reference radius, with improved logarithmic density-profile matching. This is a short damping interval, roughly 0.23 stellar dynamical times. Lower-resolution undamped test stars continue to contract substantially, so the source-derived stellar weapons remain calibrated **fictional effects** in the playable engine pending a stable multi-dynamical-time stellar model.

Workpaper: [BG-PHY-10](../workpapers/physics/M4_GRAVITY_CUBIC_SHELLS.md). Physical acceptance issue: [#8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

### BG-PHY-11: Computational scale and stellar equilibrium validation

The cold-degenerate white-dwarf branch now uses a spatial cell index for compact-support adaptive SPH density and gravitational softening corrections. It retains exact long-range pairwise gravity for controlled numerical comparison and reproducible summation. Additional physical diagnostics measure 50%, 90%, 95% and 99% enclosed-mass radii and mismatch from the continuum hydrostatic mass profile.

The dedicated high-resolution benchmark evaluates a 1,024-particle star after 3 simulated seconds of damping and approximately three subsequent stellar dynamical times without damping. Radius, mass CDF, energy and momentum are observed throughout the hold, and failed *physical* acceptance gates are recorded alongside passing software numerical invariants. A stable initial star remains a prerequisite for promoting the fictional white-dwarf shotgun and tidal-breakup effects to a calibrated physical simulation.

See [BG-PHY-11](../workpapers/physics/NEIGHBOR_MASS_PROFILE_HOLD.md) and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

### BG-PHY-12: Why the white dwarf oscillates

The stellar physics engine now measures the cold-degenerate star's *pressure* and *self-gravity* contributions separately in radial mass shells. It compares the SPH force budget against the independent Newtonian hydrostatic profile, including net radial and tangential force fluctuations. A local pressure-to-gravity force-virial ratio helps diagnose large-scale imbalance.

The experimental preparation stage offers a time-dependent damping schedule with a strong early phase and gradually weakening damping before a fully undamped hold. A matched-duration A/B benchmark records the star's radius, 99%-mass radius, radial mass-CDF mismatch, force budget and energy/momentum invariants. The finite-volume/timestep tests determine whether a schedule improves equilibrium. The story-derived stellar-shearing mechanics remain explicitly provisional until the hydrostatic acceptance gates pass.

See [BG-PHY-12](../workpapers/physics/RADIAL_FORCE_RELAX_SCHEDULE.md) and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

### BG-PHY-13: SPH smoothing scale and cold white-dwarf convergence

The stellar research engine now has a controlled parameter-sweep runner for the dimensionless SPH smoothing length eta. For each resolution it reuses the same nearly cubic particle realization, M4 Poisson self-gravity, cold electron EOS and physical times, then measures radial acceleration against the independent hydrostatic enclosed-mass field. Explicit kernel-support neighbor counts distinguish physical resolution from spatial-index candidates.

Each variant undergoes a three-dynamical-time undamped hold. Radius, enclosed-mass percentiles, radial mass-CDF mismatch, force-virial balance, and energy/momentum invariants feed a separate physical-equilibrium verdict. The benchmark explores N=256, 512 and 1,024 with reproducible CSV artifacts. The next graphics and gameplay stages retain their existing stellar destruction proxy until the [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) requirements are verified.

See [BG-PHY-13](../workpapers/physics/ETA_HYDROSTATIC_SWEEP.md).

## 11. Unresolved design decisions

- **Fidelity:** whether orbital maneuvers use realistic dynamical constants in all encounters, or curated time/length scales for playable pacing.
- **Particle scale:** maximum independent entities and whether swarms are particles, fields or clustered agents.
- **Narrative control:** single persistent protagonist or campaign switching between descendants and clades.
- **Scientific assumptions:** mechanism of industrial tiny-black-hole manufacture, white-dwarf handling at near-light-speed, radiation shielding, viable mental continuity and the extent of memetic influence.
- **Art and distribution:** no visual art direction, rights-cleared soundtrack or public code license has yet been approved.

The full fiction and historical workpapers remain primary sources; engineering changes in this GDD are proposals pending further testing.
