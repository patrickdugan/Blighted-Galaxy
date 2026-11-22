# BG-STR-06 & M3 | Episodes XVII–XXII, Causal Patricia Trie, Convoy/Halo Encounters, and Transhuman Contact Boxing

**Scope:** Playable campaign expansion (`src/patricia_trie.rs`, `src/episodes.rs`, `src/strategy.rs`, `src/mission.rs`, `src/campaign_save.rs` v3) and M3 transhuman close-contact boxing & morphogenesis (`src/morph_boxing.rs`, `Mode::Boxing` in `src/main.rs`).

## 1. Causal Patricia Trie World Simulation (`src/patricia_trie.rs`)

Implements GDD §7 with zero external cryptographic crates:

- **Per-system state (`PatriciaNode`):**
  - `records: Vec<TrieRecord>` — signed, content-addressed records (`record_id`, `signature`) derived deterministically with FNV-1a + SplitMix64.
  - `root_hash: u64` — deterministic, insertion-order-independent 64-bit Merkle-Patricia root computed over a canonical 16-way nibble-prefix compressed radix trie (`RadixNode`).
  - `policy: PolicyWeights { assimilate, defend, coexist }` — normalized simplex weights derived strictly from local records (`sum == 1.0`, each strictly positive and bounded).
  - `outgoing: Vec<InFlightTrieMessage>` — light-speed delivery queue (`c = 1 ly/yr`, `arrival_year = launch_year + distance_ly`).
- **Causal locality:** Remote nodes never mutate their trie root or policy weights until `arrival_year <= year`.
- **Memetic broadcast & Episode XXII convergence (`Shift+M`):** `broadcast_memetic_record` emits a signed `MemeticBroadcast` record together with the origin's branch history at `c = 1 ly/yr`. Upon causal arrival, each peer newly accepting a reconciliation broadcast echoes its unseen local branch records back out at `c = 1 ly/yr`. Once light cones intersect, all active systems hold the same content-addressed set and converge to `branch_divergence_count() == 1`.

## 2. Episodes XVII–XXII Campaign Spine & Weena Planner (`src/episodes.rs`)

Tracks the six continuous late-war chapters from `design/GDD.md §8` and the recovered `workpapers/lore/macro-history.txt`. Recovered narrative summaries are tagged `[CANON]`; playable progression predicates are tagged `[STAGING]`:

| Episode | Title | Playable `[STAGING]` Objective Predicate |
| --- | --- | --- |
| XVII | You Never Go Full Blight | Any Drummer system reaches `industry >= 2.5` and holds `kugelblitz >= 1` (or launches a Halo strike) |
| XVIII | Patricia Tries | At least 1 tactical or boxing victory AND `>= 1` contact or memetic record in `PatriciaTrieWorld` |
| XIX | Blowing Up The Galaxy | Defend a Drummer system with a Guf screen after an array discharge, OR capture Nebular Wastes (site 6) |
| XX | Death Star Hop-Scotch | Guf defense / Nebular Wastes secured AND a Drummer convoy delivered/escorted (or `>= 6` Drummer systems) |
| XXI | Double-barreled White Dwarf Shotguns | Win a White-Dwarf Tidal Interdiction OR a tactical Halo Strike |
| XXII | Amicus Brief | Reconcile `PatriciaTrieWorld` to `branch_divergence_count() <= 1` with `average_coexist() >= 0.45` |

### Great Aunt Weena Autonomous Planner (`Shift+N`)

`plan_weena_expansion` runs on exact 2,500-year planning boundaries when `weena_autonomy` is enabled (toggled independently of Blight autonomy `N` via `Shift+N`):
1. Expands industry at Weena worlds with `industry < 2.0` when stores permit.
2. Constructs synchronized Guf relay nodes (`Project::GufNode`) at qualified Weena worlds.
3. Causally recruits reachable Independent/Drummer systems (`claim_for_weena`, light travel + 2,000 yr delay).
4. Issues finite-slew `aim_guf` commands and fires `fire_guf` only once `>= 3` mounts are aligned and `array_charged >= 100.0`.

## 3. Convoy Escort (`Key5`) & Tactical Halo Strike (`Key6`)

- **`TacticalKind::ConvoyEscort` (`5` on overworld):** Launchable when an active Drummer `CargoConvoy` is tied to the selected system or within 2,500 ly of the player's ship. Spawns a fast bandit raider squadron (`beta = 0.028..0.055`, 2–4 nimble escorts, Flak unlocked). Victory marks `convoy.escorted = true` (absorbing the next Blight interdiction cleanly) and awards local escort credit; defeat destroys the escrowed convoy; withdrawal leaves convoy and 250-yr cooldown untouched.
- **`TacticalKind::HaloStrike` (`6` on overworld):** Launchable at a hostile Blight/Weena system when Drummers hold a stored kugelblitz. Snapshots a high-velocity `0.76c` approach (`intercept_range_m = 1,400 Mm`) with `WeaponKind::Siege` armed and causal worldlines validated. Victory consumes 1 nearest Drummer kugelblitz upon resolution and inflicts heavy `4.5` defense + `0.75` industry strike damage; withdrawal retains the kugelblitz with zero double-charging; defeat applies normal defensive loss.

## 4. M3 Transhuman Close-Contact Boxing & Morphogenesis (`src/morph_boxing.rs`)

- **Morphologies (`MorphForm`, keys `1–5` in `Mode::Boxing`):**
  - `Streamlined`: high thrust/speed, minimal hitbox (`2.8 m`).
  - `AegisGuard`: wide `5.4 m` hexagonal barrier (`65%` damage mitigation), lower speed.
  - `StrikerLimb`: extended `19.5 m` continuous-swept strike reach (`142 kg` base impact).
  - `SwarmErode`: `22 m` nanite cloud eroding enemy viable mass and matter while bypassing half of armor mitigation.
  - `SingularityCore`: `440 MW` surge wattage and piercing strikes at sustained heat and coherence cost.
- **Broad-phase & zero-tunneling narrow-phase:** `UniformGrid2D` (`BTreeMap` cells for cross-platform determinism) indexes swept body/limb AABBs; narrow-phase resolves exact times of impact in `[0, 1]` via `crate::swept_circle_hit` and `swept_aabb_hit` (Minkowski slab intersection).
- **Physical resource budgets:** Viable mass (`mass_kg`), energy buffer & wattage (`energy_j`, `wattage_w`), nanite matter (`matter_kg`), cognitive coherence (`0..1`), and radiative heat (`0..1`), with energy- and matter-limited structural regeneration.
- **Velocity-matching entry & memetic debrief:** `can_enter_contact_boxing(&SiJoust)` gates in-battle transition (`<= 80 Mm` separation and `<= 0.015c` relative speed, triggered via `7` in `Mode::Intercept`); direct launch from overworld via `7`. Pressing `F` (or `E`) fires a memetic disruption pulse that degrades opponent coherence and feeds a signed `TrieEventKind::ContactDuel` record into `PatriciaTrieWorld` on `Tab` debrief.

## 5. Save Schema v3 & v2 Migration (`src/campaign_save.rs`)

- `SAVE_VERSION = 3` persists `weena_autonomy`, `episodes`, `patricia_trie`, and `CargoConvoy.escorted`.
- Valid schema `v2` saves load and migrate cleanly (defaulting `weena_autonomy = false`, fresh `EpisodeTracker`, fresh `PatriciaTrieWorld`, and `escorted = false`).
- Invalid versions (`v1`, `v99`), tampered trie root hashes, forged signatures, invalid policy weight sums, superluminal trie messages, and out-of-bounds episode indices are rejected before mutating live state.

## 6. Verification & Read-Only QA Fixtures

```sh
cargo test --all-targets
cargo check --all-targets
cargo run --release -- --demo-mission 3 --demo-convoy-escort --snapshot-overworld bg-mission-escort.png
cargo run --release -- --demo-mission 6 --demo-halo-strike --snapshot-overworld bg-mission-halo.png
cargo run --release -- --demo-boxing --snapshot-overworld bg-m3-boxing.png
```
