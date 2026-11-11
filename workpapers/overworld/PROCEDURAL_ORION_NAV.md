# BG-OVR-01: Procedural Galactic Navigation

**Implemented:** November 11, 2026  
**Runtime:** Macroquad 0.4.14; original Rust game remains intact  
**Source evidence:** [Blighted Galaxy Macro-history](../lore/macro-history.txt), [GDD](../../design/GDD.md) and [Patricia Tries](../../stories/manuscripts/patricia-tries.txt)

## Art direction and geometric constraints

The reference is the **face-on Milky Way chart provided by the author**, with a faint polar coordinate stencil, pale logarithmic arms, a warm barred bulge and a few deeply blackened regions speckled with distant crimson points. The Blight is an occlusion of visible starlight by densely replicated Jenkins Swarms. The red star-points represent scattered occupied systems or visible machine leakage. Their appearance can suggest blood spatter on first glance while resolving to small stars on closer inspection.

The game's map textures are created from scratch by code at startup. The background contains a 1152² RGBA field of exponential disk light, four logarithmic arms, a tilted bar, a bulge and deterministic stellar speckles. The occultation layer is a separate 1152² texture generated from the causal front, direction-dependent defensive delay and sparse red source positions. The layer is refreshed in 1,000-year chronology buckets and uploaded to the existing GPU texture.

The red tendrils **never extend beyond the global causal radius**. Their fingertips vary in reach due to thousand-year-to-tens-of-thousands-of-year delays. Red patches elsewhere in the galaxy are ordinary stellar color variation and are not claimed as separate Blight outbreaks. There is no full-disk blood-red recoloring. The old, heavily darkened region within the local front is opaque enough to suppress visible disk stars, with patchy boundaries and finite survivors.

The title imitates late-1980s and early-1990s PC science-fiction graphic design through a few layered metallic text passes, a dark red underline and a **sparse spray of individual red dots** under the letters. The map remains a functional simulation display with understated labels and astronomically grounded coordinate overlays.

## World coordinates

The 2D world uses light-years in the face-on Galactic plane:

```text
Galactic Center:  (0, 0)
Approximate Sol:  (0, -26,000 ly)
Orion Blight seed: (-1,100, -25,900 ly)

Galaxy texture: [-59,000, +59,000] ly in X and Y.
Galactic disk visible to approximately 52,000 ly from the center.
```

The physical map is only a top-down strategic projection. The original Milky Way asset's plotted labels are a useful design reference, while absolute azimuthal orientation and the locations of Olympus, Bel-Air or the Nebular Wastes remain **schematic** because the recovered manuscript provides relative routes and narrative regions rather than galactocentric coordinates.

The major radial stencils are centered on the Galactic Center:

```text
r = 15,000; 30,000; 45,000; 60,000 ly
galactic longitude = 0,30,...330 degrees
```

The local Orion-region stencil uses 15-degree angular sectors, 2,000-light-year rings and a separate dashed causal envelope. A click yields both addresses: galaxy sector and Orion sector, plus distances from the Galactic Center and outbreak. The selected sector has a current modeled Blight status.

## The hard expansion law

The global *potential* wavefront is

```text
R_causal(t) = 0.10 (ly/year) * t_years.
R_causal(100,000 yr) = 10,000 ly.
```

A local direction `theta` experiences a delayed arrival:

```text
R_local(theta,t)
  = clamp(R_causal(t) - 0.10 * Delay(theta,t),
          0, R_causal(t)).

Delay(theta,t) >= 0.
```

The modeled delays combine several broad resistance belts, narrow favorable expansion corridors and a gradually mobilizing post-human defense network. These delays are illustrative campaign parameters, not new historical facts. A late-era subset of previously occupied territory is tagged reclaimed after repopulation. The status machine never assigns Held, Contested or Reclaimed outside the global causal radius.

**Invariant:** at `t = 100,000 years`, every colored occupied or contested system is at radius <=10,000 light-years from the Orion seed. No tendril, frontline, node or procedural Blight stencil extends beyond that boundary, including in the texture generator. The local sector front remains inside the circle even when resistance delays differ by many millennia.

The hypothetical Blight present at 100,000 years is strictly a *regional* infection. Macro-history episodes depicting distant effects can be represented through knowledge, transmissions, migration, weapons and later history, without implying instantaneous uniform control of the disk.

## Macroscale historical sites

| Waypoint | Recovered narrative role | Placement |
|---|---|---|
| Sol / Orion | Blight's Orion Arm origin and surviving human home neighborhood | Approximate Orion Spur |
| Amicus | Group-mind triage and post-human defensive project | Orion region |
| Olympus | Titans' encircled infrastructure and Perseus last stand | Perseus schematic |
| Great Library | K-II archival Dyson civilization and transcendence crisis | Outer-arm side of Perseus schematic |
| Outer Arm exodus | Population flight and wild-black-hole Halo Drive recovery | Outer disk schematic |
| Bel-Air | Wicker lineage and Alfonso/Tatyana's ringworld | Beyond Outer Arm |
| Nebular Wastes | Formoi/Drummers, Great Aunt Weena and Chamber of Guf | Approx. 5,000 ly north of Bel-Air |
| Perseus Front | Dimming sky, Sauronoids and delayed warnings | Perseus narrative region |
| Canis Major Exodus | Flight off Galactic plane to a distant overdensity | Projection of an outbound corridor |

The historical manuscripts give years to specific episodes unevenly. **The T+100,000 game clock is an openly parameterized outbreak scenario**, and must not be mistaken for a dated canonical event list.

## Relativistic navigation and encounters

Map time is galactic coordinate time. A straight, constant-speed flight uses

```text
distance_ly = |destination - start|
galactic_elapsed_years = distance_ly / beta
proper_elapsed_years = galactic_elapsed_years sqrt(1-beta²).
```

At the default `beta=0.76`, a 5,000-light-year leg takes about **6,579 coordinate years** and **4,275 proper years**. The manuscripts describe another ~1/1000 subjective-time process; this is not ordinary special-relativistic dilation and would need an explicit separate cognition/suspension mechanic if implemented.

Flight progress updates the ship's true galactic position and the campaign clock; paused timeline steps do not move a ship. Seeking to a past/future epoch clears active travel to avoid inconsistent trajectories. A new arrival can emit an encounter notice tied to the nearest named source waypoint. Keys 1/2 enter the existing SI intercept or tidal encounter screens, and 3/Tab returns to the map.

## Controls and determinism

- **WASD/arrows**, right drag: pan; scroll, Z/X: zoom; C: Galactic Center; O: Orion; F: frontier focus
- **G** toggles the Galactic Center stencil, **L** the Orion 15-degree/2-kly sector stencil, **B** the occultation layer
- Left click: select physical destination; **Enter**: launch at 0.76c
- **Space**: pause, Q/E: years per real second, comma/period: 1-kyr time jump, Home/End: 0/100-kyr scenario anchor
- **H**: context help; **1** and **2**: existing combat encounters; **3/Tab**: return; **Esc**: quit

There is no external texture dependency, and random fields have fixed seeds. Unit tests cover the 10,000-light-year causal cap in every azimuth, directional delay and tendril range, stable dual-frame sector addressing, relativistic travel and proper time, timeline seek behavior, and texture coordinate registration.

## Reproduction and build

```bash
cargo test --lib overworld
cargo check --bin blighted-galaxy
cargo run --release

# Old startup behavior:
cargo run --release -- --encounter
```

The existing gameplay engine, white-dwarf SPH experiments, SI intercept and tidal encounter math are unaffected. The campaign overview is deliberately an **active navigation implementation**, with additional historical event graph and fleet logistics to be added in later changes.

## Next engineering steps

Expose a scenario serialization format storing current epoch, waypoint coordinates and delayed historical knowledge. Replace the simple angular resistance-delay scenario with a deterministic event-queue engine across star systems and light-cone communication. Add ship acceleration and deceleration worldlines with time-dependent gamma, along with subjective cognition time as a separate instrument, before using the 5,000-ly Bel-Air to Nebular-Wastes journey as a scripted story encounter.

An eventual time-slice simulator should store **actual wavefront state** and causal history, rather than recompute broad angular resistance from current epoch. That would support Patricia Trie fork divergence, learned strategies, liberation campaigns, memetic warfare and the long-term strategic recovery described in the source.
