# BG-PHY-07: Moving-source Kerr optics and damped white-dwarf SPH relaxation

**Status:** Executable research milestone with tests and CI artifacts.  
**Scope:** New Rust implementations in `src/kerr_optics.rs`, `src/stellar_relax.rs` and the supporting equatorial Kerr path sampler in `src/kerr_camera.rs`.  
**Previous work:** [BG-PHY-06](KERR_CAMERA_WHITE_DWARF.md).  
**Story basis:** The surviving *Patricia Tries*, *Blowing Up The Galaxy* and *White Dwarf Super Shotgun* manuscripts motivate the need for causally delayed perception and gravity-driven encounters. The numerical models here are new research, with explicit physical assumptions.

## 1. Kerr optical information from a moving combatant

An observer receives photons emitted at an earlier event on a source's history. Near a black hole, the flight path is a curved null geodesic and its coordinate travel time differs from the line-of-sight distance divided by c. Determining a null geodesic that joins an emission event on a moving orbit to an observer is a numerical emitter–observer boundary problem.

**Primary references:** [Propagation time delay and frame dragging in Kerr spacetime](https://academic.oup.com/mnras/article/516/2/1768/6673442) and [Radiation Transport Around Kerr Black Holes](https://dspace.mit.edu/server/api/core/bitstreams/0f2653e8-c7b1-4058-a672-6913364176c3/content).

The demonstration assumes an isolated 10-solar-mass rotating black hole described by the stationary Kerr metric, with a zero-angular-momentum observer (ZAMO) at radius **85 GM/c²** on the equatorial plane. A combatant follows a circular timelike equatorial Kerr test geodesic at radius **20 GM/c²**, with angular speed computed from the metric's circular-orbit formula.

The ray tracer starts at the observer, launches an incoming **past-directed** null ray from the ZAMO local tetrad, then integrates its equatorial Kerr path toward the black hole. Every sample includes a Boyer–Lindquist coordinate-time offset, physical chart x/y coordinates, orbital radius, and azimuth. The samples are dense enough to resolve a finite source envelope, with radial chord spacing bounded at approximately 0.65 gravitational radii.

The emitter's past is a bounded record of position and velocity in shared coordinate time. The algorithm intersects each photon chord with the matching section of the stored source worldline, using the existing swept-circle collision primitive. The earliest returning path that meets the source envelope supplies its emission time, observed position and light-travel delay.

This is **curved-spacetime causal perception**, at the level of an equatorial null path and a source orbit history. In particular:

- Photon propagation is a Kerr null geodesic, and the source trajectory is a Kerr circular test-body solution.
- The finite target envelope and segment collision use a Euclidean x/y mapping of Boyer–Lindquist chart coordinates. A physical local collision volume must eventually use a local orthonormal spatial tetrad and a four-dimensional timelike worldtube.
- The search scans a finite bundle of incoming screen rays, so it can miss narrow high-order lensed images. It returns the shortest-delay match from the sampled paths.
- The observer is held at a fixed radius. General observer motion, energy/redshift of individual photons, plasma scattering and full non-equatorial moving targets remain outside this integration.
- The new optical game is an **incoming sensor/acquisition** exercise. It does not pretend that the same past-directed ray solves an outgoing laser interception problem.

### Deterministic moving-source scene

The source orbit is phased so that a chosen null ray reaches it at the required **past** orbital angle. This makes a repeatable test of the optical solver without granting the simulated sensor instantaneous information. The source's present position is calculated separately for debug display.

The demo samples a ray bundle around the central optical axis, reports the connecting ray's impact parameter and emission time, and visually separates the received light from the target's present simulated position.

**Example measurements from CI on a 10-solar-mass black hole:** 

| Signed spin a* | Kerr coordinate delay | Radial Euclidean flight baseline | Source displacement during delay |
|---:|---:|---:|---:|
| 0.0 | 3.297829 ms | 3.201671 ms | 14.624 gravitational radii |
| +0.8 | 3.297196 ms | 3.201673 ms | 14.498 gravitational radii |
| −0.8 | 3.297178 ms | 3.201673 ms | 14.747 gravitational radii |

These are **chosen benchmark geodesics** with a source envelope of one gravitational radius, not universal delays for all observers, spin orientations or orbital radii. The approximately 96-microsecond excess in the nonrotating case is consistent with a gravitational propagation delay for this ray. A high-accuracy external ray-tracing comparison and convergence with respect to sampling remain necessary before treating the optical link as a precision Kerr timing solver.

### Controls and reproduction

Run `cargo run --release --bin kerr_optical_lab`. Keys **1/2/3** select spin 0, +0.75 or −0.75, **F** toggles automatic optical acquisition, **A/D** steer the manual photon probe, and **Q/E** adjust the source clock. **Z/X** or the mouse wheel controls zoom, and **R** resets. The cyan path is the curved incoming ray; a yellow marker shows light received from the earlier source position. The red marker shows current simulated position for debugging.

The headless `optical_relax_bench` tests spin 0 and ±0.8 with recorded histories, checks the null constraints and tests that an incomplete history cannot supply an unrecorded emission event. The original flat-spacetime causal combat AI remains available in the default game. Its full upgrade to Kerr geometry will require source/observer worldtubes, a ray shooting solver and a physically specified firing/propulsion law.

## 2. Damped 3D white-dwarf SPH relaxation

The BG-PHY-06 hydrostatic solver produces an ideal cold electron-degenerate star with the desired mass, radius and radial density profile. Sampling that one-dimensional continuum as 3D SPH particles introduces discrete-density noise, imperfect pressure gradients and transient velocities.

**Primary references:** [Stellar GADGET: SPH initial conditions and relaxation](https://academic.oup.com/mnras/article/424/3/2222/977734) and [SPH simulations of white-dwarf collisions](https://academic.oup.com/mnras/article/406/4/2749/1022010). Both studies emphasize relaxation of isolated stars before binary encounters, using much larger particle counts than this interactive prototype.

The new relaxation routine evolves an isolated collection of 3D equal-mass particles with the existing symmetric pairwise self-gravity and **cold relativistic electron EOS** pressure. It temporarily applies artificial damping relative to the mass-weighted bulk velocity:

    dv_internal/dt = a_SPH - v_internal / tau_damp

In the discrete implementation, velocity Verlet performs the conservative force step, then internal velocities are scaled by exp(−dt/tau). Positions and velocities are recentered against the original bulk frame to remove floating-point center-of-mass drift. Because damping intentionally removes kinetic energy, **energy is not conserved during relaxation**. The conservative subsequent encounter code is unchanged.

The diagnostic report records root-mean-square residual acceleration, logarithmic SPH density mismatch against the one-dimensional hydrostatic density profile, maximum particle radius, residual internal velocities, mass balance and center-of-mass drift.

### Measured convergence

The target star is 0.65 solar masses with a hydrostatic reference radius of **8,402 km**. The short baseline run covers 0.75 simulated seconds using a fixed 0.001-second step and damping time of 0.15 seconds. The extended run covers 15 simulated seconds with dt=0.01 seconds and tau=0.7 seconds.

| Particle count | Relax duration | Final/initial RMS acceleration | Density log-RMS before → after | Largest particle radius |
|---:|---:|---:|---:|---:|
| 64 | 0.75 s | 0.72585 | 1.13361 → 1.04386 | 7,309 km |
| 128 | 0.75 s | 0.71469 | 0.77733 → 0.70090 | 7,550 km |
| 128 | 15 s | 0.01504 | 0.77733 → 0.24684 | 5,129 km |
| 256 | 15 s | 0.02493 | 0.59274 → 0.20446 | 5,888 km |
| 512 | 15 s | 0.02261 | 0.55482 → 0.20929 | 6,609 km |

This experiment shows a large reduction in force residual after longer damping, accompanied by **excessive numerical contraction** of the low-resolution star. Even at 512 particles the maximum particle radius is about 21% below the hydrostatic reference radius. It would be misleading to call the relaxed particle configuration a physically calibrated equilibrium white dwarf. The residual speed and density mismatch remain significant.

The underlying pressure closure is the correct ideal zero-temperature electron law. The error arises from the coarse **particle discretization, fixed smoothing length, gravitational softening, approximate sampling and relaxation protocol**. The longer run is therefore an *equilibrium-seeking numerical experiment*, not a completed high-fidelity initial condition.

The stellar lab now exposes the **H** key to regenerate the hydrostatic seed, relax it in isolation, and then place it beside the black hole. The UI shows force-reduction and density-mismatch diagnostics. The operation is explicit because a longer damped run trades a smaller force residual for a visible structural bias.

### Required gates before astrophysical use

1. Build an adaptive smoothing length and a consistent pressure/gravitational-softening treatment, with appropriate variable-h corrections.
2. Increase particle count and test radial density against the hydrostatic reference. A target acceptance gate is **relative radius bias below 5%** and a density-profile error decreasing under resolution refinement.
3. Test relaxation over several dynamical times with a decreasing damping schedule. Evolve an isolated relaxed star *without damping* for multiple dynamical times to verify that radius, internal kinetic energy and density profile remain stable.
4. Only then insert the star into a companion black-hole field and test tidal stripping, capture boundaries and pressure/energy exchange.
5. Profile the O(N²) CPU algorithm against spatial indexing or GPU compute; the current interactive demo is deliberately small.

## 3. Source tree and acceptance

- `src/kerr_camera.rs`: new time-labelled equatorial Kerr null path sampler.
- `src/kerr_optics.rs`: circular test-body orbit, bounded source histories, moving-source worldtube intersection and causal acquisition.
- `src/bin/kerr_optical_lab.rs`: equatorial optical link demonstration with time controls.
- `src/stellar_relax.rs`: isolated cold-degenerate SPH relaxation and diagnostics.
- `src/bin/white_dwarf_lab.rs`: optional damped initialization in the existing stellar demo.
- `src/bin/optical_relax_bench.rs`: numerical ray link and relaxation convergence results, with optional CSV output.

CI exercises the complete Rust tests, the Kerr photon link and the stellar relaxation benchmark, then builds Windows executables. Numerical acceptance checks cover speed-of-light causal bounds, photon radial/polar constraints, Shapiro-like delay of a radial Schwarzschild photon, source emission precedence, finite/limited history, conservation of bulk mass and center during relaxation, and repeatable convergence diagnostics.

No original recovered fiction or historical game-design source is rewritten by this milestone.
