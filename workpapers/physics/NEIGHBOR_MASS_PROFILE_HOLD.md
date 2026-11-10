# BG-PHY-11: Compact-SPH neighbor indexing and high-resolution equilibrium audit

**Research milestone:** 10 November 2026  
**Program:** Blighted Galaxy clean-room Rust physics reconstruction  
**Tracked acceptance issue:** [#8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)  
**Prior work:** [BG-PHY-10 M4 Poisson gravity](M4_GRAVITY_CUBIC_SHELLS.md)  
**Code:** `sph_neighbors.rs`, `sph_massprofile.rs`, `sph_adaptive.rs`, `bin/phy11_bench.rs`

## 1. Problem and source of truth

A cold electron-degenerate white dwarf should retain the hydrostatic reference structure while isolated, without artificial damping. The independent 1D Newtonian hydrostatic solver, using a zero-temperature relativistic electron equation of state with `mu_e=2`, gives a **0.65-solar-mass star with radius approximately 8,402.21 km** and dynamical time

    t_dyn = sqrt(R^3 / (GM)) ≈ 2.62223 s.

Earlier BG-PHY-10 trials reached 1,024 particles, but used only 0.6 s of damping (about 0.23 dynamical times). This short profile did not establish an equilibrium. For the same Poisson-matched M4 potential, the N=128 and N=256 stars contracted by 10.97% and 14.39%, respectively, when evolved undamped for three dynamical times.

The next acceptance gate requires all of the following measured on the **same** high-resolution particle realization:
- post-relaxation outer radius within 5% of the 1D hydrostatic result;
- an isolated undamped hold of at least three dynamical times with at most 5% radius change, including the outermost and 99%-enclosed-mass radii;
- improvement of the radial mass-density distribution under increasing N and a documented energy/momentum invariant audit.

The 5% radial-CDF RMS target is an additional numerical development criterion; its definition below is distinct from the original Issue #8 acceptance target.

## 2. Deterministic cell-list neighbors

**References:** Price & Monaghan (2007), [adaptive gravity and per-particle neighbor search](https://academic.oup.com/mnras/article/374/4/1347/963634); [SPH variable-h neighbor-search techniques](https://pmc.ncbi.nlm.nih.gov/articles/PMC9810768/).

For the 3D M4 density kernel, interactions have compact support at distance `r < 2h_i`. The new spatial index stores particle indices in a hash map keyed by integer grid cell coordinates:

    cell(p) = floor(p / edge)
    edge = 2 * max_i(h_i).

For each particle, checking the center cell and its 26 immediate neighbors includes every parcel inside the maximum kernel support. Candidate lists are sorted into original particle-index order to preserve deterministic summation. Negative and large positions are supported.

Adaptive smoothing lengths change during the iteration `h_i=eta*(m_i/rho_i)^(1/3)`. The cell list is therefore rebuilt for **each fixed-position density iteration**, using the *current* maximum h, and again for the final density/Omega evaluation. Using the configuration's global upper h bound as a cell edge would have made the cells so large that nearly all stellar particles would appear as neighbors. The compact M4 gravitational zeta sums use the same spatial candidate principle. Plummer softening has infinite support, so its zeta sum remains all-pairs.

A test checks each compact-support pair is present in the list and compares all indexed densities and grad-h correction factors with the original direct pairwise sum to tight floating-point tolerance. This preserves the previous physical equations. It only changes which known-zero kernel terms are skipped.

**Computational boundary:** Global stellar gravity, the full pairwise energy diagnostic, and force loops still scale as O(N²). This optimization accelerates compact-support density and zeta work; a tree/FMM or mutually symmetric tree approximation is required before large-scale stellar calculations become inexpensive. Cell-list preparation and hash lookups also carry overhead. The headless benchmark measures the average number of candidates and single-pass indexed/dense timing on its actual particle distribution.

## 3. Radial mass-profile diagnostics

An SPH star's maximum particle radius can be dominated by a single stray tracer. A radius acceptance test should still report that outer extent, while also measuring how much mass actually resides at each radius.

Let `r_i=|x_i-x_COM|`, with `x_COM = sum_i(m_i x_i) / sum_i m_i`. Sort `r_i` and construct the midpoint mass-CDF `F_i=(sum_{j<i}m_j + m_i/2)/M`. Linear interpolation in this CDF gives the mass-weighted radial percentiles `R50`, `R90`, `R95`, and `R99`.

The 1D reference profile provides the inverse mass-CDF `R_ref(f)`. The new diagnostic reports

    RMS_CDF = sqrt( (1/19) sum_{k=1..19}
                [(R_particle(k/20)-R_ref(k/20))/R_ref(1)]^2 ).

The CDF mismatch is dimensionless and normalized to the hydrostatic stellar radius. The benchmark tracks the worst mismatch throughout the undamped run. A 5% threshold here is intentionally conservative as an additional research criterion.

Unit tests verify mass and COM preservation, monotonic percentiles, and sensitivity to uniform radial expansion. These diagnostics are model-agnostic and can be reused for future neutron stars and stellar remnants.

## 4. Reproduction

The new headless program runs the independent Poisson-matched M4 gravity and nearly cubic shell initializer:

    cargo run --release --bin phy11_bench --       --n 1024 --relax-steps 300 --hold-steps 788       --csv bg-phy11-highres.csv

The default physical step is 0.01 s. The code first applies 300 damped relaxation steps (3 simulated seconds), then disables damping for 788 steps (7.88 s, approximately 3.005 dynamical times). It reports:

- initial and post-relaxation mass-CDF mismatch;
- residual accelerations and outer radius;
- density-neighbor candidate fraction and direct/indexed density agreement;
- Rmax, R99, R95, R90, R50 and CDF mismatch throughout the undamped hold;
- maximum relative radial drift during the hold, including intermediate snapshots;
- energy and momentum drift, and a separate explicit physical gate verdict.

The numerical benchmark keeps the invariant assertions (mass balance, finite/controlled energy drift, momentum). **It writes a physics-gate failure in the CSV if the star remains out of equilibrium**, while permitting the software CI job to succeed when its numerical implementation is working correctly. A failed physical gate is not proof of a software or compilation defect; it remains an open scientific issue.

## 5. CI and artifacts

The lightweight CI runs a 64-particle smoke check, all Rust tests and the Windows release build. The separate **BG-PHY-11 high-resolution** GitHub Actions workflow tests 1,024 particles over the three-dynamical-time hold and publishes `bg-phy11-highres.csv`. The high-resolution job has a fixed CI timeout and is a regression/measurement gate that can be rerun manually.

**Source provenance:** The Blighted Galaxy story archive motivates the desired stellar destruction and gravitational-shearing encounters. The numerical EOS, M4 force, cell index and mass-profile diagnostics are new mathematical and engineering contributions. The game still treats large-scale stellar weapons as explicitly fictional effects until a stable self-gravitating star model and a realistic disruption model are validated.

## 6. Follow-on work

- Measure candidate-list hit rates and the CPU cost of dynamic density iteration across N=1,024, 2,048 and 4,096 under a fixed physical window.
- Build a verified, spatially indexed long-range tree or fast multipole gravity backend, taking care with pair symmetry and total momentum accounting.
- Study stable stellar initial conditions through glass/HEALPix-like equal-volume sampling and larger N, with a damping schedule gradually reduced over several dynamical times.
- Extend undamped evaluation to radial shell oscillation modes, density-profile evolution and virial-balance diagnostics.
- Only after the independent equilibrium gates pass, evolve a self-gravitating white dwarf through a black-hole tidal encounter with an appropriate relativity/capture boundary.
