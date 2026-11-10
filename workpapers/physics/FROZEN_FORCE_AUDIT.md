# BG-PHY-15: Frozen hydrostatic force errors and finite-sampling radius bias

**Date:** 10 November 2026  
**Predecessor:** [BG-PHY-14](NEIGHBOR_GROWTH_HROOT.md)  
**Open physical gate:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)  
**Code:** `src/sph_force_audit.rs`, `src/bin/phy15_bench.rs`.

## Motivation

BG-PHY-14 tests fixed smoothing parameter eta=1.20 against
eta(N)=1.05(N/256)^(1/6). Exact implicit h-rho closure agrees with the
existing underrelaxed method to machine precision in the N256/N512 equilibrated
snapshots, although the cold white dwarf retains a substantial radius deficit.
To determine where the force mismatch comes from, BG-PHY-15 decomposes
pressure support, self-gravity and tangential particle noise at **frozen**
initial and post-damping positions.

The benchmark retains the same reference white dwarf:
0.65 solar masses, mu_e=2, R_ref approximately 8402.207 km, cold relativistic
electron-degeneracy pressure, matched cubic M4 gravity and variable-h terms.
No pressure rescaling, radius adjustment or artificial confining potential is
introduced.

## The initial particle radius is a quantile, not the continuum surface

The approximately cubic-shell sampler places antipodal equal-mass pairs at
enclosed-mass fractions

```text
f_i = (i + 0.5) / (N/2),     i=0,...,N/2-1.
```

Therefore, before any relaxation,

```text
R_max,seed = R_ref(M(<r)/M = 1 - 1/N).
b_sampling(N) = R_max,seed / R_ref - 1.
```

The outermost particle samples a mass quantile strictly below unity, so an
unperturbed exact radial seed already has a negative **last-particle radius**
bias. This term is computed directly from the independently integrated
hydrostatic inverse mass CDF, and verified against the actual oriented seed
for every run. Antipodal rotations preserve that quantile exactly.

Near a cold nonrelativistic-degeneracy surface,
P is proportional to rho^(5/3), so rho is proportional to (R-r)^(3/2).
The missing surface-shell mass is proportional to (R-r)^(5/2), giving

```text
R_ref - R_max,seed  ~  N^(-2/5)     (asymptotic near-surface law).
```

The *actual dynamical contraction* is reported separately:

```text
c_dynamics = R_max,post / R_max,seed - 1.
```

The existing absolute radius gate remains unchanged: the benchmark reports
both the bias against R_ref and contraction relative to the initialized
particle realization. The latter is a diagnostic, never a substitute
acceptance criterion.

## Vector force errors against independent 1D hydrostatics

For each particle i, let r_i be its radius about the discrete center of mass,
u_i its radial unit vector, and

```text
g_ref(r_i) = G M_ref(<r_i) / r_i^2.
```

The continuum hydrostatic pressure and gravity accelerations are
+g_ref u_i and -g_ref u_i. At fixed positions, the exact SPH
`forces` method supplies separate pressure, gravity and total accelerations,
including adaptive smoothing and gravity grad-h corrections.

Define the common mass-weighted reference norm

```text
S = sum_i m_i g_ref(r_i)^2.
```

Report dimensionless RMS vector errors:

```text
e_P   = sqrt(sum_i m_i |a_P,i - g_ref(r_i) u_i|^2 / S);
e_G   = sqrt(sum_i m_i |a_G,i + g_ref(r_i) u_i|^2 / S);
e_net = sqrt(sum_i m_i |a_P,i + a_G,i|^2 / S).
```

Also report radial mean ratios

```text
p_ratio = sum_i m_i (a_P,i dot u_i) / sum_i m_i g_ref(r_i);
g_ratio = sum_i m_i (-a_G,i dot u_i) / sum_i m_i g_ref(r_i).
```

and tangential RMS values for pressure, gravity and total accelerations,
also divided by sqrt(S). Each radial shell uses its own reference norm and
covers approximately one eighth of the **particle mass**. The shell CSV
reports the pressure/gravity ratios and vector/tangential errors to locate
bias in the interior or near the stellar surface.

Interpretation:
- Low p_ratio with g_ratio close to unity indicates insufficient resolved
  outward pressure acceleration.
- High g_ratio with p_ratio close to unity indicates excess resolved
  inward gravitational acceleration.
- Ratios near unity with large vector and tangential errors point to local
  angular noise that cancels in spherical averages.
- A small net error in a post-relaxation snapshot can coexist with large
  component errors and a contracted stellar radius, which is why the **initial
  seed snapshot** is a first-class measurement.

These tests diagnose numerical biases at the actual particle coordinates.
Their definitions do not assert that individual SPH pressure gradients
represent the continuum accurately at finite resolution.

## Reproduction

```sh
cargo test --all-targets

cargo run --release --bin phy15_bench -- \
  --n 256 --seeds 0,17 --fixed-eta 1.20 \
  --relax-steps 300 --dt 0.01 --shells 8 \
  --csv bg-phy15-n256.csv --shell-csv bg-phy15-n256-shells.csv

cargo run --release --bin phy15_bench -- \
  --n 64 --seeds 0 --relax-steps 6 --shells 8 \
  --csv bg-phy15-smoke.csv --shell-csv bg-phy15-smoke-shells.csv
```

The dedicated matrix workflow compares N256 and N512 with orientation seeds
0 and 17; N1024 uses seed 0 initially to limit all-pairs gravity cost.
The fixed and growing eta runs share exactly the same initial positions
at each N and seed. A short N64 CI smoke confirms finite diagnostics and CSV
generation. Snapshot audits preserve particle coordinates and velocities.
The relaxation phase uses the original 3.00s constant-tau preparation.

The output records absolute radius bias, known initial sampling bias,
post-damping contraction relative to the same initial seed, mass-CDF
difference, independent reference gravity RMS and all component force errors.
The experiments **do not** run a new undamped 3-dynamical-time stability hold.
BG-PHY-14 and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)
retain responsibility for full physical acceptance.

## Decisions after seeing the matrix

Inspect initial force errors before interpreting post-damping equilibrium,
where the system may already have collapsed to a different configuration.
If pressure dominates the error, compare a disorder-robust kernel such as
Wendland C4 at matched *realized* neighbor counts and evaluate the pressure
E0 error directly. If gravity dominates, examine the M4 softening and its
grad-h terms against the continuum field before modifying the EOS.
Run additional shell-orientation seeds when ensemble spread is comparable
to the measured resolution improvement.

### References

- Zhu, Hernquist & Li, *Numerical Convergence in Smoothed Particle
  Hydrodynamics* (2015), DOI 10.1088/0004-637X/800/1/6.
- Price & Monaghan, *An energy-conserving formalism for adaptive gravitational
  force softening* (2007), DOI 10.1111/j.1365-2966.2006.11241.x.
- Dehnen & Aly, *Improving convergence in smoothed particle hydrodynamics
  simulations without pairing instability* (2012), DOI 10.1111/j.1365-2966.2012.21439.x.
