# BG-PHY-16: Wendland C4 and mass-CDF-preserving angular glass

**Date:** 10 November 2026  
**Predecessors:** [BG-PHY-14](NEIGHBOR_GROWTH_HROOT.md), [BG-PHY-15](FROZEN_FORCE_AUDIT.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## 1. Hypothesis and experimental controls

BG-PHY-15 measured dominant initial pressure force errors even when the
matched M4 gravitational force agreed with the 1D hydrostatic field to
a few percent. It also identified low outward pressure gradients in the
outer mass shells, while the particlewise pressure noise is strongest
in the interior. These phenomena may come from both particle disorder and
the M4 kernel's finite-neighbor consistency / pairing limitations.

The experiment separates two controlled factors:

1. Angular geometry: original approximately cubic-shell antipodal seeds
   versus a deterministic, geometry-only angular-glass relaxation.
2. Kernel: matched M4 versus normalized Wendland C4, using the **same**
   kernel in density, pressure, grad-h and the Poisson-matched gravitational
   potential. The independent cold-electron EOS is unchanged.

Each pair is crossed with fixed eta=1.20 and
eta(N)=1.05*(N/256)^(1/6). N256 and N512 each run shell seed 0 and
independent shell-orientation seed 17. N1024 runs seed 0. Each configuration
uses the same 0.65-solar-mass, mu_e=2 hydrostatic reference.

The default experiment only reads the force field at fixed positions.
A separate N256 experiment applies the existing 3-second damping step to
all eight seed/kernel/eta combinations from the same base seed.
Unchanged EOS, gravity G, particle masses, and reference star make
the within-kernel seed comparison causal with respect to angular geometry.
Comparisons across kernels include different effective interpolation
lengths, which are recorded via mean h and realized support neighbors.

## 2. Wendland C4 density and derivatives

For q=r/h, q in [0,2), and h measured as in the existing M4:

```text
W_C4(r,h) = 495 / (256*pi*h^3)
             * (1-q/2)^6 * (1 + 3q + 35q^2/12).

dW/dr   = [495 / (256*pi*h^4)]
           * [-(7/3)*q*(2+5q)*(1-q/2)^5].

dW/dh   = -[495 / (256*pi*h^4)] * (3F + q*F'), 
           where F=(1-q/2)^6*(1+3q+35q^2/12).
```

All derivatives vanish outside compact support as applicable. The
normalization is exactly unity in the 3D spherical measure
4*pi*integral_0^2h r^2 W(r,h)dr. Gamma=5/3-like cold electron
degeneracy physics is supplied by the existing independent EOS.

The opt-in adaptive configuration `with_wendland_c4()` chooses this
kernel and keeps the existing underrelaxed h/rho closure. The Newton
solver gains a kernel argument while preserving the previous M4-only API.
The force law uses the matched Omega factor based on the same density kernel:

```text
rho_i(h_i) = sum_j m_j W_C4(r_ij,h_i)
Omega_i    = 1 + h_i/(3*rho_i) * d(rho_i)/dh_i.
```

## 3. Matching gravity via Poisson's equation

An SPH density kernel and gravitational softening must agree if a
controlled comparison is to retain the variational grad-h correction.
For a unit mass C4 kernel, define q=r/h and

```text
M_C4(<q) = q^3/16384 * [
    42240 - 59136 q^2 + 79200 q^4 - 73920 q^5
    + 30800 q^6 - 6336 q^7 + 525 q^8
],      0 <= q <= 2.
```

For 1<q<2, code uses a stable factorized complementary mass polynomial
to avoid subtractive cancellation. `M_C4(<2)=1` exactly.
The potential is continuous at q=2, with `phi=-1/r` outside support.
Within q<2 its analytical primitive is

```text
phi_C4(r,h) = (1/h)*[
    -55/32 + (165/128)q^2 -(231/256)q^4
    +(825/1024)q^6 -(165/256)q^7
    +(1925/8192)q^8 -(11/256)q^9
    +(105/32768)q^10
].

dphi/dr = M_C4(<q)/r^2.
dphi/dh = 55*(2-q)^7*(3q+2)*(7q^2+8q+4)
           /(32768*h^2).
```

Tests compare the analytical derivatives with finite differences in h
and r, verify 3D density normalization, and ensure the preexisting
M4 path reproduces its own kernel, radial derivative and potential.

## 4. Deterministic angular-glass initializer

Each pair in the approximately cubic-shell seed has physical positions

```text
x_i,+ = C + r_i*u_i,
x_i,- = C - r_i*u_i,
|u_i| = 1.
```

The new geometric optimizer adjusts unit directions u_i while holding
r_i fixed. Its pair interaction is a symmetric Gaussian repulsion

```text
U = sum_{i<j} [
   exp(-|r_i*u_i - r_j*u_j|^2 / (2*sigma_ij^2))
 + exp(-|r_i*u_i + r_j*u_j|^2 / (2*sigma_ij^2))
]

sigma_ij = (h_ref,i + h_ref,j)/2,
h_ref,i  = clamp(eta*(m_i/rho_ref(r_i))^(1/3), h_min, h_max).
```

Projected tangent-space gradient descent with monotone line search
minimizes U. It preserves the entire radial mass CDF, the center of mass,
initial momentum and kinetic velocities. Glass moves change within-shell
and inter-shell angular correlations, unlike a rigid shell rotation.

This objective is geometric; it is *not* an optimization against the SPH
force-audit metric. Its effects on pressure and gravity errors are
therefore informative even when there is no improvement.

The benchmark records its objective before and after accepted sweeps,
the maximum pair rotation in degrees, and the exact initial Rmax
sampling offset at enclosed mass fraction 1-1/N.

## 5. Metrics and acceptance

Use the BG-PHY-15 `frozen_force_audit` without advancing particles.
The summary includes pressure/gravity/net vector errors against the
independent 1D hydrostatic gravity, radial pressure/gravity ratios,
tangential error components and outer-shell pressure support.

The companion shell CSV provides all eight equal-mass-shell diagnostics.
The initial exact sampled radius, radius-CDF error, mean smoothing length,
realized neighbors within 2h_i and angular glass objective are also
reported. The optional N256 3-second damping pilot measures radius
contraction relative to the initial sampled radius and repeats the
force snapshot.

Physical star acceptance remains the same as BG-PHY-14 / Issue #8,
which requires full undamped 3-dynamical-time stability after the
damping preparation and a radius/profile accuracy gate. This experiment
is a diagnostic and a short dynamical pilot.

## 6. Reproduction

```sh
cargo test --all-targets

cargo run --release --bin phy16_bench -- \
  --n 256 --seeds 0,17 --glass-sweeps 20 \
  --shells 8 --csv bg-phy16-n256.csv \
  --shell-csv bg-phy16-n256-shells.csv

cargo run --release --bin phy16_bench -- \
  --n 256 --seeds 0 --glass-sweeps 20 \
  --relax-steps 300 --dt 0.01 --shells 8 \
  --csv bg-phy16-n256-dynamics.csv \
  --shell-csv bg-phy16-n256-dynamics-shells.csv
```

The [BG-PHY-16 Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy16.yml)
runs the N256/N512/N1024 frozen factorial and N256 dynamical pilot.
The regular CI runs the N64 short smoke and compiles the new Windows binary.

## 7. Measured reference results (10 November 2026)

The [factorial matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42765895092)
passed N256, N512 and N1024 frozen comparisons and an N256
matched-duration 3-second damping pilot. See each Actions job's
summary and shell CSV artifacts for full precision. The following
values use orientation seed 0 and fixed eta=1.20, except where noted.

| N | kernel | original pressure vector RMS | glass pressure vector RMS |
|---:|---|---:|---:|
|256|M4|0.3230|0.1798|
|512|M4|0.3817|0.2074|
|1024|M4|0.4399|0.2747|
|256|Wendland C4|0.9375|0.6686|
|512|Wendland C4|1.1413|0.8846|
|1024|Wendland C4|1.4156|1.2244|

Glass improves every case in this table, including the inner-shell
pressure error. The current constrained geometric objective reduces
pressure noise even though it does not directly optimize the SPH
equilibrium equation.

The N256 seed 0 3-second dynamical pilot measures:

| eta law | layout | kernel | contraction vs seed Rmax | absolute Rmax bias | post net/reference force RMS |
|---|---|---|---:|---:|---:|
|fixed|shell|M4|−11.66%|−17.56%|0.0457|
|fixed|glass|M4|−11.06%|−17.00%|0.0387|
|fixed|shell|C4|−7.66%|−13.83%|0.0793|
|fixed|glass|C4|−8.49%|−14.60%|0.0703|
|growth|shell|M4|−10.26%|−16.26%|0.0578|
|growth|glass|M4|−10.21%|−16.21%|0.0494|
|growth|shell|C4|−3.58%|−10.03%|0.1371|
|growth|glass|C4|−4.76%|−11.12%|0.1278|

The C4 cases produce less contraction in this short pilot, while their
post-damping net-force mismatch exceeds the M4 controls. None of these
configurations passes the independent 5% reference-radius gate.

The new closure audit reveals why the most undersampled C4 cases need
independent root checks: for N256, growth eta=1.05, seed0, the C4
underrelaxed closure reaches max uncapped residual 0.01153, compared
with 1.17e-5 for M4. At fixed eta=1.20, C4 max residual is 6.65e-5.
All tested configurations have no hard h caps or Omega-floor hits.

**Independent exact-root frozen comparison** (same positions, N256,
seed 0):

| eta | kernel | 24-pass pressure RMS | safeguarded Newton pressure RMS |
|---:|---|---:|---:|
|1.20|M4|0.323026|0.323027|
|1.20|C4|0.937546|0.937590|
|1.05|M4|0.582172|0.582176|
|1.05|C4|3.496443|3.684815|

Exact roots achieve <1e-9 dimensionless residual. The low-neighbor C4
force noise remains, and under exact root correction it increases by
about 5.4% in the growth-law case. The smallest Omega for that
configuration is about 0.178 under Newton, above the 0.02 floor.
Incomplete h/rho closure is therefore not the primary source of the
large pressure-gradient mismatch.

The C4 kernel's pairing resistance permits much larger neighbor sets,
so testing its intended high-neighbor operating regime requires a
separate eta sweep at matched effective support. Every comparison must
report realized neighbors and h, while retaining the original EOS and
white-dwarf equilibrium gates.

## 8. Scientific decision rule

Inspect pressure and gravity error, inner-shell tangential pressure RMS,
outer-shell pressure support, and post-damping radius contraction
together. A kernel lowering pressure error while maintaining a matched
gravitational field is a promising candidate for longer undamped holds.
A seed geometry lowering the same errors independently can become
a starting point for higher-resolution glass/discrepancy experiments.

As always, keep the original playable physics defaults until the
independent radius and stability criteria pass.

### References

- Dehnen & Aly (2012), *Improving convergence in smoothed particle
  hydrodynamics simulations without pairing instability*,
  DOI 10.1111/j.1365-2966.2012.21439.x.
- Zhu, Hernquist & Li (2015), *Numerical Convergence in Smoothed Particle
  Hydrodynamics*, DOI 10.1088/0004-637X/800/1/6.
- Price & Monaghan (2007), *An energy-conserving formalism for adaptive
  gravitational force softening*, DOI 10.1111/j.1365-2966.2006.11241.x.
