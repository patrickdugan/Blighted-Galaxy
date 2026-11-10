# BG-PHY-25 | One-sided free-boundary enthalpy reconstruction

**Date:** 10 November 2026  
**Parent:** [BG-PHY-24](RADIAL_HOMOLOGY_VIRIAL.md)  
**Tracking:** [Issue #29](https://github.com/patrickdugan/Blighted-Galaxy/issues/29)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Scientific question and experiment boundary

BG-PHY-24 showed that the exact cold-EOS plus Poisson-matched
gravity Hamiltonian has an **outermost-shell pressure deficit** even
when the initial star's *global* radial homology derivative prefers
expansion. At N256 with angular glass, the M4 original-density
outmost-shell pressure/reference gravity is 0.731 versus inward
gravity/reference 0.942. For C4 with the BG-PHY-20 surface-gated MLS
density, the corresponding figures are 0.661 and 0.928.

The existing density correction improves the second-outermost
shell and damped radial oscillations but allows the post-damping
star to contract 20–22.5% below the independently integrated
cold-electron hydrostatic radius.

This milestone tests an unlabeled one-sided **specific-enthalpy**
boundary reconstruction. It is initially a scalar-energy and
finite-difference virial diagnostic. It does not provide or
activate a new analytic force in the playable or research
integrators until its quality is demonstrated and its full
derivative has been checked.

## Physics: enthalpy vanishes at the stellar vacuum surface

For the unchanged zero-temperature electron-degenerate barotropic EOS,

```text
dH = dP/rho,
a_pressure = -grad H,
H(rho=0) = 0.
```

BG-PHY-19 found that the exact reference enthalpy is
reconstructed accurately near the free surface, where pressure
itself is strongly nonlinear. This motivates a linear local
enthalpy fit incorporating an estimated vacuum surface location,
rather than dividing a reconstructed pressure gradient by a
small local density.

The candidate uses **only** current particle coordinates, masses,
baseline exact-root SPH densities, the selected density kernel and
the existing EOS. The independent 1D white-dwarf hydrostatic
profile appears exclusively in scoring.

## Boundary geometry inferred from truncated kernel moments

At each particle i, using h_i from the existing safeguarded root,

```text
V_j = m_j/rho_SPH,j
q_ij = (x_j-x_i)/h_i
w_ij = V_j W(|x_i-x_j|,h_i)

S_i = sum_j w_ij
Q_i = sum_j w_ij q_ij
n_out,i = -Q_i/|Q_i|.
```

Q_i points toward the weighted center of the observed kernel
support. Near a free surface it therefore estimates the inward
side of the particle cloud. The opposite direction gives the
estimated local outward surface normal. The first-moment norm
and normal orientation are measured in the benchmark.

To estimate the distance to a locally planar vacuum boundary,
consider an ideal normalized isotropic kernel with a half-space
containing all volume on one side of the plane. A particle
inside that half-space by distance a*h sees the fraction

```text
S_flat(a) = 4*pi integral_(q=0)^2 q² W(q,1)
             * f(q,a) dq,

f(q,a)=1                     for q<=a,
f(q,a)=(1+a/q)/2            for q>a,

0<=a<=2.
```

Thus S_flat(0)=1/2 and S_flat(2)=1. BG-PHY-25 uses a
128-bin midpoint radial quadrature normalized by S_flat(2),
constructs a 65-node monotone lookup for each kernel, and
inverts the measured partition S_i to estimate a_i in units
of smoothing length. It clips S_i below 0.5 to a_i=0 and
above 1 to a_i=2.

The geometry is a *flat-surface approximation*. Particle
disorder, large radial curvature, and a density gradient
can change the measured partition independently of geometric
truncation. The inferred distance is a hypothesis to test,
with fallback and quantitative diagnostics.

## Candidate enthalpy reconstructions

With H_j=H_cold(rho_SPH,j) and the affine basis

```text
phi(q)=[1,qx,qy,qz],
M_i=sum_j w_ij phi(q_ij) phi(q_ij)^T,
b_i=sum_j w_ij phi(q_ij) H_j,
beta_i=M_i^-1 b_i.
```

The ordinary affine enthalpy fit predicts
`H_affine,i=beta_i[0]`.

A second one-sided fit adds a **virtual vacuum point**

```text
q_ghost = a_i n_out,i,
H_ghost = 0,
w_ghost = clamp(1-S_i,0,0.5).
```

Its normal equations add `w_ghost phi(q_ghost)
phi(q_ghost)^T` to M_i, with zero contribution to b_i.
This uses no 1D hydrostatic target density, and the full
kernel geometry determines all ghost parameters. The EOS
maps each positive enthalpy intercept back to a density:

```text
rho_H,i = rho_cold(H_affine,i)
rho_G,i = rho_cold(H_ghost,i).
```

BG-PHY-20's original smooth support gate is kept fixed:

```text
t_i=clamp(2(1-S_i),0,1),
alpha_i=t_i²(3-2*t_i),

rho_H,blend,i=(1-alpha_i)rho_SPH,i + alpha_i rho_H,i
rho_G,blend,i=(1-alpha_i)rho_SPH,i + alpha_i rho_G,i.
```

The screening candidates are:

| Choice | Density entering the scalar cold-EOS energy |
|---|---|
| Original | exact-root SPH rho |
| Surface MLS | BG-PHY-20 blended affine rho |
| Enthalpy affine blend | rho_H,blend |
| Vacuum ghost blend | rho_G,blend |

Any ill-conditioned, nonpositive or nonfinite enthalpy fit
uses the original baseline rho_i for that candidate and
records its fallback fraction. The 4×4 normal equations
use diagonal preconditioning and a Cholesky pivot-spread
limit of 1e8, matching the existing MLS contract.

## Follow-on test: first-to-zeroth kernel moment ratio

The first four-job N64/N128/N256 matrix passed. On N256 seed0, the
partition-inferred vacuum ghost worsened last-shell pressure support
and moved the discrete homology energy minimum inward. At fixed C4
the original density's last-shell pressure/reference-g was 0.681;
the BG-PHY-20 surface MLS yielded 0.661, ordinary affine-enthalpy
blend 0.679, and partition-only vacuum ghost 0.568. The ghost
raised last-shell density log-RMS from 0.712 (original) to 1.407
and shifted the preferred homology radius bias to -11.16%, compared
with -2.85% under original SPH. Seed17 independently reproduces
the poor ghost response.

The inferred outward normal itself aligned well with the stellar
radial direction (0.997 mean cosine in the N256 C4 last shell).
The last shell's measured volume partition was only 0.466,
already **below the 0.5 half-space ideal**. The first ghost model
therefore placed the vacuum boundary at zero distance from many
nonzero-enthalpy interior particles. This exposes a confound:
the volume-weighted partition includes errors in SPH volume
estimation as well as actual missing geometric support.

A second parameter-free geometric hypothesis uses the ratio of
the **first and zeroth kernel moments** to eliminate a common
multiplicative source-volume bias. For a locally planar boundary
at depth `a*h`, the truncated normal first-moment magnitude is

```text
Q_flat(a) = pi integral_(q=a)^2
              q^3 W(q,1) [1-(a/q)^2] dq.

R_flat(a) = Q_flat(a)/S_flat(a),
R_observed = |Q_i|/S_i.
```

If each local source volume has the same factor `c` error,
then `S_observed=c S_flat(a)` and
`|Q_observed|=c Q_flat(a)`, so their ratio is invariant.
The new candidate inverts the monotone function R_flat(a)
using a kernel-specific 65-node lookup. Its ghost weight is
the **ideal missing geometric volume**
`1-S_flat(a_ratio)`, independent of the observed
partition's multiplicative bias. The outward normal continues
to come only from the observed Q direction.

The fifth screening candidate is
`moment_ratio_ghost_blend`, with the same cold EOS,
affine enthalpy fit, positivity checks and BG-PHY-20
smooth surface exposure as before. The kernel
shape, lookup resolution and ratio inversion are fixed
across both particle-orientation seeds, M4 and C4.
A synthetic test checks that inferred depth is unchanged
under source-volume scalings 0.63, 0.81 and 1.0 for
several actual half-space depths.

This is a *new exploratory hypothesis motivated by the
observed failed partition-only ghost*, and its evidence
should be identified separately from the original
preregistered four-candidate experiment. It may still
fail if the SPH volume bias varies systematically over
the support or the physical surface departs strongly
from a plane. Its N256 holdout seed17 comparison is
retained and no reference stellar label enters the formula.

## Synthetic verification

Before white-dwarf screening, unit tests require:

- Normalized C4/M4 half-space kernel fraction 0.5 at the
  surface and 1 when the support lies fully inside.
- Monotone invertibility and accurate recovery of known
  planar support depths.
- Exact affine enthalpy reproduction with a synthetic
  irregular one-sided point cloud satisfying H=0 on the
  supplied vacuum plane, as well as ordinary constant-field
  reproduction in the absence of the ghost constraint.
- Finite, positive reconstruction for a synthetic curved
  enthalpy profile, and detection of the correct outward
  direction from the one-sided kernel first moment.
- Positive finite density fields on a real initialized
  cold white-dwarf point cloud.

A vacuum ghost is deliberately a **physical boundary
condition**, so a positive constant enthalpy in the entire
half-space is not a compatible test case for the ghost.
The interior affine fit must reproduce constants exactly.

## Energy-based screening, preserving conservation

For a frozen configuration x, each candidate defines a
physical scalar total potential,

```text
V_candidate(x) =
  sum_i m_i u_cold(rho_candidate,i(x))
  + Phi_G(x).
```

Every evaluation completely re-solves SPH h/rho, volumes,
the one-sided geometry and any MLS coefficients. The matched
Poisson gravitational pair potential is unchanged and uses
the original h(x).

**Radial force without substituting a nonconservative pair
law:** displace only one equal-mass shell k by a small physical
distance epsilon*Rref in its outward radial unit vectors:

```text
x_i^+ = x_i + epsilon Rref rhat_i   for i in shell k
x_i^- = x_i - epsilon Rref rhat_i   for i in shell k

sum_(i in k) F_P,i dot rhat_i
  = -[U(x^+) - U(x^-)]/(2 epsilon Rref).
```

The same operation on Phi_G gives inward gravity, normalized
against `sum_(i in shell) m_i G M_ref(<r_i)/r_i²`.
For each shell a second proportional-radius perturbation
measures the pressure and gravity virial contributions,
and the eight-shell virial sum is independently comparable
with BG-PHY-24.

This is a computationally expensive **finite-difference
diagnostic** of the complete cold-EOS Hamiltonian, and is
not an analytic candidate force yet. It preserves a
variational reference for every new density candidate.

The benchmark also scans a radial homology bracket
0.75<=s<=1.30 and refines negative-to-positive dV/ds
crossings, reporting an explicit unbracketed result if
none exists. Both initial star positions and the
previous solver's 3s damped positions can be evaluated.
On post-damped positions, *every candidate* is tested
at identical coordinates to distinguish direct
density effects from geometric changes.

For the old SPH and BG-PHY-20 blended density, a unit
test compares this independent shell-direction scalar
derivative against the established analytic shell force
from BG-PHY-24 at N32.

## Measurements and decision gate

The dedicated [BG-PHY-25 workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy25.yml)
tests N64 and N128 pilots, then N256 orientation seeds
0 and 17 at matched M4 eta1.20 / C4 eta1.70.
Every CSV records the eight shell pressure/gravity
ratios, last and penultimate shell support, original
kernel partition, inferred outward-normal alignment,
density log-RMS against the independent reference,
quadratic enthalpy-gradient error, homology minimum,
radius bias, EOS+gravity energy and fallback/cap
fractions.

A candidate is eligible for Stage C analytic
force derivation only if it improves the **last
mass shell** and the implied radial homology minimum
for both M4/C4 and independent orientations, without
displacing the second-outermost layer into a stronger
inward deficit or causing substantial fit fallbacks.

Stage C must derive the new candidate rho(x)'s total
analytic derivative through SPH h, volume weights,
boundary normal, planar support inversion, ghost
location, enthalpy and MLS inverse. It must then
match an independent complete scalar energy gradient
at N16/24/32 before any N256 gravity-coupled
physical dynamics. The original 3s damped and
7.88s undamped Issue #8 acceptance thresholds
remain binding, with no gameplay force substitution
or EOS tuning in the present milestone.

## References

- [BG-PHY-24 measured homology and virial](RADIAL_HOMOLOGY_VIRIAL.md).
- [BG-PHY-23 conservative gravity plus pressure Hamiltonian](GRAVITY_COUPLED_HAMILTONIAN.md).
- [BG-PHY-20 blended MLS density study](FREE_SURFACE_DENSITY.md).
- Smoothed particle hydrodynamics for free-surface and
  multiphase flows: a review (2025), DOI 10.1088/1361-6633/ada80f.
- *The Mathematics of Smoothed Particle Hydrodynamics
  Consistency* (2021), DOI 10.3389/fams.2021.797455.
