# BG-PHY-21: Reproducing kernels and a compatible-work pressure pilot

**Date:** 10 November 2026  
**Issue:** [#22](https://github.com/patrickdugan/Blighted-Galaxy/issues/22)  
**Predecessor:** [BG-PHY-20](FREE_SURFACE_DENSITY.md)  
**Physical acceptance:** [#8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Scientific scope

BG-PHY-20 demonstrated that a boundary-gated affine MLS density field
reduces surface density log-RMS and frozen quadratic enthalpy-gradient
error across 36 tested C4 cases. The surface-gated result is a scalar
reconstruction, and its live position and smoothing-length derivatives
remain to be derived.

BG-PHY-21 implements the first part of a conservative reproducing-kernel
path:

1. First-order 3D reproducing-kernel coefficients and analytic query
   derivatives, checked independently against finite differences.
2. Anti-symmetric reproducing-kernel **pair pressure forces** based on
   the original SPH volume weights, with optional surface-gated
   reconstructed pressure samples.
3. A compatible pair-work update in a gravity-free, small-N pilot,
   guaranteeing discrete conservation of **kinetic plus evolved internal
   energy** and of linear momentum for the pilot's pairwise forces.
4. Numerical diagnostics for the still-unresolved cold-barotropic
   internal-energy consistency and angular-momentum drift.

The module is opt-in and separate from the playable dynamics.
It is intentionally restricted to pressure-only pilot systems; it
must not be scored against Issue #8 as a stellar equilibrium solver.

## First-order reproducing kernel

Use q_ij=(x_i-x_j)/h_i and frozen baseline SPH volumes
V_j=m_j/rho_SPH,j. At the query point i, define

```text
m0_i = sum_j V_j W_ij
m1_i = sum_j V_j q_ij W_ij
m2_i = sum_j V_j q_ij q_ij^T W_ij

B_i = -m2_i^{-1} m1_i
A_i = 1/[m0_i + B_i · m1_i]

W^R_ij = A_i [1+B_i·q_ij] W_ij.
```

When m2 is nonsingular and its condition estimate is acceptable,
this yields the discrete reproduction identities

```text
sum_j V_j W^R_ij = 1
sum_j V_j q_ij W^R_ij = 0.
```

Both identity residuals are measured directly on irregular seed clouds.
A singular or ill-conditioned moment matrix triggers an explicitly
counted fallback to the original kernel.

## Query-coordinate derivatives

The coefficient derivatives are needed when the *interpolation query*
moves with particle source coordinates, volumes and h held fixed.
For a component gamma,

```text
d_gamma B = -m2^{-1}[d_gamma m1+(d_gamma m2)B]
D = m0+B·m1
d_gamma A = -A²[d_gamma m0+(d_gamma B)·m1+B·d_gamma m1].
```

Since q=(x_i-x_j)/h_i, d_gamma q_alpha = delta_alpha,gamma/h_i.
The reconstructed gradient is

```text
grad W^R_ij =
  A (1+B·q) grad W
 + (grad A)(1+B·q) W
 + A W [B/h + (grad B)^T q].
```

The code forms the full three-dimensional moment derivatives of m0,
m1 and m2, and tests the gradient against an independent finite
difference of the complete query reconstruction (including A and B)
at fixed source particles. Rotation covariance is checked in M4
and C4, and the derivative of the partition identity is measured
as sum_j V_j h_i grad W^R_ij ≈ 0.

**Important limitation:** the query derivative is not yet the total
Lagrangian derivative with respect to moving particle coordinates.
That would also differentiate V_j, h_i(x), the selected MLS density,
S_i and its smoothstep gate, and the neighboring-particle corrections.
The full chain rule/implicit Jacobian is a later gate.

## Anti-symmetric pair forces

Define corrected pairwise gradient

```text
G_ij = [grad_i W^R_ij-grad_j W^R_ji]/2
G_ji = -G_ij.
```

Use the same **original SPH** volumes V_i in the RK moments and
the CRK pressure pair term:

```text
F_ij = -V_i V_j [P_i+P_j] G_ij,
m_i a_i = sum_(j!=i) F_ij.
```

P_i is evaluated from the unchanged zero-temperature Fermi EOS
at either original SPH density or the BG-PHY-20 surface-gated
affine MLS density. Because F_ji=-F_ij, total linear momentum
is exactly conserved in the discrete pairwise update up to floating
summation roundoff, regardless of whether the RK gradient is central.

Angular momentum has a separate residual

```text
dL/dt = sum_(i<j) (x_i-x_j) × F_ij.
```

This residual is generally nonzero for noncentral RK gradients.
BG-PHY-21 measures pair torque and angular momentum drift directly.

## Compatible thermal work

For an impulse of duration delta_t with all net velocity kicks applied
simultaneously, define the midpoint velocity

```text
v_i,mid = [v_i,before+v_i,after]/2.
```

The exact discrete pair kinetic-energy change sums to

```text
Delta K = sum_(i<j)
  (v_i,mid-v_j,mid) · F_ij delta_t.
```

The pilot divides each pair's opposite thermal work equally:

```text
m_i Delta u_i,ij = -0.5
   (v_i,mid-v_j,mid) · F_ij delta_t

m_j Delta u_j,ij = -0.5
   (v_i,mid-v_j,mid) · F_ij delta_t.
```

Therefore, Delta(K + sum_i m_i u_i) is zero to floating precision
for every compatible kick. The pilot uses two half kicks with a
drift between them, permits negative dt for time-reversal tests,
and records mass, momentum and total kinetic plus stored thermal energy.

**Thermodynamic consistency is a distinct unresolved obligation.**
The pressure is computed from P(rho_hat) at each force evaluation,
while the compatible internal-energy reservoir evolves by pressure
work. The resulting u_i can differ from the cold EOS's
u_EOS(rho_hat_i). BG-PHY-21 records that RMS discrepancy and makes no
claim that a coupled fluid continuity/energy law has been solved.
The pilot is gravity-free, so it makes no gravitational-energy claim.

## Verification

The first lightweight matrix covers N32 and N64 angular-glass seeds
and N128 shell seeds, each comparing M4/C4 and original/surface-blended
density at matched eta=1.7. The test suite verifies:

- exact constant and linear reproducing identities at interior and
  truncated-support locations;
- analytic query derivatives against finite differences;
- rotation covariance of reconstructed gradients and boost invariance
  of pressure accelerations;
- anti-symmetric pair force and total linear-momentum closure;
- discrete compatible K+U conservation after velocity-Verlet steps;
- reverse stepping (+dt then -dt), and measured angular momentum drift;
- cold EOS internal-energy residual after pressure-only evolution.

The new executable:

```bash
cargo test --lib sph_rk -- --nocapture
cargo test --lib sph_rk_compatible -- --nocapture

cargo run --release --bin phy21_bench -- \
  --n 32 --seed 17 --layout glass --eta 1.7 \
  --steps 16 --dt 0.001 --csv bg-phy21-n32.csv
```

The [BG-PHY-21 Actions workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy21.yml)
uses a small isolated matrix and uploads CSVs. Regular Rust/Windows
CI additionally compiles and smoke-tests the new binary.

## Remaining requirements for a stellar solver

The following are future gates before connecting this pilot to the
0.65-Msun white-dwarf dynamics:

1. Differentiate the BG-PHY-20 blend through S_i, all affine-MLS
   moment coefficients, rho_SPH, and the implicit h/rho closure.
   Include cap and fallback regime transitions, and test complete
   derivative against finite-difference *total* particle coordinate
   perturbations.
2. Derive a thermodynamically compatible density/pressure/energy
   evolution so the evolved u_i agrees with the cold barotropic
   EOS as resolution improves. The compatible-work identity alone
   does not establish this.
3. Add Poisson-matched gravitational potential and variable-h
   correction under an energy-consistent total potential or
   compatible gravity work. Restore angular-momentum regression tests.
4. Test finite-mass isolated N32/64/128 systems and confirm timestep
   convergence, independent angular seed robustness and fallback
   behavior.
5. Only then run the original N256/512/1024 stellar gate:
   3s damping, >=7.88s undamped, radius bias <5%, Rmax and R99 drift
   <5%, mass-CDF RMS <5%, and mass/momentum/energy checks.

Keep the original gameplay default untouched throughout this process.

### Literature

- Frontiere, Raskin & Owen, *CRKSPH - A Conservative Reproducing Kernel
  Smoothed Particle Hydrodynamics Scheme*, JCP 332 (2017), 160–209,
  DOI 10.1016/j.jcp.2016.12.004.
- Owen, *A compatibly differenced total energy conserving form of SPH*,
  Int. J. Numer. Meth. Fluids 75 (2014), 749–774,
  DOI 10.1002/fld.3912.
- PySPH, CRKSPH and corrected-gradient equations:
  https://pysph.readthedocs.io/en/main/reference/equations.html
- [BG-PHY-20 measured density workpaper](FREE_SURFACE_DENSITY.md).

## Measured N32/64/128 compatible-work results

The dedicated [BG-PHY-21 numerical matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42801881949) completed all N32, N64 and N128 jobs. Each runs M4 and Wendland C4 with both the baseline SPH density and the BG-PHY-20 surface-gated MLS candidate. Test systems run 12 or 16 pressure-only steps at dt=0.001s, then the same number of reversed steps. Each job published a four-row CSV.

Across 12 runs:

| Invariant or diagnostic | Maximum absolute/normalized measured error |
|---|---:|
| Invalid RK correction fraction | 0 |
| Constant reproduction error | 1.78e-15 |
| Linear reproduction error | 1.18e-15 |
| Full corrected query-gradient partition residual | 6.51e-15 |
| Relative K+stored-U drift | 1.70e-15 |
| Normalized linear momentum drift | 1.02e-18 |
| Normalized angular momentum drift | 1.12e-4 |
| RMS cold-u_EOS vs compatible-work u discrepancy | 3.53e-5 |
| Maximum reverse position error / reference radius | 3.59e-15 |
| Maximum reverse velocity error / 1e6 m/s | 1.06e-12 |

The separate lib tests verified the RK moment matrix identities, analytic
query-gradient finite differences, rotation covariance, pair momentum,
Galilean boost invariance and short-step energy/reversibility assertions.

**Frozen force-vector error against the 1D hydrostatic pressure acceleration:**

| N | M4 original | M4 blended | C4 original | C4 blended |
|---:|---:|---:|---:|---:|
|32|0.5382|0.5337|0.4435|0.4444|
|64|0.4864|0.4885|0.4695|0.4455|
|128|0.5524|0.5080|0.4963|0.4530|

These small-N comparisons establish the correctness of the RK
reproducing identities and compatible-work conservation. The initial
force field is still appreciably discrepant from the continuum pressure
reference. The blend's benefit is clearer at N128 than N32, and the
finite timestep and pressure-only scope limit extrapolation.

**Interpretation:** extremely small relative total-energy drift is a
consequence of the discrete midpoint-compatible work identity, not
proof that the evolved thermodynamic u matches the cold barotropic
u_EOS(rho_hat). The measured EOS consistency residual is small only
over short 12–16 ms trajectories and must be tested under longer,
nonlinear compression and expansion. Angular-momentum drift is
nonzero due to noncentral corrected kernel forces.

The next useful numerical oracle is the expensive total discrete
energy gradient of U(x)=sum_i m_i u_EOS(rho_hat_i(x)), evaluated by
finite-difference perturbations that recompute density, h, partition
gate and MLS coefficients at each changed configuration. Comparing
this variational gradient to the fast antisymmetric RK pair force
will measure the missing density/thermodynamic chain-rule terms
without assuming their magnitude.

## BG-PHY-21B: total cold-barotropic force oracle

The small-N compatible-work pilot conserves K+evolved-U exactly by
construction, but its P(rho_hat) and U_work evolution remain only
approximately cold-EOS-compatible. BG-PHY-21B builds an independent,
more expensive reference from the **total discrete cold-barotropic
internal energy**, with no independently evolved thermal reservoir.

For the same original SPH density or surface-gated MLS field define

```text
U(x) = sum_i m_i u_Fermi(rho_hat_i(x))
F_var,k = -dU/dx_k.
```

The EOS identity `du_Fermi/drho = P_Fermi/rho²` gives the exact
chain-rule expression

```text
F_var,k = -sum_i m_i [P_i/rho_hat_i²] d(rho_hat_i)/dx_k.
```

For the BG-PHY-20 blend,

```text
rho_hat_i = (1-alpha_i)rho_SPH_i + alpha_i rho_MLS_i,
d rho_hat_i =
   (1-alpha_i) d rho_SPH_i
 + alpha_i d rho_MLS_i
 + (rho_MLS_i-rho_SPH_i) d alpha_i.
```

In the partial-support regime `0<t_i<1`,

```text
t_i = 2(1-S_i)
alpha_i=t_i²(3-2t_i)
d alpha_i = 6*t_i*(1-t_i)*(-2 d S_i).
```

The partition factor `S_i=sum_j (m_j/rho_SPH,j) W_ij(h_i)`
also depends on particle positions and on the h/rho closure.
For an uncapped smoothing length defined implicitly by

```text
h_i³ rho_SPH,i(h_i,x)=eta³ m_i,
```

implicit differentiation gives

```text
d h_i/dx_k = -[h_i/(3 rho_SPH,i Omega_i)]
             * (partial rho_SPH,i/partial x_k) at fixed h,
Omega_i=1+(h_i/(3 rho_SPH,i))*partial_h rho_SPH,i.
```

At a hard h cap the derivative is piecewise constrained, and
changes of cap or MLS fallback introduce non-smooth transitions.
For the affine-MLS coefficient vector `beta_i=M_i^{-1} b_i`,
the full coefficient derivative is

```text
d beta_i = M_i^{-1}[d b_i-(d M_i) beta_i].
```

**What is implemented now:** central differences of the entire
scalar U(x) with respect to each physical particle coordinate.
For every +/− perturbation the code reconstructs the exact
safeguarded h/rho roots, original kernel volumes, partition S,
affine MLS and surface gate, then recomputes U from the cold EOS.
This acts as a full finite-difference oracle for the displayed chain
rule without prematurely claiming an analytic implementation.

BG-PHY-21B measures the RMS difference of this total barotropic
force from the fast antisymmetric RK pair force, net force and
net torque of the unprojected variational force, finite-difference
step-halving consistency, and short pressure-only velocity-Verlet
energy drift and reverse-step errors.

Because U is a scalar invariant under translation and rotation of
the entire particle constellation, its **exact** gradient is
momentum- and angular-momentum-conserving. A finite-difference
approximation has errors controlled by step size and floating
roundoff. This oracle still excludes gravity, caps/fallback transition
regularity proofs and any scalable closed-form derivatives. Its
purpose is to quantify the missing terms, guide their analytic
derivation and reject spurious energy-conservation claims.

```bash
cargo test --lib sph_rk_variational -- --nocapture

cargo run --release --bin phy21_variational_bench -- \
  --n 16 --seed 17 --layout shell --eta 1.7 \
  --fd-fraction 0.0001 --steps 2 --dt 0.002 \
  --csv bg-phy21b-n16.csv
```

The [N16/N24 variational matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy21b.yml)
is isolated from the fast N32/N64/N128 compatible-work workflow.
It is explicitly a research-grade O(N coordinate probes × full
density recomputation) reference, not a playable solver.

## Measured BG-PHY-21B total-variation and thermodynamic controls

The [BG-PHY-21B N16/N24 matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42803311423) passed all tests, measuring four kernel/density combinations at each resolution. The complete reference energy uses the existing Fermi EOS and recomputes all particle dependent h/rho/MLS/partition terms on every coordinate perturbation.

| Full-energy oracle diagnostic | Worst of eight tested cases |
|---|---:|
| Central-difference step-halving relative force RMS | 2.03e-8 |
| Unprojected normalized total force | 1.17e-11 |
| Unprojected normalized torque | 1.80e-9 |
| Pressure-only Hamiltonian relative energy drift | 6.19e-12 |
| Relative momentum drift during variational steps | 1.09e-14 |

**The conventional SPH force validates the total-energy gradient.**
When the same exact-root *original SPH density* is used for the Fermi
internal energy and for the ordinary adaptive-h SPH pressure formula,
their relative force RMS discrepancy is 8.7e-9 (N16/M4), 1.17e-8
(N16/C4), 8.6e-9 (N24/M4), and 1.84e-8 (N24/C4).
The regression test now enforces agreement below 1e-6.

**Measured force discrepancy from the full cold-EOS gradient:**

| N | Kernel | Density | Fast CRK pressure RMS discrepancy | Original SPH formula RMS discrepancy |
|---:|---|---|---:|---:|
| 16 | M4 | original | 1.5612 | 8.73e-9 |
| 16 | M4 | surface blended | 1.6033 | 0.04056 |
| 16 | C4 | original | 1.0987 | 1.17e-8 |
| 16 | C4 | surface blended | 1.0912 | 0.03869 |
| 24 | M4 | original | 1.3456 | 8.60e-9 |
| 24 | M4 | surface blended | 1.3877 | 0.06843 |
| 24 | C4 | original | 0.4901 | 1.84e-8 |
| 24 | C4 | surface blended | 0.4507 | 0.06734 |

For the surface-blended density, the old symmetric SPH pressure
formula is being evaluated at its *original* SPH rho and does not
contain the blend's extra position/h/MLS derivatives; the 4–7%
discrepancy measures the resulting missing thermodynamic terms
relative to the full blended-U gradient. The fast reproducing-kernel
pressure force is substantially farther from the exact discrete
energy gradient at these resolutions, despite its exact momentum
and compatible-work account.

The strong separation between the conventional SPH control
(~1e-8 force mismatch) and CRK pair force (~0.45–1.6 force
mismatch) validates the measurement method and gives a concrete
target for an analytic conservative derivation. A CRK force can
use compatible work to preserve K+evolved-U exactly yet fail to
correspond to the cold-EOS Hamiltonian U(rho_hat). The energy
identity in isolation is insufficient to certify the new
pressure law as a physically faithful white-dwarf model.

**Next gate:** derive the full analytic derivative of the
surface-gated rho_hat via
`d rho_SPH/dx`, `d h/dx` from the implicit closure,
`d S/dx`, `d beta_MLS/dx` from the affine normal matrix,
and their cross-particle couplings. Compare each term with the
BG-PHY-21B complete finite-difference reference before designing
a scalable pairwise/compatible energy implementation. The existing
ordinary SPH formula remains the exact barotropic control.

