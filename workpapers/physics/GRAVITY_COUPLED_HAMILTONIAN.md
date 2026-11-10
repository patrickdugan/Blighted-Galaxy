# BG-PHY-23: Poisson-matched gravity with the analytic MLS pressure Hamiltonian

**Date:** November 10, 2026  
**Predecessor:** [BG-PHY-22](ANALYTIC_SURFACE_VARIATION.md)  
**Tracking:** [Issue #25](https://github.com/patrickdugan/Blighted-Galaxy/issues/25), [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Mathematical objective

The previous BG-PHY-22 proof-by-numerical-oracle established that the
analytic force of a cold-Fermi internal-energy potential with
surface-gated MLS density agrees with independent finite differences
at relative RMS below 4e-8. It also preserves exact momentum and
angular-momentum symmetries in a pressure-only Hamiltonian.

The present step couples that *unchanged* pressure force to exactly
the Poisson-matched, variational variable-h self-gravity already in
`sph_adaptive`. We define

```text
H(x,v) = T(v) + U_cold(x) + Phi_G(x),

T = sum_i m_i |v_i|² / 2,

U_cold = sum_i m_i u_Fermi(rho_hat_i(x)),

Phi_G = G/2 sum_{i<j} m_i m_j [
  phi_K(r_ij,h_i(x)) + phi_K(r_ij,h_j(x))
].
```

The pressure density is either the exact original SPH density or
BG-PHY-20's partition-gated affine MLS rho_hat. In *both* variants,
the smoothing length h_i and gravity kernel are determined from the
**original** exact-root SPH rho_i and its `Omega_i`, preserving the
matched gravity potential and established closure.

The pair potential `phi_K` is the Poisson integral of the selected
M4 or Wendland C4 density kernel, is negative, and equals `-1/r`
beyond 2h. Each half-pair uses the corresponding particle's
smoothing length, precisely as the existing `energy_j_cached`
and matched gravity acceleration.

## Force derivation and implementation

The pressure force is BG-PHY-22's independently audited

```text
F_i^P = -partial U_cold/partial x_i,
```

including the implicit h root, source volumes, MLS moment inverse,
the partition `S_i`, smoothstep gate, and all cross-particle
coordinate dependence.

The gravitational force has the ordinary pair component

```text
F_i,ij^G,pair = G m_i m_j * 0.5[
  phi'_r(r_ij,h_i)+phi'_r(r_ij,h_j)
] * (x_j-x_i)/r_ij
```

and the implicit-h component from the potential's h dependence:

```text
zeta_i = -(h_i/(3 rho_i)) sum_(j ne i)
          m_j partial_h phi_K(r_ij,h_i),

Omega_i = 1 + (h_i/(3 rho_i)) partial_h rho_i,

F_i,ij^G,h = -G m_i m_j / 2 *
  [ (zeta_i/Omega_i) grad_i W_ij(h_i)
  + (zeta_j/Omega_j) grad_i W_ij(h_j) ].
```

This is the existing `sph_adaptive::forces().gravity_accelerations_m_s2`,
so the new research pilot reuses its measured and tested implementation
rather than introducing a second gravity discretization. The code
checks the combined force against that original exact-root SPH
force when the density correction is disabled. The expected
agreement is at floating-point precision.

The new `sph_surface_gravity` module exposes:
- a coupled force/energy snapshot with separate pressure, gravity,
  cold internal and gravitational potential components;
- an independent total-potential finite-difference oracle that
  **re-solves** smoothing, density, partition and MLS for each
  coordinate perturbation;
- an opt-in gravity-coupled Hamiltonian `kick-drift-kick` integrator,
  plus a separate damped preparation operation.

This reuse provides exact conservation by differentiating one
time-independent scalar Hamiltonian, subject to finite integration
step error and piecewise-active h/MLS branches.

## Independent numerical verification

At N16/24/32/64, M4 and C4 are each crossed with original and
surface-blended cold EOS density, with uniform mass seeds and
a dedicated unequal-mass asymmetric seed. Each case compares

```text
F_analytic_i = F_P,i + F_G,i

F_fd,i = -[H_potential(x+eps e_i)-H_potential(x-eps e_i)]
         / (2eps).
```

The scalar FD path uses BG-PHY-21B's independently defined cold-EOS
potential and a **fresh** exact h/rho re-solve for the gravitational
potential. That makes the gradient check independent of the analytic
pressure force, and uses the exact same scalar as the coupled
Hamiltonian. The root and kernel definitions themselves are shared
by scientific necessity.

The numerical gate is relative RMS force discrepancy below 2e-5,
normalized unprojected net force and torque below 1e-9, and
short 4-step Hamiltonian energy drift below 2e-5 with reversible
positions/velocities. The original-density control must reproduce
the existing matched SPH total acceleration within 1e-10 relative
force RMS. All variants report the relative gravitational potential,
cold-EOS potential, h caps, minimum Omega, runtime and FD step
fraction.

N128 is tested without the expensive full coordinate-difference
oracle; it checks analytic pressure-plus-gravity evolution and
momentum/energy/torque invariants.

## Physical acceptance pilot

The [BG-PHY-23 hold workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy23-hold.yml)
runs matched N128 and N256 pilot systems, each with:
- the same 0.65-solar-mass, mu_e=2 independent hydrostatic reference
  and angular-glass initializer;
- 3.00s constant-rate damping at tau=0.7s, dt=0.01s;
- at least 7.88s (3 t_dyn) *undamped* Hamiltonian evolution;
- original matched M4/eta=1.20 SPH density control, and
  matched Wendland C4/eta=1.70 surface-gated MLS candidate.

The acceptance is independent of CI success:

```text
abs(post_Rmax/Rref-1) < 0.05
max_undamped_abs(Rmax/Rmax_post-1) < 0.05
max_undamped_abs(R99/R99_post-1) < 0.05
max_radial_mass_quantile_RMS < 0.05

energy_relative_drift < 1e-3
momentum_relative_drift < 1e-8
angular_relative_drift < 1e-8
undamped_duration >= 3 * sqrt(Rref^3/(G Mref)).
```

The runner writes a summary CSV with physical-gate Boolean, a
longitudinal 32+snapshot trace, separate EOS/grav energies and
the exact seed's pre-existing sampling-radius offset. A numerical
success with a physical gate failure is an informative result;
the physical thresholds are not tuned to the outcome.

## Reproduction

```sh
cargo test --lib sph_surface_gravity -- --nocapture
cargo run --release --bin phy23_bench -- \
  --n 24 --seed 0 --layout glass --eta 1.7 \
  --steps 4 --dt 0.002 --fd-fraction 0.0001 \
  --csv bg-phy23-n24.csv

cargo run --release --bin phy23_hold -- \
  --n 256 --seed 0 --layout glass --kernel c4 --density blend \
  --eta 1.7 --relax-steps 300 --hold-steps 788 \
  --dt 0.01 --damping-tau 0.7 \
  --csv bg-phy23-n256-hold.csv \
  --trace-csv bg-phy23-n256-hold-trace.csv
```

## Measured results: total energy and four physical holds

The [final dedicated full-potential matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42825830058) passed N16, N24, asymmetric/unequal-mass N32, N64 and analytic-only N128. Each of the first four resolutions evaluated original SPH density and surface-gated MLS with the same M4/Wendland matched Poisson gravitational potential. The direct coupled force reproduces the existing exact-root SPH force in the original-density control to floating-point precision.

| Numerical check | Worst observed |
|---|---:|
| Analytic vs independently re-solved total-potential force, N16–64 | **5.104e−7** relative RMS |
| Original-density coupled force vs existing adaptive SPH force | **4.635e−15** relative RMS |
| Analytic net force ratio, unprojected | 1.04e−15 |
| Analytic net torque ratio, unprojected | 1.14e−15 |
| Short Hamiltonian relative energy drift | 4.59e−12 |
| Independent total potential scalar discrepancy | 0 at CSV precision |

The new force is the gradient of the full *cold EOS internal plus
self-gravitational energy*, including MLS/partition and variable-h
terms. Numerical success is independent of the physical radius and
stable-hold acceptance.

The [four full-duration stellar holds](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42825829940) completed successfully after raising the BG-PHY-22 research-only particle guard from N128 to N256. Every hold includes 3.00 seconds of tau=0.7s damping and 7.88 seconds of undamped Hamiltonian evolution.

| Metric | N128 M4 original | N128 C4 blended | N256 M4 original | N256 C4 blended |
|---|---:|---:|---:|---:|
| Post-damping Rmax (km) | 6875.95 | 6569.18 | 6973.46 | 6512.16 |
| Rmax bias relative to Rref | −18.16% | −21.82% | −17.00% | −22.49% |
| Contraction from sampled seed radius | −10.24% | −14.24% | −11.06% | −16.94% |
| Maximum undamped Rmax drift | 14.84% | 10.85% | 10.95% | **5.43%** |
| Maximum undamped R99 drift | 14.84% | 10.85% | 9.36% | **4.62%** |
| Worst radial mass-quantile RMS | 6.06% | 6.36% | 4.35% | **3.81%** |
| Maximum total Hamiltonian energy drift | 2.21e−8 | 1.23e−7 | 1.56e−8 | 6.06e−8 |
| Maximum normalized momentum drift | 1.09e−16 | 1.96e−16 | 1.43e−16 | 2.12e−16 |
| Maximum normalized angular-momentum drift | 8.12e−16 | 2.39e−15 | 3.76e−16 | 1.16e−15 |
| Full duration checked | yes | yes | yes | yes |
| **Issue #8 physical acceptance** | **FAIL** | **FAIL** | **FAIL** | **FAIL** |

All four simulations pass the numerical conservation checks over
the full duration. The C4 surface blend produces substantially
smaller radial oscillations at N256, with R99 excursion and
radial mass-CDF RMS individually meeting the original 5% thresholds.
Its outermost-particle drift remains slightly above 5%, and its
post-damping stellar radius is 22.49% below the independent 1D
reference.

### Matched-kernel density ablation at N256

The first physical comparison varied the density estimator and kernel
together. A second [six-job hold matrix](https://github.com/patrickdugan/Blighted-Galaxy/actions/runs/42827024133)
isolated the effect of the BG-PHY-20 surface correction **at fixed
kernel, eta, particle seed and timestep**. All six jobs completed
their original 3.00s damping / 7.88s undamped protocol; the
two N128 controls are retained. Every N256 run has identical
initial glass particle positions and cold EOS.

| N256 metric | M4 original, eta1.20 | M4 surface MLS, eta1.20 | C4 original, eta1.70 | C4 surface MLS, eta1.70 |
|---|---:|---:|---:|---:|
| Post-damped radius bias | −17.00% | −20.13% | −17.26% | −22.49% |
| Maximum Rmax drift | 10.95% | **5.34%** | 15.10% | **5.43%** |
| Maximum R99 drift | 9.36% | **4.68%** | 13.22% | **4.62%** |
| Worst radial CDF RMS | 4.35% | **3.20%** | 5.29% | **3.81%** |
| Maximum relative Hamiltonian energy drift | 1.56e−8 | 6.70e−8 | 1.21e−8 | 6.06e−8 |
| Relative angular-momentum drift | 3.76e−16 | 7.34e−16 | 2.98e−16 | 1.16e−15 |
| Issue #8 physical gate | FAIL | FAIL | FAIL | FAIL |

**Controlled result:** at fixed M4 the density correction reduces
maximum undamped Rmax excursions by approximately 51%, while at fixed
Wendland C4 it reduces them by approximately 64%. Both kernels gain
R99 stability and improve the radial mass-CDF error.

The same correction shifts the post-damping equilibrium inward
by approximately 3.13 percentage points of reference radius under
M4 and 5.23 percentage points under C4. These are numerical
outcomes of the complete Hamiltonian, already verified against
independent energy differentiation, rather than consequences
of a force or EOS normalization change.

All six physical acceptance Booleans remain **false**.
For both corrected N256 systems, R99 and CDF pass their 5%
individual gates; maximum Rmax remains narrowly above 5%
(5.34% for M4; 5.43% for C4), while the independent
hydrostatic radius remains 20–22.5% too small.

The density correction is therefore the causal source of
improved stability **within the tested kernel and eta values**.
The larger equilibrium contraction remains the next blocking
hydrostatic discrepancy. A radial virial/outer-shell force
decomposition at common frozen positions would isolate whether
this is surface support truncation or a broader pressure-density
functional bias.

### Scientific implication

The analytic variational correction now makes it possible to
compare different density reconstructions without losing the
cold-barotropic energy contract. It does **not yet** recover the
correct stellar equilibrium. The better N256 dynamical stability
of the C4 surface blend coexists with stronger radial contraction,
suggesting that the main remaining defect involves the equilibrium
surface/pressure-support model, free-boundary treatment, or
particle discretization at fixed N.

This experiment keeps the cold EOS and the reference hydrostatic
integrator unchanged. Any next intervention should be evaluated
against the *independent radius* and outer-shell force support,
while preserving the successful Hamiltonian conservation and the
existing Issue #8 acceptance gates.

A practical next discriminating step is a controlled free-surface
support closure with boundary-aware volume and gradient moments,
testing whether the radial force balance at the outer 1–2
mass shells recovers the reference at N256 before a large N512/1024
dynamics campaign. Retain the exact EOS energy potential as the
derivation and conservation target.

## Limitations and next gate

The solver is research-only. It currently has no compact binary
dynamical coupling to the game engine, no central black hole,
radiative transfer, nuclear burning, hydrodynamic shocks,
general relativity, or aspherical tidal destruction validation.
The exact MLS correction has piecewise branches where
h bounds become active or the affine moment solver falls back.
The coupled Hamiltonian is correct within a smooth active region;
detect and report switching when validating large movements.

For gameplay and white-dwarf attack sequences, test rotationally
asymmetric initial conditions, larger SPH resolutions, and
time-step convergence before enabling any default behavior.
All original physics paths remain unchanged. Issue #8 remains
the independent hydrostatic acceptance source.

### Research sources

- [BG-PHY-22 analytic surface density and EOS variation](ANALYTIC_SURFACE_VARIATION.md).
- [BG-PHY-21 total discrete energy-gradient oracle](CONSERVATIVE_RK_PILOT.md).
- [BG-PHY-20 density and enthalpy diagnostics](FREE_SURFACE_DENSITY.md).
- Price & Monaghan (2007), *An energy-conserving formalism for adaptive
  gravitational force softening*, DOI 10.1111/j.1365-2966.2006.11241.x.
- Frontiere, Raskin & Owen (2017), *CRKSPH*, DOI
  10.1016/j.jcp.2016.12.004.
