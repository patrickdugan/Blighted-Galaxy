# BG-PHY-22: Analytic derivative of surface-gated MLS cold-EOS energy

**Date:** November 10, 2026  
**Predecessor:** [BG-PHY-21](CONSERVATIVE_RK_PILOT.md)  
**Tracking:** [Issue #22](https://github.com/patrickdugan/Blighted-Galaxy/issues/22)  
**Physical equilibrium:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## Aim and exact energy

BG-PHY-21B evaluated a complete cold-barotropic force through expensive
coordinate finite differences. At N=16/24, the existing adaptive-h SPH
pressure force matches that exact energy derivative to approximately
1e-8 *when using original SPH density*. The BG-PHY-20 surface-gated
density introduces additional implicit coordinate dependencies. Simply
applying the old pressure force to that density missed 4–7% of the
variational force in those small-N cases.

BG-PHY-22 implements the total analytic derivative of the
specific BG-PHY-20 surface-gated density, giving

```text
U(x) = sum_i m_i u_cold(rho_hat_i(x))
F_k = -dU/dx_k
lambda_i = dU/d rho_hat_i
         = m_i P_cold(rho_hat_i)/rho_hat_i².
```

The EOS, 0.65 solar mass hydrostatic reference, particle masses and
matched density/gravity kernel definitions remain unchanged.
The original gameplay integrator is untouched.

## SPH baseline and implicit smoothing-length derivative

The baseline is the exact safeguarded root of

```text
rho_i(h_i,x) = sum_j m_j W(|x_i-x_j|,h_i)
h_i³ rho_i = eta³ m_i.

Omega_i = 1 + h_i/(3 rho_i) * partial_h rho_i.
```

At fixed h, define the vector derivative

```text
partial_(x_k) rho_i |_h
  = sum_j m_j grad_i W_ij (delta_ik - delta_jk).
```

For an uncapped root,

```text
dh_i = -h_i/(3 rho_i Omega_i) * d rho_i|_h
d rho_i = (1/Omega_i) * d rho_i|_h.
```

At a hard smoothing-length cap, dh_i=0 and the density follows
the partial position derivative. The density-gradient floor Omega
and discrete transitions between cap states remain a recognized
nonsmooth branch; the benchmark reports their incidence and does
not treat a root active-set transition as differentiable.

## Surface partition, affine MLS and smoothstep

The BG-PHY-20 field is

```text
V_j = m_j/rho_j
q_ij = (x_j-x_i)/h_i
phi_ij = [1, qx, qy, qz]

w_ij = V_j W_ij
M_i = sum_j w_ij phi_ij phi_ij^T
b_i = sum_j w_ij phi_ij rho_j
beta_i = M_i^(-1) b_i
rho_MLS,i = beta_i[0]
S_i = sum_j w_ij.

t_i = clamp(2(1-S_i),0,1)
alpha_i = t_i²(3-2 t_i)
rho_hat_i=(1-alpha_i)rho_i+alpha_i rho_MLS,i.
```

When the 4×4 MLS matrix is conditioned and its intercept positive,
its exact differential is

```text
d beta_i = M_i^(-1)(d b_i - d M_i beta_i).
```

A Cholesky solver computes beta_i and
`z_i = M_i^(-T)e0` with the same diagonal preconditioning as
BG-PHY-20. Then

```text
d rho_MLS,i = z_i·(d b_i - d M_i beta_i).
```

For each j, put `Z_ij=z_i·phi_ij`,
`B_ij=beta_i·phi_ij`, and `R_ij=rho_j-B_ij`.
The differential becomes

```text
d rho_MLS,i =
sum_j [
  Z_ij R_ij d w_ij
  + w_ij Z_ij d rho_j
  + w_ij (R_ij z_i,grad - Z_ij beta_i,grad) · d q_ij
].
```

The other terms are

```text
d w_ij = V_j d W_ij - (w_ij/rho_j) d rho_j

d q_ij = (d x_j-d x_i)/h_i-q_ij d h_i/h_i

d W_ij = grad_i W_ij·(d x_i-d x_j)
       + (partial_h W_ij) d h_i.
```

In the interior smoothstep region `0<t_i<1`,

```text
d alpha_i = -12 t_i(1-t_i) d S_i,
d S_i = sum_j d w_ij.
```

At the clamped endpoints, the smoothstep derivative vanishes.
The total corrected-density differential is therefore

```text
d rho_hat_i =
 (1-alpha_i) d rho_i + alpha_i d rho_MLS,i
 + (rho_MLS,i-rho_i) d alpha_i.
```

## Reverse accumulation

The analytic implementation accumulates `dU/d rho_i`,
`dU/d h_i` and direct gradients `dU/d x_k` in one N² pass.

For one i,j, the local scalar coefficient of `dw_ij` is

```text
A_ij = lambda_i[
  alpha_i Z_ij R_ij
  + (rho_MLS,i-rho_i) d alpha_i/d S_i
].

D_ij = lambda_i alpha_i w_ij Z_ij
       - A_ij w_ij/rho_j

C_ij = lambda_i alpha_i w_ij
       [ R_ij z_i,grad-Z_ij beta_i,grad].
```

The contributions are

```text
adj_rho_j += D_ij
adj_h_i += A_ij V_j partial_h W_ij
           - (C_ij·q_ij)/h_i

grad_x_i += A_ij V_j grad_i W_ij - C_ij/h_i
grad_x_j -= A_ij V_j grad_i W_ij - C_ij/h_i.
```

Additionally, `adj_rho_i += lambda_i(1-alpha_i)`.

Finally, for each uncapped root j,

```text
K_j = (adj_rho_j - adj_h_j h_j/(3 rho_j))/Omega_j

grad_x_j += K_j sum_k m_k grad_j W_jk
grad_x_k -= K_j m_k grad_j W_jk.
```

At an h cap, the `adj_h` contribution vanishes and
`Omega=1`. The exact analytic force is `F=-grad_x U`.
The code records the RMS magnitudes of its direct-coordinate,
density-chain and implicit-h-chain components separately.

This is a barotropic Hamiltonian pressure formulation for the
specified corrected density field. Its scalar energy construction
makes the exact force translation and rotation invariant; a
finite-difference reference is used to test the implementation,
including the original SPH control.

## Verification gates

The dedicated [BG-PHY-22 workflow](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy22.yml)
tests N16/N24/N32 at both M4 and matched Wendland C4, with
original SPH or surface-blended MLS density. For each case it
recomputes BG-PHY-21B's full energy under +/− coordinate
perturbations and compares the analytic force component-by-component
with that independent finite-difference result.

For the small-N tests, required normalized analytic/FD force RMS
disagreement is below 5e-5. The exact scalar EOS energy must also
agree, and net force/torque symmetry residuals must be below 1e-9.
N64 exercises the analytic formulation without the expensive
N64 full-coordinate finite-difference oracle.

The pilot uses a gravity-free velocity Verlet step with
`F=-grad U(x)` and recomputes U from the EOS for the total
Hamiltonian energy check. Momentum, angular momentum, energy drift
and the reverse-stepping position/velocity residual are measured.
It remains opt-in with no automatic pressure-force replacement.

```bash
cargo test --lib sph_surface_variation -- --nocapture
cargo run --release --bin phy22_bench -- \
  --n 16 --seed 17 --layout shell --eta 1.7 \
  --steps 3 --dt 0.002 --fd-fraction 0.0001 \
  --csv bg-phy22-n16.csv

cargo run --release --bin phy22_bench -- \
  --n 64 --seed 0 --layout glass --eta 1.7 \
  --steps 4 --no-oracle --csv bg-phy22-n64.csv
```

## Remaining distinctions and next use

This milestone differentiates the **full internal energy for the
specific BG-PHY-20 surface-blended density functional**. It uses a
compact kernel with smoothly vanishing support gradients and a
piecewise smooth MLS fallback. Near cap/fallback thresholds,
the validity of a continuous derivative must be checked against
one-sided perturbations and active-set discontinuities.

The pressure-only Hamiltonian test excludes self-gravity, viscosity,
discrete shocks and other stellar gameplay forces. Its analytic
N² accumulation is still a research implementation that must be
profiled before high-resolution adoption. A later integrator must
combine it with Poisson-matched variable-h gravitational potential
derivatives and the cold EOS without introducing an independently
inconsistent thermal reservoir.

The earlier RK-compatible-work force remains a useful
momentum/kinetic-work baseline and a source of RK moment tests.
BG-PHY-22's internal-energy variation gives a separately derived
target for its missing thermodynamic terms.

Once small-N analytic-vs-oracle checks pass, compare N32/64/128
timesteps and support configurations, add gravity with energy
consistency, and only then return to the original Issue #8
white-dwarf physical radius and undamped stability gate.

### References

- Frontiere, Raskin & Owen (2017), *CRKSPH*, JCP 332,
  DOI 10.1016/j.jcp.2016.12.004.
- Owen (2014), *A compatibly differenced total energy conserving
  form of SPH*, DOI 10.1002/fld.3912.
- [BG-PHY-20 density workpaper](FREE_SURFACE_DENSITY.md).
- [BG-PHY-21 total-energy oracle](CONSERVATIVE_RK_PILOT.md).
