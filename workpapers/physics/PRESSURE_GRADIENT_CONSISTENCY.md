# BG-PHY-18: Frozen gradient consistency and free-surface pressure bias

**Date:** 10 November 2026  
**Dependencies:** [BG-PHY-17](HIGH_NEIGHBOR_C4_SWEEP.md), [BG-PHY-16](WENDLAND_GLASS_SEEDS.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## 1. Experimental reason

BG-PHY-17 has completed a 100–200-neighbor Wendland C4 convergence sweep.
At N1024 with angular-glass seed 0, raising support from 100 to 200
neighbors lowers the initial normalized pressure-vector RMS from 0.329
to 0.120, and reduces inner E0 noise. Simultaneously, the outermost
equal-mass shell pressure/reference radial ratio falls from 0.830 to
0.750. The complete physical hold measured:

| N1024 glass seed 0 | M4 eta 1.20 | C4 150 neighbors | C4 200 neighbors |
|---|---:|---:|---:|
| post-damping absolute Rmax bias | -13.68% | -14.26% | -14.66% |
| worst undamped Rmax excursion | 8.10% | 10.03% | 11.80% |
| worst undamped R99 excursion | 3.45% | 5.43% | 7.10% |

The simulation conserves numerical energy and momentum in all cases,
but fails the independent hydrostatic radius and stability gates.

BG-PHY-18 tests a narrower hypothesis: **the free boundary breaks
discrete first-order pressure-gradient consistency while the SPH
density estimator introduces a separate field bias**. The experiment
runs on frozen seeded stars. It never substitutes a non-variational
gradient into the momentum or energy equations.

## 2. A first-moment correction for diagnostics only

Write V_j=m_j/rho_j and q_ij=x_j-x_i. For the actual SPH density
kernel and current smoothing length, define

```text
M_i = sum_(j != i) V_j grad_i W_ij(h_i) ⊗ q_ij.
D_i(f) = sum_(j != i) V_j [f_j-f_i] grad_i W_ij(h_i).
```

The direct difference gradient D_i exactly differentiates a constant
field to zero, but a linear field f=f0+b·x yields D_i=M_i b.
For a nonsingular moment matrix, the **linear-consistent**
reconstruction is therefore

```text
grad_i f = M_i^(-1) D_i(f).
```

The implementation analytically inverts the symmetric 3x3 moment
matrix after checking positive determinant, a trace-scaled cutoff and
condition number <=1e7. Missing or singular matrices are excluded from
the paired error metrics, and the usable fraction is recorded. Direct
and corrected gradients share exactly the same valid-particle subset.

The corrected pressure acceleration is evaluated as a **frozen
diagnostic**

```text
a_P,direct(i) = -D_i(P_num) / rho_num,i
a_P,corr(i)   = -M_i^(-1) D_i(P_num) / rho_num,i.
```

This directed estimator is not a conservative pair-force law. Its use
for live dynamics would require new pair-symmetry and energy/variational
derivations, plus independent momentum and angular-momentum tests.
Neither the existing gameplay default nor the research integrator is
changed in this milestone.

## 3. An independent pressure-field oracle

We also evaluate the same D_i and M_i using pressure and density
sampled from the one-dimensional hydrostatic reference:

```text
rho_ref,i = rho_1D(|x_i-C|)
P_ref,i   = P_Fermi(rho_ref,i)
a_P,ref,direct(i) = -D_i(P_ref) / rho_ref,i
a_P,ref,corr(i)   = -M_i^(-1) D_i(P_ref) / rho_ref,i.
```

The volume factors V_j remain the *same SPH m_j/rho_num,j* in
both numerical and reference oracle estimators; the only changes are
the sampled scalar field and the local denominator. This separates
errors in the pressure samples from the geometry's ability to resolve
the reference field at fixed spatial coordinates.

The independent expected outward hydrostatic acceleration is

```text
a_ref,i = [G M_1D(<r_i)/r_i^2] * rhat_i.
```

Each vector error uses the same reference normalization

```text
epsilon(a) =
sqrt(sum_i m_i |a_i-a_ref,i|^2 /
     sum_i m_i |a_ref,i|^2).
```

In addition to four epsilon values, the benchmark records signed
radial pressure/reference ratios, tangential RMS, the moment trace
and condition number, as well as the log-density discrepancy
ln(rho_num/rho_ref). Both global and inner/outer enclosed-mass
quartile values are emitted.

The oracle estimates interpolation/gradient error on the actual
irregular point cloud using the correct continuum pressure samples.
The numerical-field estimates retain the kernel-density sampling
error. A large difference between oracle and numerical-field errors
would implicate the SPH density/pressure field. If the oracle and
numerical field both fail near the surface, the kernel moment
truncation and nonlinear free-boundary gradient remain primary
suspects.

## 4. Frozen factorial and guardrails

The CI matrix uses N256, N1024 and N2048, with independent shell
orientations seed0/17 where affordable, and seed0 at N2048.
At each N, angular-shell and angular-glass coordinates are frozen
across kernels. M4 eta=1.20 is the standard baseline; Wendland C4
is compared at eta=1.20 and exact-root-calibrated 100,150,200
physical support neighbors.

The calibration and all force measurements use the safeguarded
Newton density root from BG-PHY-14. The new matrix correction
depends only on kernel gradients at frozen particle coordinates.
The existing particle initializer, cold-electron EOS, gravity
constants, conserved pair forces, variable-h corrections and
star-acceptance gates remain unchanged.

## 5. Reproduction

```sh
cargo test --all-targets

cargo run --release --bin phy18_bench -- \
  --n 1024 --seeds 0,17 --targets 100,150,200 \
  --glass-sweeps 20 --calibration-steps 12 \
  --csv bg-phy18-1024.csv

# Fast verification smoke
cargo run --release --bin phy18_bench -- \
  --n 64 --seeds 0 --targets 20,32 \
  --glass-sweeps 4 --calibration-steps 6 \
  --csv bg-phy18-smoke.csv
```

See [BG-PHY-18 GitHub Actions](https://github.com/patrickdugan/Blighted-Galaxy/actions/workflows/bg-phy18.yml).
A full matrix CSV allows paired comparisons by N, seed, layout, kernel
and target, and includes the baseline uncorrected conservative
SPH pressure-force error for context.

## 6. Scientific decision

If the oracle-gradient correction restores the outer-shell pressure
ratio to approximately 1 and sharply reduces its RMS error, the
dominant defect in that diagnostic is the discrete first moment.
If the numerical-pressure corrected gradient remains much worse,
inspect density reconstruction and kernel normalization near the
vacuum interface. Large density log-error and nonunit moment trace
may coexist.

The next candidate for an actual dynamical implementation must obey
conservative flux/variational constraints and demonstrate stable
momentum, energy, and angular momentum before it is considered
against the original M4 force law. BG-PHY-18 makes no gameplay
switch and keeps Issue #8 open.

## References

- J. Bonet and T.-S. L. Lok (1999), *Variational and momentum
  preservation aspects of Smooth Particle Hydrodynamic
  formulations*, Computer Methods in Applied Mechanics and
  Engineering 180, 97–115,
  https://doi.org/10.1016/S0045-7825(99)00051-1.
- D. García-Senz, R. M. Cabezón and J. A. Escartín (2012),
  *Improving smoothed particle hydrodynamics with an integral
  approach to calculating gradients*, Astronomy & Astrophysics
  538, A9, https://doi.org/10.1051/0004-6361/201117939.
- W. Dehnen and H. Aly (2012), *Improving convergence in smoothed
  particle hydrodynamics simulations without pairing instability*,
  MNRAS 425, 1068, https://doi.org/10.1111/j.1365-2966.2012.21439.x.
