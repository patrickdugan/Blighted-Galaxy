# BG-PHY-14 — Neighbor growth, implicit h–ρ roots, and white-dwarf seed ensembles

**Date:** 10 November 2026  
**Source:** Blighted Galaxy research reconstruction  
**Predecessor:** [BG-PHY-13](ETA_HYDROSTATIC_SWEEP.md)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)  
**Code:** \`src/sph_hroot.rs\`, \`src/sph_adaptive.rs\`, \`src/sph_sampling.rs\`, \`src/bin/phy14_bench.rs\`.

## 1. Research question

BG-PHY-13 compared fixed values of the adaptive SPH smoothing parameter eta at N=256,512,1024. Every test used the same mu_e=2 cold-electron equation of state, the Poisson-matched cubic M4 gravitational kernel and its grad-h correction, the same 3-second damped preparation at each count, and a 7.88-second undamped hold (3.005 reference dynamical times).

All tested stars failed at least one independent equilibrium gate. Eta=0.85 reduced the radius deficit at N1024 to 8.78%, with a worst outer-radius excursion of 5.77%, and had only about 21 neighbors inside the 2h kernel support. Holding eta fixed as N increases keeps the *physical* neighbor count nearly constant.

This milestone tests three distinct hypotheses:

1. **Neighbor count growth:** a resolution-scaled eta(N) can improve SPH consistency, with the physical neighbor count increasing as N grows.
2. **Implicit h–rho closure:** finite underrelaxed smoothing-length iterations may leave a measurable residual in h_i^3 rho_i = eta^3 m_i, affecting the recovered density and forces.
3. **Angular seed sensitivity:** independent shell orientation realizations may change low-resolution SPH force noise and stability even when radial mass quantiles remain identical.

The experiment retains the cold electron EOS and does not assume any hypothesis is true before seeing the benchmark results.

## 2. Scientific basis of the neighbor-growth law

Zhu, Hernquist & Li (2015), [*Numerical Convergence in Smoothed Particle Hydrodynamics*](https://doi.org/10.1088/0004-637X/800/1/6), derive conditions for formal SPH consistency:

    N_total → infinity,
    h → 0,
    N_neighbors → infinity.

For typical density kernels they identify the compromise scaling N_neighbors proportional to N_total^(1/2).

For the 3D M4 kernel with support r<2h and equal-mass particles,

    h_i = eta * (m_i/rho_i)^(1/3),

and a uniform-medium estimate gives

    N_neighbors ≈ (32*pi/3) eta³.

Therefore the compromise scaling becomes

    eta(N) = eta_ref * (N/N_ref)^(1/6).

The benchmark uses eta_ref=1.05 at N_ref=256 as an initial test, paired at every resolution with **fixed eta=1.20**. The resulting nominal growth schedule is:

| N | Fixed eta | Growth eta | Approximate growth-law support count |
|---:|---:|---:|---:|
| 256 | 1.2000 | 1.0500 | 38.8 |
| 512 | 1.2000 | 1.1786 | 54.9 |
| 1024 | 1.2000 | 1.3229 | 77.6 |

The support counts are uniform-density estimates. The benchmark separately measures actual neighbors within each particle's 2h_i support. Strong density stratification, kernel caps and free-surface particles can cause departures.

The chosen eta_ref is an exploratory initial scale, not an astrophysically calibrated value. The M4 kernel also has neighbor-number limits from particle pairing and E0 force noise. A nominally correct growth exponent alone does not establish a convergent or stable white dwarf.

## 3. An implicit smoothing-length root solver

Reference: Price & Monaghan (2007), [*An energy-conserving formalism for adaptive gravitational force softening*](https://doi.org/10.1111/j.1365-2966.2006.11241.x); Price (2008), [*Modelling discontinuities and Kelvin–Helmholtz instabilities in SPH*](https://users.monash.edu/~dprice/pubs/kh/kh.pdf).

BG-PHY-08 used a fixed number of underrelaxed iterations to update

    rho_i(h_i) = sum_j m_j W(r_ij,h_i),
    h_i = eta*(m_i/rho_i)^(1/3).

BG-PHY-14 adds a separate safeguarded Newton–bisection root, solving the dimensionless equation

    F_i(h) = h³ rho_i(h)/(eta³ m_i) - 1 = 0.

Its derivative is

    F'_i(h) = [3h² rho_i(h)+h³ d(rho_i)/dh] / (eta³ m_i)
            = 3h² rho_i(h) Omega_i / (eta³ m_i),

where Omega_i = 1 + h/(3rho) d(rho_i)/dh is the same variable-h correction used by the momentum-conserving pressure and gravitational grad-h forces.

The code maintains a bracket over the configured smoothing-length interval and uses a Newton step only when its derivative is finite, positive, and the proposed next h lies strictly inside the bracket. Otherwise it bisects. Residuals are measured in a dimensionless form; **roots outside the permitted h interval are reported as capped solutions** instead of being silently counted as converged roots.

A root result records h, rho, Omega, d(rho)/dh, residual, iteration count and cap/interior classification. Root candidates include the complete permitted maximum support, including configurations in which the root lies at an h bound.

The original underrelaxed method remains the default in the playable engine. An opt-in \`with_newton_density()\` mode allows numerical comparisons. At identical frozen particle positions, the benchmark audits the old closure residual, computes the safeguarded root and measures the RMS difference in SPH accelerations when the root solver is used. This distinguishes smoothing-length integration error from a physically incorrect initial stellar profile.

### Precision and cost limitations

The safeguarded root may take many more density evaluations than the original fixed 24-pass underrelaxed method, particularly with low neighbor counts. The root dynamics therefore receive a small-N CI pilot first. The high-resolution comparison keeps the established integration scheme, using the root solver for **snapshot audits only**. The code records residuals at capped and uncapped particles separately. An Omega floor and h cap can still affect exact variational dynamics.

A Newton root of the discrete equation guarantees **closure consistency at the represented particle positions**. It does not guarantee that the discrete density matches the continuum, that the numerical particle distribution is isotropic, or that pressure balances gravity at the correct radius.

## 4. Independent orientation ensembles with fixed radial mass CDF

The new \`seed_cubic_shells_oriented\` implementation starts from the approximately cubic-shell layout introduced in BG-PHY-10. Each radial shell receives an independent, deterministically generated three-dimensional proper rotation using a seeded splitmix64 stream and Euler angles. Each radial pair is rotated together with its antipode.

That construction preserves the hydrostatic radial mass quantile for every parcel, particle mass, initial bulk velocity and center of mass. Its within-shell angular arrangement is unchanged up to a rotation; the *relative orientations of different shells* vary. A single global rigid rotation would be physically redundant in an isotropic self-gravitating star, so shell-specific rotations are required to expose numerical sensitivity.

The baseline seed \`0\` reproduces the prior particle layout. Seed \`17\` supplies the independent layout in N256 and N512 validation. The first N1024 matrix is restricted to seed 0 to keep the exact all-pairs gravitational force comparison within a practical CI budget; a second N1024 seed is an explicit follow-up if the low-N ensemble shows significant spread.

## 5. Matched-duration experiment

The new \`phy14_bench\` compares fixed eta=1.20 with the growth eta(N) law, keeping the EOS, M4 gravity, constant tau=0.7s damping and physical time interval identical between cases at each resolution.

| N | Orientation seeds | Damped preparation | Undamped hold | Dynamic density method |
|---:|---|---:|---:|---|
| 256 | 0, 17 | 3.00s | 7.88s | Existing underrelaxed; Newton audit |
| 512 | 0, 17 | 3.00s | 7.88s | Existing underrelaxed; Newton audit |
| 1024 | 0 | 3.00s | 7.88s | Existing underrelaxed; Newton audit |

Additional smoke tests run N64 using the original dynamic density method and a very small N32 opt-in Newton-dynamics simulation. These have shorter holds and are explicitly disqualified from passing the physical equilibrium gate.

Each case emits:

- realized support-neighbor count and post-relaxation hydrostatic radius bias;
- post-relaxation radial mass-CDF RMS, independent hydrostatic gravity-field RMS discrepancy and global force-residual ratio;
- worst undamped Rmax and R99 excursions across intermediate snapshots, worst CDF mismatch, final energy/momentum drift;
- maximum/RMS residual of the original h–rho closure, residual of the safeguarded Newton solution, cap fraction and maximum change in h;
- relative acceleration difference between the two density closures at exactly the same particle positions.

The benchmark scores the physical gate independently:

    complete_hold = duration >= 3*t_dyn
    radius_ok = |R_post/R_hydro-1| < 0.05
    stable = max|Rmax(t)/Rmax(0)-1| < 0.05
             AND max|R99(t)/R99(0)-1| < 0.05
    profile_ok = max_t(CDF_RMS(t)) < 0.05
    pass = complete_hold AND radius_ok AND stable AND profile_ok.

The mass/momentum invariants and bounded integration-energy error remain software assertions. No artificial radius adjustment or change in the physical electron pressure coefficient is permitted.

## 6. Reproduction

    cargo test --all-targets

    cargo run --release --bin phy14_bench -- \
      --n 512 --seeds 0,17 \
      --fixed-eta 1.20 --growth-base-eta 1.05 \
      --growth-base-n 256 --relax-steps 300 --hold-steps 788 \
      --csv bg-phy14-n512.csv

    cargo run --release --bin phy14_bench -- \
      --n 32 --seeds 0 --root-dynamics \
      --relax-steps 4 --hold-steps 8

The dedicated **BG-PHY-14 growing-neighbor consistency** Actions workflow runs one job per particle count, and publishes numeric CSVs. It runs N256 and N512 with two seeds; N1024 with a baseline seed. Regular CI verifies all new Rust tests, both short smoke programs, and Windows release compilation.

**Scientific decision:** BG-PHY-13 already showed that reducing eta can produce a larger radius while worsening the local force residual. BG-PHY-14 will permit a decision only after the measured neighbor counts, root closure residual, independent angular-seed variance, and long undamped stability are known. [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) remains open until both the original radius and undamped-stability thresholds pass.
