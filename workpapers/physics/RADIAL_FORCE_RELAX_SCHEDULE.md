# BG-PHY-12: Radial force balance and tapered white-dwarf relaxation

**Status:** Implemented experimental numerical study, 10 November 2026  
**Prior milestone:** [BG-PHY-11](NEIGHBOR_MASS_PROFILE_HOLD.md)  
**Tracked physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)  
**Code:** src/sph_adaptive.rs, src/stellar_balance.rs, src/stellar_relax.rs, src/bin/phy12_bench.rs.

The earlier BG-PHY-11 1,024-particle cold white dwarf retained an 8.42% maximum excursion of the outer radius during a three-dynamical-time undamped hold. Its 99%-mass radius moved 3.25%, and its post-damping outer radius remained 13.37% below the 8,402 km hydrostatic reference. This motivated a radial breakdown of pressure and gravity and a controlled experiment with a damping schedule.

The original Blighted Galaxy fiction motivates a physically useful white-dwarf disruption engine. The SPH model and diagnostic techniques here are newly developed and remain an experimental approximation to cold degenerate stellar dynamics.

## 1. Reference physics

The independent hydrostatic model integrates a cold, uniform-composition electron-degenerate EOS at mu_e=2 and the Newtonian spherical stellar-structure equations:

    dP/dr = -rho(r)*G*M_enclosed(r)/r²
    dM/dr = 4*pi*r²*rho(r)

Equivalently, the outward pressure acceleration of a hydrostatic fluid element balances its inward self-gravity:

    a_P(r) = -(1/rho)*dP/dr = G*M_enclosed(r)/r².

The numerical star uses the BG-PHY-10 Poisson-matched M4 spline gravitational kernel and its grad-h correction, with cold relativistic electron pressure and fixed-mass particles. The SPH forces are a discrete approximation to the continuum equations. The derived reference pressure/gravity balance is independent of the SPH force solver.

**References:** Pakmor et al., [Stellar GADGET](https://academic.oup.com/mnras/article/424/3/2222/977734), Section 5.2 on long relaxation and time-dependent damping; [Review of numerical common-envelope simulations](https://link.springer.com/article/10.1007/s41115-023-00017-x), on hydrostatic residual and Mach diagnostics; Price & Monaghan (2007), [adaptive gravitational softening](https://doi.org/10.1111/j.1365-2966.2006.11241.x).

## 2. Pairwise force decomposition

The original adaptive force evaluator \`sph_adaptive::forces\` now returns, for every particle i, its exact pairwise **pressure acceleration** and **gravity acceleration**, alongside the existing total acceleration. The pressure component includes the density-dependent Omega grad-h correction; the gravitational component includes the M4 potential and its zeta grad-h correction. The massive central body's field is included in the gravity budget when present, and the isolated-star diagnostics set the central mass to zero.

    a_total_i = a_pressure_i + a_gravity_i

The decomposition leaves the existing total acceleration and integrator unchanged. Unit tests check the reconstruction, and verify that each internal force component preserves net linear momentum within roundoff for all three supported gravitational softening laws.

For a star with center of mass x_COM, define radial unit vector \`e_r=(x_i-x_COM)/r_i\`. The diagnostics use

    a_pressure,r,i = e_r dot a_pressure_i
    a_gravity,in,i = -e_r dot a_gravity_i
    a_net,r,i = e_r dot a_total_i.

A positive a_pressure,r is outward while a positive a_gravity,in is inward. The shell mean uses particle masses as weights; the shell also records **RMS of each particle's net radial acceleration** to expose force fluctuations that a shell mean would otherwise cancel.

The star is divided into radial bins by **enclosed particle mass fraction**. Each bin reports mean radius, particle count, mean pressure and gravity accelerations, continuum reference gravity \`G M_enclosed(r)/r²\`, signed net radial acceleration, RMS radial net acceleration, and RMS tangential net acceleration.

The global diagnostics include

    RMS_net = sqrt[sum_i m_i |a_total_i|² / M]
    RMS_grav = sqrt[sum_i m_i |a_gravity_i|² / M]
    force_residual_ratio = RMS_net / RMS_grav

and a *force-virial ratio*,

    V_P = sum_i m_i (x_i-x_COM) dot a_pressure_i
    V_G = sum_i m_i (x_i-x_COM) dot a_gravity_i
    ratio = V_P / |V_G|.

Hydrostatic equilibrium gives a small total force residual and a force-virial ratio near one. The force virial is a balance diagnostic; the gas additionally has kinetic and internal-energy terms in the full virial theorem. The finite volume, softened gravity and barotropic closure limit quantitative interpretation of this ratio.

## 3. Matched-duration damping schedule test

The BG-PHY-11 reference relaxed the star under a constant damping time \`tau=0.7s\` before disabling the artificial drag for an undamped hold. This milestone preserves that exact constant-damping path.

A new optional schedule uses

    dv_internal/dt = a_SPH - v_internal / tau(t)

with a piecewise damping time:

- Strong initial damping \`tau_initial=0.15s\` until 40% of the relaxation interval.
- Geometric growth of tau to \`tau_final=1.5s\` between 40% and 90%.
- Damping off for the final 10% of the scheduled relaxation, followed by a fully undamped stability hold.

The discrete velocity factor is \`exp(-dt/tau)\` while damping is active, and exactly 1 afterward. Both schedules evolve the **same initial particle coordinates and masses** with the same SPH M4 gravitational law, fixed dt, relaxation duration, and subsequent hold duration. No change is made to the electron EOS, radial reference profile, the equilibrium acceptance gates, or the interaction potential to improve the outcome artificially.

The schedule and parameter validation, exact correspondence of the constant-schedule wrapper with the old implementation, finite tapered relaxation and center/mass preservation have independent Rust unit tests.

The long-run CI workflow uses a **512-particle star**, 600 relaxation steps at dt=0.01s (6s), followed by 788 undamped steps (7.88s or three reference dynamical times), for each schedule. The lightweight CI uses 64 particles and short smoke steps. Both write a time series and radial shell CSVs.

## 4. Reproduction

    cargo test --all-targets
    cargo run --release --bin phy12_bench -- \
      --n 128 --relax-steps 300 --hold-steps 262 \
      --csv bg-phy12-timeseries.csv \
      --shell-csv bg-phy12-shells.csv

CI publishes both in the **bg-phy12-512-reference** artifact from the dedicated \`BG-PHY-12 radial stellar equilibrium A/B\` workflow.

**Time-series fields:** schedule name, phase, coordinate time, maximum radius and R99, radius error relative to hydrostatic, radial-CDF mismatch, pressure/gravity force-virial ratio, RMS net/gravity force ratio, maximum shell radial imbalance, net-force RMS, energy drift and momentum drift.

**Radial shell fields:** schedule and phase (initial, after damping, after undamped hold), shell mass boundaries, count, mean radius, outward pressure acceleration, inward stellar gravity, independent hydrostatic gravity, signed radial residual, radial residual RMS and tangential residual RMS.

The benchmark reports three independent physical criteria: relaxed maximum-radius error within 5% of the hydrostatic star, maximum Rmax/R99 excursion below 5% throughout the undamped hold, and normalized radial-CDF RMS below 5% (an additional research gate). Passing the code's energy/momentum tests cannot substitute for the physical criteria.

## 5. Experimental findings and decisions

The BG-PHY-12 comparison and the next physics-engine choice should be based on measured benchmark results. A smaller RMS force residual may coexist with too-small radius, large tangential force noise or a mismatch in the mass profile. A better schedule should improve the full set of diagnostics, especially the unforced radial and mass-CDF evolution.

The numerical solver still uses exact all-pairs self-gravity O(N²), with an accelerated compact-support density/grad-h neighbor index. A 3D self-gravitating white dwarf with 512 to 1,024 fluid particles is a low-resolution scientific prototype. The next major engineering steps are more faithful hydrostatic particle placement, a validated gravity tree/FMM for greater resolution, proper time-step convergence and the long isolated stability experiment. Promotion into stellar weapon gameplay requires stable equilibrium on the corresponding physical scales.
