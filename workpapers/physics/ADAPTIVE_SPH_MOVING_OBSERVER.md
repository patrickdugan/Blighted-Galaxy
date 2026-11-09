# BG-PHY-08: Moving Kerr observers, local worldtubes and adaptive-h SPH

**Date:** 10 November 2026  
**Status:** Experimental code with CI benchmarks, pending fidelity gates for white-dwarf stability.  
**Prior milestones:** [BG-PHY-06](KERR_CAMERA_WHITE_DWARF.md), [BG-PHY-07](KERR_CAUSAL_RELAX.md).  
**Repo issue:** [#8, converge cold white-dwarf SPH](https://github.com/patrickdugan/Blighted-Galaxy/issues/8).

This work extends the clean-room Blighted Galaxy reconstruction. The surviving stories require relativistic delayed perception and stellar shearing; the original lost Rust engine has not been recovered. The physically derived operations below remain distinct from the fictional weapons, oversized shields and gameplay collision envelopes.

## 1. Moving observers in stationary Kerr

**References:** [Bardeen–Press–Teukolsky ZAMO construction](https://www.astro.sunysb.edu/lattimer/AST390/gr_kerr.pdf); [null geodesic constants in the locally nonrotating frame](https://doi.org/10.1088/0067-0049/207/1/6); [circular Kerr emitter optics](https://academic.oup.com/mnras/article/475/3/3829/4795326).

BG-PHY-07 reconstructed incoming equatorial Kerr photons from a stationary observer at radius 85 gravitational radii. The current solver gives the observer an equatorial **circular timelike Kerr geodesic** and rotates the ray launch event as the observer moves in Boyer–Lindquist azimuth. A second observer mode follows a local zero-angular-momentum observer (ZAMO) at the same radius for comparison.

Let r_g=GM/c² and a=Jc/(GM²). In the equatorial plane, dimensionless Boyer–Lindquist radius x, Kerr Delta=x²−2x+a², and A=(x²+a²)²−a² Delta:

    lapse       = sqrt(x² Delta/A)
    g_phiphi    = A/x²
    omega_ZAMO  = (2ax/A) c/r_g
    Omega_orbit = s(c/r_g)/(x^(3/2)+s a)

The orbiting observer's local azimuthal velocity relative to the ZAMO tetrad is

    beta = (Omega_orbit-omega_ZAMO) sqrt(g_phiphi) r_g/(c lapse)
    gamma = 1/sqrt(1-beta²)
    d(tau)/dt = lapse/gamma.

For a **past-directed** photon, the local azimuth direction n_phi observed in the moving camera transforms into the ZAMO direction through

    n_phi,ZAMO = (n_phi,camera - beta)/(1 - beta n_phi,camera).

The solver maps this direction back into the camera's equatorial pinhole alpha coordinate. It then uses the BG-PHY-06 ZAMO null-ray integrator in the stationary Kerr background, rotating its entire photon path by the actual observer azimuth at reception. Kerr axial symmetry makes this rotation legitimate.

This is a local Lorentz aberration transform for a circular observer, and it does not extend to an accelerating arbitrary observer's full transported tetrad. The source, however, can follow any **recorded equatorial timelike track**, provided its sampled worldline passes the chronology and subluminal continuity checks.

## 2. Target collision volume in the local emitter metric

Earlier BG-PHY-07 ray matching employed a circular target envelope measured in a Euclidean map of Boyer–Lindquist x/y coordinates. This milestone replaces the chord intersection with a **local emitter-rest-frame spatial metric estimate**.

In the equatorial plane:

    g_tt=-[1-2/x]
    g_tphi=-2a/x
    g_rr=x²/Delta
    g_phiphi=A/x².

For a recorded source event, decompose the chart velocity into the radial rate dx/dt_geometric and azimuthal rate dphi/dt_geometric. Normalize its timelike four-velocity:

    U^t = 1/sqrt(-[g_tt + 2g_tphi w + g_phiphi w² + g_rr v_r²])
    U^r=v_r U^t
    U^phi=w U^t

The source covector has U_r=g_rr U^r and U_phi=(g_tphi+g_phiphi w)U^t. The rest-space projection tensor is h_mu,nu=g_mu,nu+U_mu U_nu. At equal chart time, a sufficiently small displacement (dr,dphi) receives a local squared distance

    ds_rest² ≈ [g_rr+U_r²]dr² + 2U_r U_phi dr dphi +
               [g_phiphi+U_phi²]dphi².

The photon/source matcher minimizes this local distance over each short, past-directed photon chord, using the recorded source position and velocity at the corresponding coordinate time. The collision envelope is thus more closely connected to the source's physical rest space. Static weak-field tests recover ordinary Euclidean distances; near Kerr strong fields, circular-orbit source targets acquire a metric-dependent distortion.

**Scope limitation:** The interval is a tangent-plane approximation to a finite source body. Integrating the physical worldtube boundary through four-dimensional Kerr spacetime is future work. The source envelope in the gameplay benchmarks spans roughly one gravitational radius at orbital radius 20r_g, so finite-size corrections and ray shooting resolution can still affect the recovered impact time.

## 3. Numerical moving-observer benchmark

Central mass: 10 solar masses. Observer: equatorial circular test-body at 85r_g, prograde orbit. Source: circular Kerr test-body at 20r_g, recorded in coordinate time. Target rest envelope: 1r_g. Rays are selected from a finite local angular bundle. The test independently checks the observer proper-time factor and spin-dependent aberration. The photon history records the earlier emission time and the source location at reception.

| Signed spin | Local observer beta | Proper / coordinate time | Curved photon delay | Source apparent displacement |
|---:|---:|---:|---:|---:|
| 0.00 | +0.109764 | 0.982194 | 3.382140 ms | 14.979 r_g |
| +0.75 | +0.109454 | 0.982228 | 3.381212 ms | 14.857 r_g |
| −0.75 | +0.110084 | 0.982160 | 3.381753 ms | 15.098 r_g |

These results are scenario-specific. A separate test injects a smooth accelerated perturbation into a historical source track and confirms that the solver finds an earlier event without requiring a circular-orbit source. Incomplete source histories yield no acquisition.

The visual program `cargo run --release --bin kerr_moving_lab` compares ZAMO and moving observers. Keys 1–3 select spin, W switches frame, F toggles automatic optical acquisition, A/D steers a manual local photon probe, Q/E changes speed, Z/X zooms and R resets. Cyan draws the curved incoming ray, yellow draws the received emission, and red marks the source's present position for debugging. This is a *perception model*, not a general-purpose ballistic shooting solution.

## 4. Variable smoothing-length cold white-dwarf SPH

**References:** [Price (2012), Smoothed particle hydrodynamics and magnetohydrodynamics](https://www.sciencedirect.com/science/article/pii/S0021999110006753), especially the variational grad-h correction; [Price & Federrath, PHANTOM smoothing-length formulation](https://users.monash.edu.au/~dprice/pubs/turbulence/pricefederrath2010.pdf).

The existing cold electron-degenerate EOS and 1D Newtonian hydrostatic reference are unchanged. The new `sph_adaptive.rs` module introduces a particle-dependent smoothing length, iterated jointly with its kernel density:

    rho_i = sum_j m_j W(|r_i-r_j|,h_i)
    h_i = eta [m_i/rho_i]^(1/3), eta=1.2.

The 3D cubic kernel has support 2h_i. The density and h are solved by underrelaxed fixed-point iteration, clamped to documented bounds. The self-consistent pressure-gradient correction is

    Omega_i = 1 + [h_i/(3rho_i)] sum_j m_j dW_ij(h_i)/dh_i

and the symmetric pair-force formula is

    a_i,pressure = - sum_j m_j [
        P_i/(Omega_i rho_i²) grad_i W_ij(h_i) +
        P_j/(Omega_j rho_j²) grad_i W_ij(h_j)].

These pressure forces conserve **total linear momentum** by pair symmetry. Gravity is Plummer-softened with a single fixed softening parameter across all pairs; it is **not yet** a fully adaptive variational gravitational-softening scheme. An Omega floor guards extreme numerical singularities, while benchmarks report the minimum correction and any smoothing-length cap saturation. This is an experimental solver suitable for controlled comparison.

For relaxation the internal velocity is damped relative to the stellar bulk frame. A subsequent undamped hold measures physical free evolution. The code retains the fixed-h method as the baseline, including its original tests and rendering.

## 5. Results against the hydrostatic reference

Cold electron EOS: mu_e=2, M=0.65 solar masses, hydrostatic R=8,402.21 km. Each row uses the **same damping interval** for the two methods at a particular N. Data were measured in Linux GitHub Actions using the BG-PHY-08 headless benchmark.

| Particles | Relax duration | Fixed-h radius | Adaptive-h radius | Adaptive radius bias |
|---:|---:|---:|---:|---:|
| 64 | 5.0 s | 4,380.84 km | 6,023.65 km | −28.3% |
| 128 | 5.0 s | 5,158.92 km | 6,297.11 km | −25.1% |
| 256 | 2.5 s | 7,047.02 km | 7,194.72 km | −14.4% |
| 512 | 2.5 s | 7,882.52 km | 7,706.85 km | −8.3% |

Adaptive smoothing improves the radius at N=64, 128 and 256 for these fixed-duration comparisons; at N=512 the fixed-h star has a smaller radius error (−6.2% versus −8.3%). At N=128 the adaptive run reaches minimum Omega≈0.815, smoothing lengths of 2.06–3.33 Mm and no hard-cap saturation. Its final RMS force is 5.76% of the initial force, compared with 7.84% for fixed h. At N=512 adaptive h reduces the force residual to 5.61% of the initial value, compared with 36.14% for fixed h. The density-profile error and stellar radius remain fidelity blockers.

The density log-RMS remains nonmonotonic. At N=128 the fixed-h result is approximately 0.2525 while adaptive-h is 0.3368 for the same five-second interval. At N=512 the corresponding values are 0.3327 and 0.4826 over 2.5 seconds. A smaller residual force at a displaced or smoothed density configuration is insufficient to certify hydrostatic equilibrium.

**Free-evolution gate:** The adaptive relaxed stars at N=64 and N=128 were evolved without damping for approximately 7.87 s, or three hydrostatic dynamical times. The outer radius expanded by **12.34% at N=64** and **8.11% at N=128**. Relative energy drift was respectively 3.19×10^−7 and 6.92×10^−8, and total momentum remained conserved to numerical precision. The stars are dynamically out of equilibrium despite the small numerical energy drift. The 5% stability gate has not been reached.

### Acceptance status for Issue #8

The previously requested **radius bias ≤5%** has not been reached. The free-evolution stability test also fails a reasonable ±5% radius-change gate. Keep the issue open. Follow-up work should include higher-resolution particle sampling and convergence, adaptive gravity softening with matching grad-epsilon corrections, refined hydrostatic particle initial conditions and a damped-to-undamped transition that stabilizes over multiple stellar dynamical times.

There is no justification to retune the physical electron EOS merely to force the radius to match. The numerical discrepancy belongs in the particle discretization and solver implementation.

## 6. Reproduction and artifact outputs

    cargo test --all-targets
    cargo run --release --bin phy08_bench -- --csv bg-phy08-adaptive.csv
    cargo run --release --bin kerr_moving_lab
    cargo run --release --bin white_dwarf_lab

In the updated white-dwarf lab, H initializes the old fixed-h relaxation and J initializes the new density-adaptive grad-h run, each as an explicitly selected method. Subsequent live stellar evolution uses the same force family that produced its initial particle state. The numerical HUD reports the chosen scheme, residual-force ratio and density mismatch. The original GDD and fiction remain in their dedicated folders, untouched.
