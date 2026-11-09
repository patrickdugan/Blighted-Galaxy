# BG-PHY-06: Kerr camera, cold electron EOS and hydrostatic white dwarfs

**Status:** Testable Rust milestone, 9 November 2026.  
**Modules:** src/kerr_camera.rs, src/white_dwarf.rs, modified src/sph3.rs  
**Runnables:** kerr_camera_lab, kerr_camera_bench, white_dwarf_lab, white_dwarf_bench

The original Blighted Galaxy fiction describes spinning black-hole weapons, relativistic flight, transhuman fights and tidal destruction. These modules explore the corresponding physical models. The simulations remain a clean-room reconstruction, with deliberately visible assumptions and test boundaries.

## 1. Non-equatorial null propagation around Kerr

**Sources:** the separability of Kerr null geodesics and Carter's constant are described in [Kerr geodesic motion in Mino time](https://academic.oup.com/ptep/article/2013/6/063E01/1625420), and the equatorial reference case is in [BG-PHY-05](KERR_PHOTON_SPH3.md).

The camera advances a ray backward from a local observer at Boyer–Lindquist radius r_o=85 GM/c². The observer is a zero-angular-momentum observer (ZAMO), and its local tetrad supplies a past-directed photon with unit local energy. The local pinhole image plane uses alpha and beta, and the inclination can be varied.

Use dimensionless r_g=GM/c² and a=Jc/(GM²). Define:

    Sigma = r² + a² cos²(theta)
    Delta = r² - 2r + a²
    A     = (r²+a²)² - a² Delta sin²(theta)
    lapse = sqrt(Sigma Delta/A)
    omega = 2 a r/A
    sqrt(g_phiphi) = sqrt(A sin²(theta)/Sigma)

The ZAMO local ray has p^(time)=-1 and spatial direction (n_r,n_theta,n_phi). Its signed Killing energy and angular momentum are

    E_past = -lapse + omega * L
    L      = n_phi * sqrt(g_phiphi)

With b=L/E_past, normalized non-equatorial separation constant is

    q = (Sigma * n_theta² / E_past²)
        + cos²(theta) * (b²/sin²(theta)-a²).

For a null ray with E_past<0, use the Mino-time radial and polar separation equations (normalized to |E_past|):

    R(r) = [(r²+a²)-a b]² - Delta * [(b-a)²+q]
    Theta(theta) = q+a² cos²(theta)-b² cot²(theta).

For numerical stability at the observer's finite large radius, the integrator uses u=1/r and a second-order radial ODE:

    F(u) = 1+(a²-b²-q)u²+2[(b-a)²+q]u³-a² q u⁴
    (du/dlambda)² = F(u)
    d²u/dlambda² = (a²-b²-q)u+3[(b-a)²+q]u²-2a²q u³.

The polar equation is

    (dtheta/dlambda)² = Theta(theta)
    d²theta/dlambda² =
        -a² sin(theta)cos(theta)+b² cos(theta)/sin³(theta).

The azimuth and coordinate-time rates carry a negative sign for past-directed tracing:

    P = r²+a²-a b
    dphi/dlambda = -[a P/Delta+b/sin²(theta)-a]
    dt/dlambda = -[(r²+a²)P/Delta+a(b-a sin²(theta))].

A fourth-order Runge–Kutta integrator advances six quantities (u,u',theta,theta',phi,t). Radial and polar constraint errors are measured, with a reduced step approaching the Boyer–Lindquist horizon cutoff. Polar coordinate reflection is used at the axis. This is a numerical chart treatment, pending a regular tetrad implementation for trajectories through a pole.

The rays terminate at the horizon-coordinate cutoff, at the observer radius after scattering, or on crossing an opaque equatorial thin disc between the spin-appropriate ISCO and 28 r_g. The disc illustrates gravitational redshift and circular geodesic Doppler factors through

    g = 1/[u^t (1-Omega*b)].

Its color gradient is illustrative; the model has no radiative-transfer equation, optically thin plasma emissivity, physically calibrated luminosity or polarization. The distant background is a synthetic angular star texture that becomes warped by geodesic mapping. Photon geodesics are non-equatorial and support q!=0; the live combat scene still uses a flat-spacetime causal worldline solver until its null-ray propagation is explicitly integrated.

### Camera controls and acceptance

Run:

    cargo run --release --bin kerr_camera_lab

Keys 1–4 select a=0,+0.45,+0.85,-0.85. Up/Down changes observer inclination, Z/X adjusts field of view, R recomputes. The initial frame is rendered at 240×144 pixels and enlarged on screen. CPU ray integration builds the texture; the GPU displays it.

The testable image renderer reports captured/disc/sky/unresolved pixel counts, RK4 step counts and radial/polar constraint residuals. It verifies a Schwarzschild central silhouette, off-equatorial rays, Carter constants, equatorial reduction, disc redshift, and spin-distorted output.

## 2. Fully relativistic electron degeneracy pressure

**Sources:** [cold white dwarf EOS overview](https://academic.oup.com/mnras/article/464/4/4349/2562621), [white-dwarf structure equations](https://dca.ue.ucsc.edu/system/files/dca/1188/1188.pdf).

The hydrostatic solver models a zero-temperature, fully degenerate and ideal electron gas, supported by ions contributing baryonic mass density. For carbon/oxygen plasma, mean molecular weight per electron mu_e≈2. The electron rest mass m_e, atomic mass unit m_u, Planck constant hbar, c and G are SI constants.

    n_e = rho/(mu_e*m_u)
    x = hbar*(3 pi² n_e)^(1/3)/(m_e c)
    P(rho) = [m_e^4 c^5/(24 pi² hbar³)]
             * [x(2x²-3)sqrt(1+x²)+3asinh(x)].

A small-x series evaluates the same pressure with controlled precision. The pressure smoothly transitions from nonrelativistic P∝rho^(5/3) toward ultrarelativistic P∝rho^(4/3).

The specific enthalpy H, its inverse and the derivative of pressure are

    H = [m_e c²/(mu_e m_u)] * (sqrt(1+x²)-1)
    rho(H) = rho0 [y(2+y)]^(3/2),
        y=H/[m_e c²/(mu_e m_u)],
        rho0=mu_e*m_u*(m_e*c/hbar)^3/(3pi²)
    dP/drho = [m_e c²/(3mu_e m_u)] * x²/sqrt(1+x²).

The electron kinetic internal energy per baryonic mass is

    u = H - P/rho,

which supplies a thermodynamically consistent cold barotropic energy function to the SPH energy diagnostic.

## 3. Newtonian hydrostatic equilibrium

**References:** [Lane–Emden stellar structure](https://aprsa.villanova.edu/files/pols.pdf); [white-dwarf hydrostatic equations](https://dca.ue.ucsc.edu/system/files/dca/1188/1188.pdf).

The integrator solves the physical-radius initial value problem:

    dH/dr = -G m(r)/r²
    dm/dr = 4pi r² rho(H).

Near the regular center, the code initializes at a small positive radius using the series m≈4pi rho_c r³/3 and H≈H_c−2pi G rho_c r²/3. RK4 steps of 10 km continue to the H=0 surface. Linear interpolation across the final step estimates the stellar surface and total mass.

The solver stores radial samples of density and enclosed mass for a given central density. A separate bracketed central-density search obtains a chosen stellar mass and an inverse mass-CDF to place equal-mass 3D SPH parcels at the correct stratified radii.

For mu_e=2, the zero-temperature Newtonian reference gives approximately:

| Central density (kg/m³) | Mass (solar masses) | Radius (km) |
|---:|---:|---:|
| 1e8 | 0.149 | 16,369 |
| 1e9 | 0.401 | 10,958 |
| 3e9 | 0.588 | 8,953 |
| 1e10 | 0.815 | 7,085 |
| 1e11 | 1.180 | 4,328 |

A target of **0.65 solar masses** gives rho_c≈4.19e9 kg/m³ and R≈8,402 km. Test assertions verify these results at a practical relative tolerance and check the mass-CDF inversion.

The numerical values assume a spherical, Newtonian, nonrotating, zero-temperature white dwarf with uniform mu_e=2. They omit ionic Coulomb physics, temperature, magnetic fields, inverse beta decay, general-relativistic structure corrections and relativistic instability.

## 4. Integrating the EOS into 3D fluid debris

The existing SPH3 module now exposes a configurable pressure law: the earlier gamma=5/3 polytropic gameplay proxy or the cold relativistic electron EOS. For either case, its pairwise pressure interactions continue to conserve linear momentum by construction, and the energy diagnostic includes the matching barotropic internal energy.

The method seed_from_hydrostatic populates antipodal particle pairs with radii drawn from the physical enclosed-mass profile. Equal masses, symmetric displacement and common bulk velocity preserve the initial total mass, center of mass and linear momentum to rounding tolerance. They form a stratified approximation to a continuum equilibrium star.

Run:

    cargo run --release --bin white_dwarf_lab

Keys 1–3 select 0.35, 0.65 and 0.95 solar masses. The lab shows a mass–radius curve alongside a rotating projection of a 192-particle degenerate SPH white dwarf near a massive black hole. Space evolves/pauses, Q/E changes rate, A/D rotates projection, and Z/X adjusts zoom.

The 3D fluid model remains coarse. Its initial radial density is hydrostatic, but particle pressure gradients and smoothing lengths do not guarantee discrete SPH force equilibrium. The lab reports energy and momentum drift for investigation. A physical disruption model additionally requires degenerate EOS composition gradients, shock heating and a proper GR capture boundary.

## 5. Numerical and rendering reproduction

The new headless programs use no graphics display:

    cargo test --all-targets
    cargo run --release --bin kerr_camera_bench
    cargo run --release --bin white_dwarf_bench

Optional generated outputs:

    cargo run --release --bin kerr_camera_bench -- --output kerr-camera.ppm
    cargo run --release --bin white_dwarf_bench -- --output wd-mass-radius.csv

CI packages a small PPM Kerr camera image and a white-dwarf mass–radius CSV as build artifacts, as well as the Windows game binaries. All generated camera colors are explicitly synthetic; the nonlinear ray paths supply the physically derived part of the image.

## 6. Follow-on implementation questions

- General-relativistic causal optical histories for moving enemies, including source/observer worldlines and Shapiro delays.
- Adaptive step control and regular Cartesian/tetrad coordinates near polar axes and horizons.
- Kerr disc radiative transfer with polarization and emitter velocity distributions.
- Resolution-dependent hydrodynamic convergence studies with appropriate degenerate thermodynamics and stable initial SPH equilibrium.
- Resource-scaled GPU photon shader or compute ray tracing, with benchmarked runtime on the RTX 5080.
