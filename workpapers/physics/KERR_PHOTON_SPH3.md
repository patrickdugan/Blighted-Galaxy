# BG-PHY-05: Kerr photons, 3D stellar fluid, GPU particle batches

**Status:** New clean-room physics code, 9 November 2026.
**Prior foundation:** BG-PHY-04 in workpapers/physics/STRONG_FIELD_CAUSAL_SPH.md.
**Sources:** The surviving Blighted Galaxy stories and game-design treatment describe relativistic combat, mass shearing, posthuman duels and stellar explosions. The equations below provide external physical reference models and explicitly approximate the fictional engineering.

## 1. Equatorial Kerr frame dragging and circular orbits

Reference: [SageManifolds Kerr geodesics](https://sagemanifolds.obspm.fr/kerrgeodesic_gw/reference/kerr_spacetime.html) and [Photon Emissions from Kerr Equatorial Geodesic Orbits](https://link.springer.com/article/10.1140/epjc/s10052-024-13113-w).

In geometric units, set r_g=GM/c² and x=r/r_g. The signed dimensionless spin is a=Jc/(GM²), restricted in this prototype to |a|<=0.998. The Kerr functions are:

    Delta = x² - 2x + a²
    r_outer / r_g = 1 + sqrt(1-a²)
    r_equatorial_ergosphere / r_g = 2
    A = (x²+a²)² - a² Delta
    ZAMO_omega = (2 a x / A)(c/r_g)

For a circular equatorial timelike orbit, with direction s=+1 (positive angular momentum) or s=-1 (negative):

    Omega = s (c/r_g) / (x^(3/2)+s a)
    E = [x^(3/2)-2 sqrt(x)+s a] /
        [x^(3/4) sqrt(x^(3/2)-3 sqrt(x)+2 s a)]
    L/(c r_g) = s[x²-2s a sqrt(x)+a²] /
        [x^(3/4) sqrt(x^(3/2)-3 sqrt(x)+2 s a)]

The ISCO radius follows Bardeen–Press–Teukolsky:
Z1=1+(1-a²)^(1/3)[(1+a)^(1/3)+(1-a)^(1/3)],
Z2=sqrt(3a²+Z1²), and x_ISCO=3+Z2-s sqrt((3-Z1)(3+Z1+2Z2)).

For zero spin the two branches return the Schwarzschild ISCO x=6.

## 2. Equatorial null trajectories

The code uses E=1, equatorial Carter constant Q=0, signed impact parameter b=L/E and P=x²+a²-a b. The radial potential is:

    R = P² - Delta (b-a)²
    (dx/dlambda)² = R/x⁴
      = 1 + (a²-b²)/x² + 2(b-a)²/x³

A second-order radial ODE continues through turning points:

    d²x/dlambda² = (b²-a²)/x³ - 3(b-a)²/x⁴
    dphi/dlambda = [b-a + a P/Delta]/x²
    d(ct/r_g)/dlambda = [a(b-a)+(x²+a²)P/Delta]/x²

RK4 integrates these equations outside the horizon, while the code monitors the radial constraint residual. The equatorial photon circular radii are 2[1+cos((2/3)arccos(-s a))]. Setting R=0 there gives the critical positive prograde and negative retrograde impact parameters. For a=0, photon orbit radius is 3 r_g and the two critical b values are +/-3sqrt(3) r_g.

Each ray stops at the Boyer–Lindquist coordinate cutoff r_outer+0.02r_g or on return to its original outer radius. These are null geodesics of a stationary single-hole background in the equatorial plane. Kerr photon trajectories are not yet wired into the default battle's causal ledger; that still uses flat-spacetime light propagation. A true image of the sky around a spinning black hole requires 3D null rays, Carter's constant, radiative transfer and a camera model.

## 3. Three-dimensional SPH-like fluid

References: [PySPH kernel documentation](https://pysph.readthedocs.io/en/main/reference/kernels.html); [Abaqus SPH analysis](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEANLRefMap/simaanl-c-sphanalysis.htm).

The module sph3.rs uses a three-dimensional cubic-spline kernel with q=r/h:

    W3 = (1/(pi h³)) * F(q)
    F(q) = 1-1.5q²+0.75q³       for 0<=q<1
         = 0.25(2-q)³          for 1<=q<2
         = 0                     for q>=2

Density rho_i=sum_j m_j W_ij is measured in kg/m³. A fixed-entropy polytropic closure describes pressure:

    P = K rho^gamma
    gamma = 5/3
    u_internal = K rho^(gamma-1)/(gamma-1)

The prototype uses K=3e8 in SI-compatible units, h=3e6 m and Plummer gravitational softening of 1.2e6 m. These values are research placeholders, and the polytropic model is not an equilibrated white-dwarf equation of state.

Pairwise pressure and gravity exchange equal and opposite momentum, including gravitational recoil on a central massive body. A velocity-Verlet integrator evolves all 3D coordinates. Tests monitor total mass, linear momentum, kinetic + gravitational + polytropic internal energy and short-step convergence. Mass and momentum conservation are algebraic symmetries of the chosen numerical scheme; energy shows finite timestep drift.

The existing 2D tide-to-self-gravity threshold is retained as a transition trigger. At crossing, the new lab initializes 128 symmetric 3D parcels inside a white-dwarf radius, preserving initial mass and bulk momentum. It then runs the 3D fluid. This is a toy stellar disruption; true hydrodynamic equilibrium, electron degeneracy, shock heating, nuclear burning, gravitational radiation and relativistic capture are outstanding.

## 4. GPU-backed triangle batching

Reference: [macroquad Mesh and Vertex API](https://docs.rs/macroquad/latest/macroquad/models/struct.Mesh.html).

Particle billboards use 4 vertices, 2 triangles and 6 indices each. The renderer groups up to 4,096 sprites into one Mesh call, staying inside the u16 vertex-index limit. A deterministic 9,500-sprite background uses three mesh draws; fluid sprites use an additional dynamic mesh.

The GPU draws the indexed triangles. SPH gravity, animation state, 3D projection and dynamic mesh creation still happen on the CPU. This stage is GPU-backed rendering, not GPU compute/SPH or hardware ray tracing. Reusable buffers, instancing and compute are subsequent performance experiments after measuring the CPU and GPU costs separately.

## 5. Reproduction

Run the following in the repository root:

    cargo test --all-targets
    cargo check --all-targets
    cargo run --release --bin kerr_sph3_bench
    cargo run --release --bin kerr_sph3_lab

The headless benchmark prints Kerr frame dragging, critical impact parameters, photon radial constraint residuals, 3D fluid energy/momentum drift and particle mesh batching. The lab has four modes: Schwarzschild photons (1), Kerr +0.85 (2), Kerr -0.85 (3), and 3D tidal debris (4). Use Q/E to adjust rate, A/D to rotate the debris projection, Z/X or wheel to zoom, R to reset.

A successful build or unit test does not establish rendered-playability. Inspect the Windows binary on the target RTX 5080 and profile draw calls, CPU frame time and SPH particle count separately.

## 6. Follow-on physics

- General Kerr null rays with nonzero Carter constant and a camera/sky imaging model.
- Retarded optical histories in curved spacetime, including Shapiro delay and gravitational bending.
- A hydrostatically equilibrated 3D degenerate-matter white dwarf, convergence with respect to h/N/dt, and energy-momentum accounting around capture.
- Persistently mapped particle buffers, actual GPU instancing or compute SPH where benchmarks justify it.
