# BG-PHY-26 | Three-moment free-surface depth and nonuniform SPH volume bias

**Date:** 21 November 2026  
**Parent:** [BG-PHY-25](FREE_BOUNDARY_ENTHALPY_GHOST.md)  
**Tracking:** [Issue #31](https://github.com/patrickdugan/Blighted-Galaxy/issues/31)  
**Physical acceptance:** [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8)

## 1. Scientific question and motivation

[BG-PHY-25](FREE_BOUNDARY_ENTHALPY_GHOST.md) showed that while the kernel first moment $Q_i$ accurately recovers the outward radial normal on the spherical $0.65\,M_\odot$ white dwarf (mean outermost-shell cosine $0.9870$ for Cubic M4 and $0.9973$ for Wendland C4 at $N=256$), both single-ratio planar depth estimators severely underestimate the distance to the independent 1D hydrostatic surface ($0.5402\,h$ for M4, $0.3850\,h$ for C4):

- Half-space partition $S_i$ infers $0.0616\,h$ (M4) and $0.0217\,h$ (C4), with $62.5\%$ of C4 outermost-shell particles having $S_i < 0.5$.
- The uniform-bias-invariant ratio $|Q_i| / S_i$ infers only $0.1202\,h$ (M4) and $0.0351\,h$ (C4), because spatially varying SPH volume weights $V_j = m_j / \rho_{\mathrm{SPH},j}$ break the cancellation of common multiplicative bias.

BG-PHY-26 tests whether jointly inverting the **first three normal kernel moments** $(M_0, M_1, M_2)$ under a local linear normal volume-bias model $V_{\mathrm{est}}(z) = c (1 + \beta z) V_{\mathrm{true}}(z)$ recovers accurate free-surface depths and improves outermost-shell variational pressure support and radial homology equilibrium.

All equations use only particle positions, masses, exact-root SPH smoothing lengths and densities, the chosen isotropic kernel, and the unchanged cold Fermi electron EOS. The independent 1D hydrostatic white-dwarf profile enters **exclusively in post-hoc scoring**.

## 2. Exact 1D quadratures and three-moment inverse

At particle $i$ with safeguarded Newton smoothing length $h_i$, let $q_{ij} = (x_j - x_i)/h_i$, $w_{ij} = (m_j / \rho_{\mathrm{SPH},j}) W(\|x_j - x_i\|, h_i)$, and

```text
Q_i = sum_j w_ij q_ij,
n_out,i = -Q_i / |Q_i|,
z_ij = n_out,i · q_ij.
```

The three measured normal moments on the same neighbor cloud are

```text
M0 = sum_j w_ij,
M1 = sum_j w_ij z_ij,
M2 = sum_j w_ij z_ij².
```

For an isotropic unit-normalized kernel truncated by a planar free surface at $z = a$ ($0 \le a \le 2$, occupied fluid half-space $z \le a$), the ideal moments are

```text
I_k(a) = integral_{|q|<2, q·n_out <= a} (q·n_out)^k W(|q|, 1) d³q,    k = 0, 1, 2, 3.
```

Integrating analytically over azimuth $\phi \in [0, 2\pi)$ and polar cosine $\mu \in [-1, \min(1, a/q)]$ reduces each $I_k(a)$ to a 1D piecewise-polynomial radial integral:

- For $0 \le q \le a$: $I_k$ receives $4\pi q^{k+2} W(q,1)/(k+1)$ for even $k$ and $0$ for odd $k$.
- For $a < q < 2$: $I_k$ receives $2\pi W(q,1) [a^{k+1} q - (-1)^{k+1} q^{k+2}] / (k+1)$.

Because $W(q,1)$ is piecewise cubic on $[0,1] \cup [1,2]$ for M4 and degree-8 polynomial on $[0,2]$ for Wendland C4, the integrand on every smooth subinterval partitioned at $\{0, \min(a,1), \max(a,1), 2\}$ has polynomial degree at most $13$. An **8-point Gauss-Legendre quadrature** per smooth segment is exact for polynomials up to degree $15$ and reproduces the exact rational endpoint limits to IEEE-754 precision ($< 10^{-15}$):

| Kernel | $I_0(0)$ | $I_1(0)$ | $I_2(0)$ | $I_3(0)$ | $I_0(2)$ | $I_1(2)$ | $I_2(2)$ | $I_3(2)$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Cubic M4 | $1/2$ | $-31/140$ | $3/20$ | $-127/1008$ | $1$ | $0$ | $3/10$ | $0$ |
| Wendland C4 | $1/2$ | $-35/192$ | $4/39$ | $-15/208$ | $1$ | $0$ | $8/39$ | $0$ |

By Leibniz's rule across the planar cap $z = a$ ($\rho_\perp d\rho_\perp = q\,dq$), the exact derivatives with respect to boundary depth are

```text
dI_k/da = 2*pi * a^k * integral_a^2 q W(q, 1) dq,    k = 0, 1, 2, 3.
```

Under $V_{\mathrm{est}}(z) = c (1 + \beta z) V_{\mathrm{true}}(z)$,

```text
M0 = c [I0(a) + beta I1(a)],
M1 = c [I1(a) + beta I2(a)],
M2 = c [I2(a) + beta I3(a)].
```

Eliminating $c$ via $r_1 = M_1 / M_0$ and $r_2 = M_2 / M_0$ gives

```text
beta(a) = [I1(a) - r1 I0(a)] / [r1 I1(a) - I2(a)],
F(a)    = [I2(a) + beta(a) I3(a)] / [I0(a) + beta(a) I1(a)] - r2 = 0.
```

Clearing the common algebraic factor $r_1 I_1(a) - I_2(a)$ yields the removable-pole-free rational form

```text
r2_hat(a) - r2 = [(I2² - I1 I3) + r1 (I0 I3 - I1 I2)] / (I0 I2 - I1²) - r2,
```

whose denominator $I_0(a) I_2(a) - I_1(a)^2 > 0$ is strictly positive on $[0,2]$ by Cauchy-Schwarz on the positive half-space measure. `src/sph_three_moment.rs` scans 256 brackets across $[0,2]$, polishes every sign change with safeguarded Newton-bisection, verifies the candidate directly on the original $\beta(a)$ and $F(a)$ equations, and requires all of:

1. $c > 0$ and $|r_1 I_1(a) - I_2(a)| > 10^{-12}$;
2. strictly positive reconstructed source volume $1 + \beta z > 10^{-4}$ across the occupied compact-support half-space $z \in [-2, a]$ and across all realized neighbor coordinates $z_{ij}$;
3. uniqueness on $[0,2]$ (multiple valid roots are marked `ambiguous_multiple_roots` and rejected);
4. Jacobian conditioning $|F'(a)| \ge 10^{-3}$ (interior particles near $a \to 2$ where $(2-a)^8 \to 0$ and $|F'(a)| \to 0$ are flagged `ill_conditioned_interior`).

Whenever a particle is not uniquely `identified`, **no vacuum $H=0$ ghost is permitted** and the estimator falls back identically to the baseline exact-root SPH density $\rho_{\mathrm{SPH},i}$.

## 3. Stage A results: synthetic verification, curvature bias, and white-dwarf depth scoring

### 3.1 Synthetic planar recovery and failure modes

Full Stage-A numerical outputs are archived in [`workpapers/physics/data/bg-phy26/bg-phy26-synthetic.csv`](data/bg-phy26/bg-phy26-synthetic.csv):

- **Exact planar linear-bias recovery:** Across $a \in [0, 2]$, $c \in [0.5, 1.5]$, and $\beta \in [-0.4, 0.4]$ for both M4 and C4, the solver recovers unique mathematical roots with $|a - a_{\mathrm{true}}| \le 6 \times 10^{-15}$ on $a \in [0, 1.4]$ and $|F'(a)| \in [0.0038, 0.349]$ (depth sensitivity $|\partial a / \partial r_2| \in [2.86, 264]$). Near full interior support ($a \ge 1.8$ for C4, $a = 2.0$ for M4), $|F'(a)|$ drops below $10^{-3}$ and is cleanly classified as `ill_conditioned_interior`. At $a = 0.5, \beta = -0.35$, $|Q|/S$ shifts by $> 0.04$ while three-moment inversion recovers $a = 0.5$ to $< 10^{-15}$.
- **Irregular 3D finite-particle half-spaces:** On jittered 3D half-space clouds with $c = 1.25, \beta = -0.30, a_{\mathrm{true}} = 0.55$, the solver identifies $a$ within $0.03\,h$, reproduces constant enthalpy to $< 10^{-11}$ (de-biased affine fit), and preserves positive accurate affine and curved ($H = 10 d + 2.5 d^2$) extrapolations.
- **Spherical surface curvature (`spherical_cap_moments`):** A sphere of radius $R$ bends inward below its tangent plane by $\Delta z \approx -q_\perp^2 / (2 R/h)$. At $R/h = 2.8$ (comparable to $N=256$ white-dwarf smoothing scales) with $c = 0.70, \beta = -0.45$, pure spherical curvature shifts the planar three-moment depth inward by **$-0.097\,h$ to $-0.101\,h$** under M4 and **$-0.068\,h$ to $-0.071\,h$** under C4, decaying monotonically as $O(h/R)$ to $\le 0.015\,h$ at $R/h = 20$.
- **Quadratic normal volume curvature (`quadratic_volume_bias_moments`):** Adding negative normal curvature $\gamma = -0.30$ in $V_{\mathrm{est}}/V_{\mathrm{true}} = c(1 + \beta z + \gamma z^2)$ shifts the inferred depth inward by a further **$-0.050\,h$** (M4) and **$-0.039\,h$** (C4) while absorbing curvature into an apparent steeper $\beta$.

### 3.2 Score-only white-dwarf surface depth audit

Measured outermost-shell ($1/8$ mass) depth comparisons against the independent 1D cold-electron hydrostatic radius $R_{\mathrm{ref}} = 8402.207\;\mathrm{km}$ (archived in [`bg-phy26-n64-seed0-summary.csv`](data/bg-phy26/bg-phy26-n64-seed0-summary.csv), [`bg-phy26-n128-seed17-summary.csv`](data/bg-phy26/bg-phy26-n128-seed17-summary.csv), and [`bg-phy26-n256-initial-summary.csv`](data/bg-phy26/bg-phy26-n256-initial-summary.csv)):

| Configuration | Kernel | True ref depth | Partition depth | `\|Q\|/S` depth | **3-moment depth** | `\|Q\|/S` abs err | **3-moment abs err** | Identified | Mean $c$ | Mean $\beta$ | Geometry gate ($\le 0.15h$, $\ge 90\%$) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| N64 seed0 | M4 | 0.3555 h | 0.0029 h | 0.0278 h | **0.2718 h** | 0.3277 h | **0.1849 h** | 50.0% | 0.522 | −1.318 | FAIL |
| N64 seed0 | C4 | 0.2484 h | 0.0000 h | 0.0000 h | **0.1557 h** | 0.2484 h | **0.1401 h** | 75.0% | 0.454 | −2.254 | FAIL (id < 90%) |
| N128 seed17 | M4 | 0.4392 h | 0.0154 h | 0.0655 h | **0.2331 h** | 0.3736 h | **0.2518 h** | 87.5% | 0.602 | −0.851 | FAIL |
| N128 seed17 | C4 | 0.3114 h | 0.0021 h | 0.0001 h | **0.1447 h** | 0.3113 h | **0.2071 h** | 87.5% | 0.526 | −1.605 | FAIL |
| N256 seed0 | M4 | **0.5402 h** | 0.0616 h | 0.1202 h | **0.2217 h** | 0.4201 h | **0.3572 h** | 93.8% | 0.729 | −0.488 | FAIL |
| N256 seed0 | C4 | **0.3850 h** | 0.0217 h | 0.0351 h | **0.1472 h** | 0.3499 h | **0.2643 h** | 93.8% | 0.614 | −1.077 | FAIL |
| N256 seed17 | M4 | **0.5406 h** | 0.0707 h | 0.1259 h | **0.2192 h** | 0.4147 h | **0.3703 h** | 93.8% | 0.763 | −0.437 | FAIL |
| N256 seed17 | C4 | **0.3852 h** | 0.0253 h | 0.0384 h | **0.1516 h** | 0.3467 h | **0.2754 h** | 93.8% | 0.636 | −1.032 | FAIL |

**Physical interpretation of the Stage-A audit:**
1. The three-moment solver confirms BG-PHY-25's hypothesis that strong negative normal volume gradients ($\beta \approx -0.46$ on M4, $\beta \approx -1.05$ on C4, alongside $c \approx 0.61\text{–}0.76$) exist in the outer white-dwarf shell because wide kernel support smooths over denser interior shells so $\rho_{\mathrm{SPH}} > \rho_{\mathrm{true}}$ near $r_{\max}$, causing $V_{\mathrm{est}}/V_{\mathrm{true}} = \rho_{\mathrm{true}}/\rho_{\mathrm{SPH}}$ to decrease sharply outward.
2. Accounting for $(c, \beta)$ increases inferred surface depth by **$1.8\times$ on M4** and **$4.2\times$ on C4** relative to $|Q|/S$, with $93.8\%$ ($30/32$) identifiability and zero ambiguous multi-root cases at $N=256$.
3. Nevertheless, the inferred depth ($0.22\,h$ M4, $0.15\,h$ C4) remains **$59\%\text{–}62\%$ below** the continuum reference depth ($0.54\,h$ M4, $0.385\,h$ C4) and **fails the Stage-A geometry-quality gate**. Beyond spherical curvature ($\sim 0.07\text{–}0.10\,h$) and quadratic volume curvature ($\sim 0.04\text{–}0.05\,h$), the discrete shell seed places its outermost midpoint quantile at $r_{\max} = 7840.7\;\mathrm{km}$ ($6.68\%$ inside $R_{\mathrm{ref}}$): no particles exist in the outer half-quantile tail where degenerate density falls as $(R_{\mathrm{ref}} - r)^{3/2}$, so local particle moments reflect half a local shell thickness rather than the continuum polytropic tail out to $R_{\mathrm{ref}}$.

## 4. Stage B diagnostic screening and the variational ghost-spring obstruction

Although the Stage-A geometry gate already disqualifies any new vacuum-ghost promotion, evaluating the re-solved conservative scalar potential $E(x) = \sum_i m_i u_{\mathrm{cold}}(\hat{\rho}_i(x)) + E_{\mathrm{grav}}(x)$ exposes a second, structural variational mechanism that explains why moment-inferred vacuum ghosts weaken outermost-shell pressure support even when static surface densities improve.

Three research-only three-moment candidates are evaluated alongside `original`, `surface_mls`, and BG-PHY-25's `moment_ratio_ghost_blend`:

- `three_moment_debiased_affine_blend`: de-biased neighbor weights $w_{ij} / [c_i (1 + \beta_i z_{ij})]$ in affine enthalpy MLS **without** a vacuum ghost;
- `three_moment_ghost_blend`: scale-normalized weights $w_{ij}/c_i$ plus a virtual $H=0$ ghost at $a_i n_{\mathrm{out},i}$ with weight $\mathrm{clamp}(1 - I_0(a_i), 0, 0.5)$ when `identified` (strict fallback to $\rho_{\mathrm{SPH},i}$ otherwise);
- `three_moment_debiased_ghost_blend`: both de-biased neighbor weights and the three-moment $H=0$ ghost.

### 4.1 Initial $N=256$ angular-glass comparison (seeds 0 and 17)

| Candidate | M4 s0 last $P/G_{\mathrm{ref}}$ | M4 s17 last $P/G_{\mathrm{ref}}$ | C4 s0 last $P/G_{\mathrm{ref}}$ | C4 s17 last $P/G_{\mathrm{ref}}$ | M4 s0 penult $P/G_{\mathrm{ref}}$ | C4 s0 penult $P/G_{\mathrm{ref}}$ | M4 s0 $R^*$ bias | C4 s0 $R^*$ bias | C4 s0 last $\rho$ log-RMS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `original` | **0.7312** | **0.7067** | **0.6807** | **0.6684** | 0.8806 | 0.8451 | **−2.27%** | **−2.85%** | 0.7119 |
| `surface_mls` | 0.7286 | 0.6871 | 0.6609 | 0.6349 | 0.9660 | 0.9800 | −3.85% | −5.44% | **0.2895** |
| `moment_ratio_ghost_blend` | 0.5188 | 0.4946 | 0.3895 | 0.3705 | 0.8050 | 0.8269 | −8.00% | −10.95% | 1.3209 |
| `three_moment_debiased_affine_blend` | 0.6509 | 0.6078 | 0.5259 | 0.4885 | 1.0842 | 0.8840 | −3.35% | −4.20% | 0.6557 |
| `three_moment_ghost_blend` | 0.1177 | 0.0614 | **−0.2105** | **−0.3271** | **1.0610** | **0.9848** | −7.13% | −8.02% | 0.6418 |
| `three_moment_debiased_ghost_blend` | 0.1108 | 0.0558 | **−0.2286** | **−0.3334** | 0.9806 | 0.9896 | −7.14% | −8.02% | 0.6476 |

Every Cholesky solve succeeded across all configurations (`fit_fallback_fraction = 0`).

### 4.2 Analytical explanation: inter-shell variational ghost-spring coupling

At fixed positions, `three_moment_ghost_blend` cuts C4 last-shell density log-RMS error in half relative to `moment_ratio_ghost_blend` ($1.3209 \to 0.6418$, slightly better than `original` $0.7119$) and moves the global homology minimum $R^*$ outward by $+2.9\%$ on C4 ($-10.95\% \to -8.02\%$) and $+0.9\%$ on M4 ($-8.00\% \to -7.13\%$). Under uniform homology scaling $x \mapsto s x$, all inter-shell distances and $h_i$ scale proportionally so dimensionless moments $(r_1, r_2)$ and inferred depth $a_i$ are scale-invariant ($\partial a_i / \partial s = 0$ away from $h$ caps).

Under an **independent outward displacement of the outermost shell alone** ($\delta r_{\mathrm{last}} > 0$ with interior shells fixed), however:

1. Separating the outermost shell from the penultimate shell stretches $|z_{ij}|$ for interior neighbors and increases the dimensionless second moment ratio $r_2 = M_2 / M_0$ ($\partial r_2 / \partial r_{\mathrm{last}} > 0$).
2. Because $F'(a_i) > 0$, larger $r_2$ increases the inferred boundary depth ($\partial a_i / \partial r_{\mathrm{last}} = (\partial r_2 / \partial r_{\mathrm{last}}) / F'(a_i) > 0$).
3. Increasing $a_i$ moves the $H=0$ ghost outward and reduces both ghost weight $1 - I_0(a_i)$ and surface gate $\alpha(I_0(a_i))$, which **raises** the reconstructed enthalpy and density ($\partial \hat{\rho}_i / \partial a_i > 0$).
4. In the conservative cold-EOS internal energy $U = \sum_i m_i u_{\mathrm{cold}}(\hat{\rho}_i)$ with $\partial u_{\mathrm{cold}} / \partial \hat{\rho}_i = P_i / \hat{\rho}_i^2 > 0$, the positive chain-rule term

```text
d rho_hat_i / d r_last = (partial rho_hat_i / partial r_last)|_{a_i}
                       + (partial rho_hat_i / partial a_i) * (partial a_i / partial r_last)
```

opposes and overwhelms standard SPH expansion cooling $(\partial \hat{\rho}_i / \partial r_{\mathrm{last}})|_{a_i} < 0$, turning $F_{P,\mathrm{last}} = -dU / d r_{\mathrm{last}}$ strongly inward on the last shell while applying an equal-and-opposite outward reaction to the penultimate shell (`penult P/G` rises to $0.98\text{–}1.15$).

### 4.3 Post-3s-damping confirmation ($100\%$ identified smooth regime)

On the post-3s-damped $N=256$ seed-0 M4 `original` geometry (`tm_d = 0.1288h` vs contracted continuum gap `0.7245h`), **100% (`32/32`)** of outermost-shell particles are uniquely identified with zero active-set transitions during the central-difference probe, yet `three_moment_ghost_blend` still produces `last P/G = -0.3406` and shifts $R^*$ from $-16.48\%$ (`original`) to $-24.92\%$. The matched C4 `original` post-damped geometry reproduces the same effect at **100% (`32/32`)** identifiability (`tm_d = 0.0585h` vs `0.5511h`, `last P/G = -0.7133`, $R^* = -28.89\%$ vs $-17.25\%$ `original`). This confirms two independent conclusions:

- The inward last-shell force is an intrinsic smooth derivative of $\hat{\rho}_i(a_i(x))$, not an active-set artifact.
- Once a star has contracted, local cloud geometry reflects the *current* contracted particle boundary rather than the continuum hydrostatic radius $R_{\mathrm{ref}}$, so any local vacuum-ghost $H=0$ constraint locks zero enthalpy onto the contracted surface.

## 5. Decision and scope boundary

1. **Stage A geometry gate:** **FAIL** ($0.26\text{–}0.37\,h$ mean absolute depth error at $N=256$ vs $\le 0.15\,h$ requirement, driven by spherical cap curvature, nonlinear degenerate volume-ratio curvature, and discrete mass-quantile shell truncation).
2. **Stage B conservative screening:** **FAIL** across both orientation seeds (`0` and `17`) and both kernels (`M4` and `C4`); three-moment de-biasing without a ghost (`three_moment_debiased_affine_blend`) avoids the ghost-spring collapse and improves $R^*$ over `surface_mls` ($-3.35\%$ vs $-3.85\%$ M4; $-4.20\%$ vs $-5.44\%$ C4), but remains below `original` in both $R^*$ and last-shell pressure support.
3. **Stage C authorization:** **Not authorized.** No analytic derivative or $N=256$ dynamical hold is warranted for these rejected candidates. The validated BG-PHY-23 Hamiltonian, cold-electron EOS, Poisson-matched gravity, playable game physics, and [Issue #8](https://github.com/patrickdugan/Blighted-Galaxy/issues/8) acceptance gates remain unchanged.

### Reproducible commands and archived datasets

```sh
cargo test --lib sph_three_moment -- --nocapture
cargo test --lib sph_boundary_potential -- --nocapture
cargo run --release --bin phy26_bench -- --n 256 --seeds 0,17 --layout glass --stage initial --kernel both \
  --synthetic-csv bg-phy26-synthetic.csv --summary bg-phy26-summary.csv --shells-csv bg-phy26-shells.csv
```

Archived CSV datasets in this branch:
- [`workpapers/physics/data/bg-phy26/bg-phy26-synthetic.csv`](data/bg-phy26/bg-phy26-synthetic.csv)
- [`workpapers/physics/data/bg-phy26/bg-phy26-n64-seed0-summary.csv`](data/bg-phy26/bg-phy26-n64-seed0-summary.csv)
- [`workpapers/physics/data/bg-phy26/bg-phy26-n128-seed17-summary.csv`](data/bg-phy26/bg-phy26-n128-seed17-summary.csv)
- [`workpapers/physics/data/bg-phy26/bg-phy26-n256-initial-summary.csv`](data/bg-phy26/bg-phy26-n256-initial-summary.csv)
- [`workpapers/physics/data/bg-phy26/bg-phy26-n256-initial-shells.csv`](data/bg-phy26/bg-phy26-n256-initial-shells.csv)
- [`workpapers/physics/data/bg-phy26/bg-phy26-n256-seed0-post-summary.csv`](data/bg-phy26/bg-phy26-n256-seed0-post-summary.csv)
