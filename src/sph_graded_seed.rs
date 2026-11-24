//! BG-PHY-27: radially graded unequal-mass antipodal shell + angular-glass
//! white-dwarf initializers for Price-Monaghan variable-h SPH.
//!
//! Equal-mass shell sampling (`alpha = 0`) places the outermost N=256 shell at
//! r_max = 7840.7 km (-6.68% inside the 8402.207 km cold-electron reference)
//! with wide surface smoothing lengths `h_i` that average over dense interior
//! shells. Power-law number-density grading `n_alpha(r) ~ rho(r)^(1 - alpha)`
//! assigns exact continuum shell masses `m_i = Delta M_i / 2` to antipodal
//! pairs while increasing surface radial resolution as `alpha -> 1` (equal
//! volume). All conservative Price-Monaghan variable-h pressure and
//! Poisson-matched gravity equations remain unmodified.

use crate::sph3::{Fluid3, V3};
use crate::white_dwarf::WhiteDwarfProfile;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct GradedShell {
    pub first_pair: usize,
    pub pair_count: usize,
    pub lower_radius_m: f64,
    pub upper_radius_m: f64,
    pub thickness_over_tile: f64,
}

#[derive(Debug, Clone)]
pub struct GradedSeedReport {
    pub alpha: f64,
    pub shells: Vec<GradedShell>,
    pub radius_min_m: f64,
    pub radius_max_m: f64,
    pub initial_rmax_bias_pct: f64,
    pub mass_ratio_max_to_min: f64,
    pub max_adjacent_pair_mass_ratio: f64,
    pub mass_relative_error: f64,
    pub center_error_m: f64,
}

struct GradedCdf {
    radii_m: Vec<f64>,
    cum_frac: Vec<f64>,
}

impl GradedCdf {
    fn build(profile: &WhiteDwarfProfile, alpha: f64) -> Self {
        assert!((0.0..=1.0).contains(&alpha) && alpha.is_finite());
        let n = profile.shells.len();
        assert!(n >= 4);
        let mut radii_m = Vec::with_capacity(n + 1);
        let mut cum = Vec::with_capacity(n + 1);
        radii_m.push(0.0);
        cum.push(0.0);
        let mut prev_r = 0.0_f64;
        let mut prev_rho = profile.rho_c_kg_m3;
        let mut total = 0.0_f64;
        let power = 1.0 - alpha;
        for s in &profile.shells {
            let r = s.radius_m;
            let rho = s.rho_kg_m3.max(0.0);
            if r > prev_r {
                let r_mid = 0.5 * (prev_r + r);
                let rho_mid = 0.5 * (prev_rho + rho);
                let weight = if power == 0.0 {
                    1.0
                } else if rho_mid > 0.0 {
                    rho_mid.powf(power)
                } else {
                    0.0
                };
                // Exact spherical shell volume (4/3)*pi*(r^3 - prev_r^3) times density^power.
                let dv = (4.0 * PI / 3.0) * (r.powi(3) - prev_r.powi(3));
                let _ = r_mid;
                total += dv * weight;
            }
            radii_m.push(r);
            cum.push(total);
            prev_r = r;
            prev_rho = rho;
        }
        assert!(total > 0.0 && total.is_finite());
        for v in &mut cum {
            *v /= total;
        }
        *cum.last_mut().unwrap() = 1.0;
        Self {
            radii_m,
            cum_frac: cum,
        }
    }

    fn radius_at_fraction(&self, f: f64) -> f64 {
        let f = f.clamp(0.0, 1.0);
        if f <= 0.0 {
            return 0.0;
        }
        if f >= 1.0 {
            return *self.radii_m.last().unwrap();
        }
        let idx = self.cum_frac.partition_point(|&c| c < f);
        if idx == 0 {
            return self.radii_m[0];
        }
        if idx >= self.cum_frac.len() {
            return *self.radii_m.last().unwrap();
        }
        let (c0, c1) = (self.cum_frac[idx - 1], self.cum_frac[idx]);
        let (r0, r1) = (self.radii_m[idx - 1], self.radii_m[idx]);
        let t = if c1 > c0 {
            ((f - c0) / (c1 - c0)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        r0 + t * (r1 - r0)
    }
}

fn radical_inverse(mut n: usize, base: usize) -> f64 {
    let mut f = 1.0 / base as f64;
    let mut result = 0.0;
    while n > 0 {
        result += (n % base) as f64 * f;
        n /= base;
        f /= base as f64;
    }
    result
}

fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn random_unit(x: &mut u64) -> f64 {
    ((splitmix64(x) >> 11) as f64) / ((1_u64 << 53) as f64)
}

fn rotate_shell_vector(
    v: V3,
    cosa: f64,
    sina: f64,
    cosb: f64,
    sinb: f64,
    cosc: f64,
    sinc: f64,
) -> V3 {
    let x1 = cosa * v.x - sina * v.y;
    let y1 = sina * v.x + cosa * v.y;
    let x2 = cosb * x1 + sinb * v.z;
    let z2 = -sinb * x1 + cosb * v.z;
    V3::new(cosc * x2 - sinc * y1, sinc * x2 + cosc * y1, z2)
}

/// Build an antipodal approximately-cubic graded-mass white-dwarf seed.
/// At `alpha = 0.0` this reproduces equal-mass shell sampling; at `alpha = 1.0`
/// it produces equal-volume shells extending close to `profile.radius_m`.
pub fn seed_graded_shells_oriented(
    profile: &WhiteDwarfProfile,
    center: V3,
    velocity: V3,
    count: usize,
    alpha: f64,
    orientation_seed: u64,
) -> (Vec<Fluid3>, GradedSeedReport) {
    assert!(count >= 8 && count % 2 == 0);
    assert!((0.0..=1.0).contains(&alpha) && alpha.is_finite());
    let cdf = GradedCdf::build(profile, alpha);
    let pairs = count / 2;
    let golden = PI * (3.0 - 5.0_f64.sqrt());
    let mut ordered = Vec::<(usize, Fluid3, Fluid3)>::with_capacity(pairs);
    let mut shells = Vec::<GradedShell>::new();
    let mut first = 0_usize;

    while first < pairs {
        let remaining = pairs - first;
        let lower = cdf.radius_at_fraction(first as f64 / pairs as f64);
        let maximum = remaining.min(96);
        let mut best = (f64::INFINITY, 1_usize, 0.0, 0.0, 0.0);
        for k in 1..=maximum {
            if k == 1 && remaining > 1 {
                continue;
            }
            if remaining - k == 1 {
                continue;
            }
            let upper = cdf.radius_at_fraction((first + k) as f64 / pairs as f64);
            let middle = cdf.radius_at_fraction((first as f64 + 0.5 * k as f64) / pairs as f64);
            let thickness = upper - lower;
            let tile = middle * (2.0 * PI / k as f64).sqrt();
            let ratio = thickness / tile.max(1.0);
            let cost = ratio.ln().abs() + if k == 2 { 0.05 } else { 0.0 };
            if cost < best.0 {
                best = (cost, k, upper, middle, ratio);
            }
        }
        let (_, k, upper, _, ratio) = best;
        let shell_index = shells.len();
        shells.push(GradedShell {
            first_pair: first,
            pair_count: k,
            lower_radius_m: lower,
            upper_radius_m: upper,
            thickness_over_tile: ratio,
        });
        let phase = 2.0 * PI * radical_inverse(shell_index + 1, 3);
        let tilt = PI * radical_inverse(shell_index + 1, 5);
        let (sp, cp) = phase.sin_cos();
        let (st, ct) = tilt.sin_cos();

        for local in 0..k {
            let i = first + local;
            let r_lo = cdf.radius_at_fraction(i as f64 / pairs as f64);
            let r_hi = cdf.radius_at_fraction((i + 1) as f64 / pairs as f64);
            let m_lo = profile.enclosed_mass_at(r_lo);
            let m_hi = profile.enclosed_mass_at(r_hi);
            let pair_mass = (m_hi - m_lo).max(1e-30 * profile.mass_kg);
            let each_mass = 0.5 * pair_mass;
            // Place pair at the exact enclosed-mass midpoint of its interval [r_lo, r_hi].
            let f_mass_mid = ((m_lo + 0.5 * pair_mass) / profile.mass_kg).clamp(0.0, 1.0);
            let radius = profile.radius_at_mass_fraction(f_mass_mid);

            let z = 1.0 - 2.0 * (local as f64 + 0.5) / k as f64;
            let phi = golden * local as f64;
            let xy = (1.0 - z * z).sqrt();
            let nx = xy * phi.cos();
            let ny = xy * phi.sin();
            let rx = ct * nx + st * z;
            let rz = -st * nx + ct * z;
            let dir = V3::new(cp * rx - sp * ny, sp * rx + cp * ny, rz);
            let d = dir * radius;
            ordered.push((
                i,
                Fluid3 {
                    pos_m: center + d,
                    vel_m_s: velocity,
                    mass_kg: each_mass,
                },
                Fluid3 {
                    pos_m: center - d,
                    vel_m_s: velocity,
                    mass_kg: each_mass,
                },
            ));
        }
        first += k;
    }

    ordered.sort_by_key(|(i, _, _)| *i);
    let mut particles = Vec::with_capacity(count);
    for (_, a, b) in ordered {
        particles.push(a);
        particles.push(b);
    }

    // Normalize tiny quadrature residual to exact profile.mass_kg.
    let raw_total: f64 = particles.iter().map(|p| p.mass_kg).sum();
    let scale_m = profile.mass_kg / raw_total;
    for p in &mut particles {
        p.mass_kg *= scale_m;
    }

    if orientation_seed != 0 {
        let mut generator = orientation_seed;
        for shell in &shells {
            let a = 2.0 * PI * random_unit(&mut generator);
            let b = (1.0 - 2.0 * random_unit(&mut generator)).acos();
            let c = 2.0 * PI * random_unit(&mut generator);
            let (sa, ca) = a.sin_cos();
            let (sb, cb) = b.sin_cos();
            let (sc, cc) = c.sin_cos();
            for pair in shell.first_pair..shell.first_pair + shell.pair_count {
                let offset = particles[2 * pair].pos_m - center;
                let rotated = rotate_shell_vector(offset, ca, sa, cb, sb, cc, sc);
                particles[2 * pair].pos_m = center + rotated;
                particles[2 * pair + 1].pos_m = center - rotated;
            }
        }
    }

    let total: f64 = particles.iter().map(|p| p.mass_kg).sum();
    let mean = particles
        .iter()
        .fold(V3::default(), |sum, p| sum + p.pos_m * p.mass_kg)
        / total;
    let rmax = particles
        .iter()
        .map(|p| (p.pos_m - center).len())
        .fold(0.0_f64, f64::max);
    let rmin = particles
        .iter()
        .map(|p| (p.pos_m - center).len())
        .fold(f64::INFINITY, f64::min);
    let m_max = particles
        .iter()
        .map(|p| p.mass_kg)
        .fold(0.0_f64, f64::max);
    let m_min = particles
        .iter()
        .map(|p| p.mass_kg)
        .fold(f64::INFINITY, f64::min);
    let mut max_adj = 1.0_f64;
    for pair in 0..pairs - 1 {
        let ma = particles[2 * pair].mass_kg;
        let mb = particles[2 * (pair + 1)].mass_kg;
        max_adj = max_adj.max((ma / mb).max(mb / ma));
    }

    let report = GradedSeedReport {
        alpha,
        shells,
        radius_min_m: rmin,
        radius_max_m: rmax,
        initial_rmax_bias_pct: 100.0 * (rmax / profile.radius_m - 1.0),
        mass_ratio_max_to_min: m_max / m_min.max(1e-99),
        max_adjacent_pair_mass_ratio: max_adj,
        mass_relative_error: (total - profile.mass_kg).abs() / profile.mass_kg,
        center_error_m: (mean - center).len(),
    };
    (particles, report)
}

fn pair_energy_and_gradient(directions: &[V3], radii: &[f64], hs: &[f64]) -> (f64, Vec<V3>) {
    let count = directions.len();
    let mut grads = vec![V3::default(); count];
    let mut objective = 0.0;
    for i in 0..count {
        for j in (i + 1)..count {
            let ri = radii[i];
            let rj = radii[j];
            let dot = directions[i].dot(directions[j]).clamp(-1.0, 1.0);
            let sigma = 0.5 * (hs[i] + hs[j]);
            let sigma_sq = sigma * sigma;
            let d0 = ri * ri + rj * rj - 2.0 * ri * rj * dot;
            let d1 = ri * ri + rj * rj + 2.0 * ri * rj * dot;
            let same = (-d0 / (2.0 * sigma_sq)).exp();
            let flipped = (-d1 / (2.0 * sigma_sq)).exp();
            objective += same + flipped;
            let g = (ri * rj / sigma_sq) * (same - flipped);
            grads[i] += directions[j] * g;
            grads[j] += directions[i] * g;
        }
    }
    (objective, grads)
}

/// Unequal-mass-aware antipodal angular-glass relaxation using each pair's own
/// local smoothing scale `h_i = eta * (m_i / rho(r_i))^(1/3)`. Preserves every
/// particle's exact radius, mass, antipodal symmetry, and velocity.
pub fn relax_antipodal_angles_graded(
    source: &[Fluid3],
    profile: &WhiteDwarfProfile,
    center: V3,
    eta: f64,
    requested_sweeps: usize,
) -> Vec<Fluid3> {
    assert!(source.len() >= 8 && source.len() % 2 == 0);
    let count = source.len() / 2;
    let mut initial = Vec::with_capacity(count);
    let mut radii = Vec::with_capacity(count);
    let mut hs = Vec::with_capacity(count);
    for pair in source.chunks_exact(2) {
        let v = pair[0].pos_m - center;
        let r = v.len();
        assert!(r > 0.0 && r < profile.radius_m);
        assert!((pair[1].pos_m - center + v).len() < 1e-7 * profile.radius_m);
        assert!((pair[0].mass_kg - pair[1].mass_kg).abs() / pair[0].mass_kg < 1e-12);
        initial.push(v * (1.0 / r));
        radii.push(r);
        let rho = profile.density_at(r).max(profile.rho_c_kg_m3 * 1e-8);
        hs.push(
            (eta * (pair[0].mass_kg / rho).cbrt())
                .clamp(0.015 * profile.radius_m, 0.90 * profile.radius_m),
        );
    }
    let mut directions = initial.clone();
    let (mut objective, _) = pair_energy_and_gradient(&directions, &radii, &hs);
    for _ in 0..requested_sweeps {
        let (_, gradients) = pair_energy_and_gradient(&directions, &radii, &hs);
        let tangents: Vec<V3> = directions
            .iter()
            .zip(&gradients)
            .map(|(u, g)| *g - *u * g.dot(*u))
            .collect();
        let mut step = 0.24;
        let mut next = Vec::with_capacity(count);
        let mut accepted = false;
        for _ in 0..12 {
            next.clear();
            for (u, g) in directions.iter().zip(&tangents) {
                let scaled = step / (1.0 + g.len());
                next.push((*u - *g * scaled).normalized());
            }
            let (trial, _) = pair_energy_and_gradient(&next, &radii, &hs);
            if trial + 1e-13 * objective.abs().max(1.0) < objective {
                objective = trial;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if !accepted {
            break;
        }
        std::mem::swap(&mut next, &mut directions);
    }
    let mut p = source.to_vec();
    for i in 0..count {
        let offset = directions[i] * radii[i];
        p[2 * i].pos_m = center + offset;
        p[2 * i + 1].pos_m = center - offset;
    }
    p
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph3::{Massive3, Sph3Config};
    use crate::sph_adaptive::{self, AdaptiveConfig, AdaptiveState};
    use crate::sph_massprofile;
    use crate::white_dwarf::{profile_for_mass, DegenerateElectrons};

    fn wd() -> WhiteDwarfProfile {
        profile_for_mass(
            0.65 * SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),
            10_000.0,
        )
        .unwrap()
    }

    #[test]
    fn graded_seeds_preserve_mass_center_antipodes_and_extend_surface_radius() {
        let star = wd();
        let center = V3::new(1e8, -2e7, 5e6);
        let mut prev_rmax = 0.0_f64;
        for alpha in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let (p, rep) = seed_graded_shells_oriented(&star, center, V3::default(), 256, alpha, 17);
            assert!(rep.mass_relative_error < 1e-13);
            assert!(rep.center_error_m < 1e-5);
            assert!(rep.radius_max_m > prev_rmax, "alpha={alpha} should push rmax outward");
            prev_rmax = rep.radius_max_m;
            let diag = sph_massprofile::profile_diagnostics(&p, &star);
            assert!(diag.radial_quantile_rms < 0.03, "alpha={alpha} cdf rms {}", diag.radial_quantile_rms);
            let glass = relax_antipodal_angles_graded(&p, &star, center, 1.2, 12);
            for (a, b) in p.iter().zip(&glass) {
                assert_eq!(a.mass_kg, b.mass_kg);
                assert!(((a.pos_m - center).len() - (b.pos_m - center).len()).abs() < 1e-5);
            }
        }
        assert!(prev_rmax / star.radius_m > 0.985, "equal-volume seed should reach >98.5% R_ref");
    }

    #[test]
    fn unequal_mass_price_monaghan_force_matches_exact_energy_gradient() {
        let star = wd();
        let center = V3::new(1e8, 0.0, 0.0);
        let (p, _) = seed_graded_shells_oriented(&star, center, V3::default(), 32, 0.6, 17);
        let phys = Sph3Config::cold_white_dwarf(2.0).unwrap();
        let empty = Massive3 {
            pos_m: V3::default(),
            vel_m_s: V3::default(),
            mass_kg: 0.0,
            radius_m: 0.0,
        };
        for (k, eta) in [(false, 1.2), (true, 1.7)] {
            let mut cfg = AdaptiveConfig::from_star_radius(star.radius_m)
                .with_cubic_gravity()
                .with_newton_density();
            cfg.min_h_m = 0.015 * star.radius_m;
            if k {
                cfg = cfg.with_wendland_c4();
            }
            cfg.eta = eta;
            let mut state = AdaptiveState::new(&p, star.radius_m, cfg);
            let f = sph_adaptive::forces(&empty, &p, &mut state, phys, cfg);
            // Total linear momentum conservation for unequal masses: sum m_i a_i == 0.
            let net_f = p
                .iter()
                .zip(&f.accelerations_m_s2)
                .fold(V3::default(), |s, (x, a)| s + *a * x.mass_kg);
            let scale_f = p
                .iter()
                .zip(&f.accelerations_m_s2)
                .map(|(x, a)| x.mass_kg * a.len())
                .sum::<f64>();
            assert!(net_f.len() / scale_f < 1e-13, "momentum violation: {}", net_f.len() / scale_f);

            // Exact variational check: F_i = -dE/dx_i across core, middle, and surface particles.
            let eps = 1e-4 * star.radius_m;
            for idx in [0usize, 14, 30] {
                let mut fd_force = V3::default();
                for axis in 0..3 {
                    let mut plus = p.clone();
                    let mut minus = p.clone();
                    let delta = match axis {
                        0 => V3::new(eps, 0.0, 0.0),
                        1 => V3::new(0.0, eps, 0.0),
                        _ => V3::new(0.0, 0.0, eps),
                    };
                    plus[idx].pos_m += delta;
                    minus[idx].pos_m -= delta;
                    let mut sp = AdaptiveState::new(&plus, star.radius_m, cfg);
                    let mut sm = AdaptiveState::new(&minus, star.radius_m, cfg);
                    let ep = sph_adaptive::adaptive_energy_j(&empty, &plus, &mut sp, phys, cfg);
                    let em = sph_adaptive::adaptive_energy_j(&empty, &minus, &mut sm, phys, cfg);
                    let comp = -(ep - em) / (2.0 * eps);
                    match axis {
                        0 => fd_force.x = comp,
                        1 => fd_force.y = comp,
                        _ => fd_force.z = comp,
                    }
                }
                let ana = f.accelerations_m_s2[idx] * p[idx].mass_kg;
                let rel = (ana - fd_force).len() / ana.len().max(fd_force.len()).max(1.0);
                assert!(
                    rel < 5e-5,
                    "unequal-mass variational mismatch (c4={k}, idx={idx}): rel={rel}"
                );
            }
        }
    }
}
