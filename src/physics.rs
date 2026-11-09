//! SI-unit reference calculations for Blighted Galaxy.
//!
//! Units are explicitly SI: metres, kilograms, seconds, joules, watts and kelvin.
//! This module is intentionally distinct from the original normalized gameplay
//! sandbox in lib.rs. No general relativity, hydrodynamics, material strength,
//! Hawking greybody spectrum or quantum-gravity endpoint is simulated here.
//!
//! Constants: NIST 2022 CODATA (https://physics.nist.gov/cuu/Constants/).
//! Hawking lifetime model: uncharged nonrotating isolated hole, constant species
//! factor; e.g. https://www.physik.uzh.ch/groups/jetzer/notes/PHY529_ART_III.pdf
use crate::V2;
use std::f64::consts::PI;

/// 2022 CODATA gravitational constant, m^3 kg^-1 s^-2.
pub const G_SI: f64 = 6.674_30e-11;
/// Exact vacuum light speed, m/s.
pub const C_SI: f64 = 299_792_458.0;
/// Reduced Planck constant, J s (2022 CODATA).
pub const HBAR_SI: f64 = 1.054_571_817e-34;
/// Exact Boltzmann constant, J/K.
pub const KB_SI: f64 = 1.380_649e-23;
pub const AU_M: f64 = 149_597_870_700.0;
pub const JULIAN_YEAR_S: f64 = 365.25 * 86_400.0;
pub const LIGHT_YEAR_M: f64 = C_SI * JULIAN_YEAR_S;
pub const SOLAR_MASS_KG: f64 = 1.988_47e30;
pub const EARTH_MASS_KG: f64 = 5.972_2e24;

/// A finite Newtonian gravitating body in a shared barycentric coordinate frame.
/// Origin and axes must be consistent for all bodies in an integration.
#[derive(Clone, Copy, Debug)]
pub struct SiBody {
    pub mass_kg: f64,
    pub radius_m: f64,
    pub pos_m: V2,
    pub vel_m_s: V2,
}
impl SiBody {
    pub fn new(mass_kg: f64, radius_m: f64, pos_m: V2, vel_m_s: V2) -> Self {
        assert!(mass_kg.is_finite() && mass_kg > 0.0, "positive finite mass required");
        assert!(radius_m.is_finite() && radius_m >= 0.0, "nonnegative finite radius required");
        assert!(pos_m.x.is_finite() && pos_m.y.is_finite(), "finite position required");
        assert!(vel_m_s.x.is_finite() && vel_m_s.y.is_finite(), "finite velocity required");
        Self { mass_kg, radius_m, pos_m, vel_m_s }
    }
}

/// Pair-symmetric Newtonian acceleration with optional Plummer softening.
/// Forces satisfy Newton's third law to roundoff; r=0 and softening=0 are
/// singular, so overlapping coincident points are skipped (invalid input).
pub fn nbody_accelerations_si(bodies: &[SiBody], softening_m: f64) -> Vec<V2> {
    assert!(softening_m.is_finite() && softening_m >= 0.0);
    let mut acceleration = vec![V2::default(); bodies.len()];
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let delta = bodies[j].pos_m - bodies[i].pos_m;
            let r2 = delta.len_sq() + softening_m * softening_m;
            if r2 == 0.0 { continue; }
            let inverse_r3 = 1.0 / (r2 * r2.sqrt());
            let pull = delta * (G_SI * inverse_r3);
            acceleration[i] += pull * bodies[j].mass_kg;
            acceleration[j] -= pull * bodies[i].mass_kg;
        }
    }
    acceleration
}

/// Second-order, time-reversible velocity Verlet, fixed time step in seconds.
/// Choose a step small enough to resolve the shortest close-approach timescale.
/// Reject GR-dominated trajectories from this Newtonian approximation.
pub fn velocity_verlet_si(bodies: &mut [SiBody], dt_s: f64, softening_m: f64) {
    assert!(dt_s.is_finite() && dt_s >= 0.0, "nonnegative finite dt required");
    if dt_s == 0.0 || bodies.is_empty() { return; }
    let before = nbody_accelerations_si(bodies, softening_m);
    for (b, a) in bodies.iter_mut().zip(before) {
        b.pos_m += b.vel_m_s * dt_s + a * (0.5 * dt_s * dt_s);
        b.vel_m_s += a * (0.5 * dt_s);
    }
    let after = nbody_accelerations_si(bodies, softening_m);
    for (b, a) in bodies.iter_mut().zip(after) {
        b.vel_m_s += a * (0.5 * dt_s);
    }
}

/// Subdivide a desired interval using a maximum substep. This is deterministic
/// for a given input and controls integration error but does not adapt to shocks.
pub fn evolve_nbody_si(bodies: &mut [SiBody], interval_s: f64, max_step_s: f64, softening_m: f64) {
    assert!(interval_s.is_finite() && interval_s >= 0.0);
    assert!(max_step_s.is_finite() && max_step_s > 0.0);
    let steps = (interval_s / max_step_s).ceil().max(1.0) as usize;
    let step_s = interval_s / steps as f64;
    for _ in 0..steps {
        velocity_verlet_si(bodies, step_s, softening_m);
    }
}

/// Mechanical energy in joules for a Newtonian Plummer potential.
/// Appropriate only where Newtonian gravity and v << c hold.
pub fn newtonian_energy_j(bodies: &[SiBody], softening_m: f64) -> f64 {
    assert!(softening_m.is_finite() && softening_m >= 0.0);
    let kinetic: f64 = bodies.iter()
        .map(|b| 0.5 * b.mass_kg * b.vel_m_s.len_sq()).sum();
    let mut potential = 0.0;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let delta = bodies[i].pos_m - bodies[j].pos_m;
            let r = (delta.len_sq() + softening_m * softening_m).sqrt();
            if r > 0.0 { potential -= G_SI * bodies[i].mass_kg * bodies[j].mass_kg / r; }
        }
    }
    kinetic + potential
}

pub fn newtonian_momentum_kg_m_s(bodies: &[SiBody]) -> V2 {
    bodies.iter().fold(V2::default(), |sum, body|
        sum + body.vel_m_s * body.mass_kg)
}

/// Circular Newtonian two-body orbit with barycenter at origin.
/// Both body velocities are tangential; no central mass is artificially fixed.
pub fn circular_binary_si(m1_kg: f64, m2_kg: f64, separation_m: f64) -> [SiBody; 2] {
    assert!(m1_kg > 0.0 && m2_kg > 0.0 && separation_m > 0.0);
    let mass = m1_kg + m2_kg;
    let tangential = (G_SI * mass / separation_m).sqrt();
    let d1 = separation_m * m2_kg / mass;
    let d2 = separation_m * m1_kg / mass;
    [
        SiBody::new(m1_kg, 0.0, V2::new(-d1, 0.0),
                    V2::new(0.0, -tangential * m2_kg / mass)),
        SiBody::new(m2_kg, 0.0, V2::new(d2, 0.0),
                    V2::new(0.0, tangential * m1_kg / mass))
    ]
}

pub fn kepler_period_s(total_mass_kg: f64, separation_m: f64) -> f64 {
    assert!(total_mass_kg > 0.0 && separation_m > 0.0);
    2.0 * PI * (separation_m.powi(3) / (G_SI * total_mass_kg)).sqrt()
}

/// Contact velocity and energy models below are special-relativistic
/// *kinematics*, not a relativistic dynamical integrator.
pub fn lorentz_gamma(speed_m_s: f64) -> Option<f64> {
    let beta2 = (speed_m_s / C_SI).powi(2);
    if !speed_m_s.is_finite() || beta2 >= 1.0 { return None; }
    Some(1.0 / (1.0 - beta2).sqrt())
}

/// Stable computation of (gamma - 1) mc^2 avoiding cancellation when beta<<1.
pub fn kinetic_energy_relativistic_j(mass_kg: f64, speed_m_s: f64) -> Option<f64> {
    if !mass_kg.is_finite() || mass_kg < 0.0 || !speed_m_s.is_finite() { return None; }
    let beta2 = (speed_m_s / C_SI).powi(2);
    if beta2 >= 1.0 { return None; }
    let root = (1.0 - beta2).sqrt();
    let gamma_minus_one = beta2 / (root * (1.0 + root));
    Some(mass_kg * C_SI * C_SI * gamma_minus_one)
}
pub fn proper_time_constant_speed_s(coordinate_time_s: f64, speed_m_s: f64) -> Option<f64> {
    if !coordinate_time_s.is_finite() || coordinate_time_s < 0.0 { return None; }
    Some(coordinate_time_s / lorentz_gamma(speed_m_s)?)
}
pub fn light_travel_time_s(distance_m: f64) -> Option<f64> {
    if distance_m.is_finite() && distance_m >= 0.0 { Some(distance_m / C_SI) }
    else { None }
}

/// Collinear Einstein composition (u+v)/(1+u*v/c²), both signed.
pub fn add_collinear_velocities_m_s(u: f64, v: f64) -> Option<f64> {
    lorentz_gamma(u)?;
    lorentz_gamma(v)?;
    let velocity = (u + v) / (1.0 + u * v / (C_SI * C_SI));
    if velocity.abs() < C_SI { Some(velocity) } else { None }
}

#[derive(Debug, Clone, Copy)]
pub struct RetardedEmission {
    pub emission_time_s: f64,
    pub delay_s: f64,
    pub emission_pos_m: V2,
}
/// Position *seen* by a stationary observer at observer_time_s, for a source
/// moving at constant vector velocity from source_pos_at_zero_m.
/// Solve |source(t_emit)-observer| = c*(t_observer-t_emit).
pub fn retarded_emission_si(
    source_pos_at_zero_m: V2,
    source_velocity_m_s: V2,
    observer_pos_m: V2,
    observer_time_s: f64,
) -> Option<RetardedEmission> {
    if !observer_time_s.is_finite() ||
       !source_velocity_m_s.len_sq().is_finite() ||
       source_velocity_m_s.len() >= C_SI { return None; }
    let now = source_pos_at_zero_m + source_velocity_m_s * observer_time_s;
    let r = now - observer_pos_m;
    let c2_minus_v2 = C_SI * C_SI - source_velocity_m_s.len_sq();
    let rv = r.dot(source_velocity_m_s);
    let discriminant = rv * rv + c2_minus_v2 * r.len_sq();
    if !discriminant.is_finite() { return None; }
    let delay = (discriminant.sqrt() - rv) / c2_minus_v2;
    if delay < 0.0 || !delay.is_finite() { return None; }
    let emit_time = observer_time_s - delay;
    Some(RetardedEmission {
        emission_time_s: emit_time,
        delay_s: delay,
        emission_pos_m: source_pos_at_zero_m + source_velocity_m_s * emit_time,
    })
}

pub fn schwarzschild_radius_m(mass_kg: f64) -> Option<f64> {
    if !mass_kg.is_finite() || mass_kg <= 0.0 { return None; }
    Some(2.0 * G_SI * mass_kg / (C_SI * C_SI))
}
pub fn hawking_temperature_k(mass_kg: f64) -> Option<f64> {
    if !mass_kg.is_finite() || mass_kg <= 0.0 { return None; }
    Some(HBAR_SI * C_SI.powi(3) / (8.0 * PI * G_SI * mass_kg * KB_SI))
}

/// Constant-species blackbody approximation coefficient for dm/dt=-alpha/m².
/// Greybody, spin, charges and particle thresholds change the physical rate.
pub fn hawking_alpha_kg3_s() -> f64 {
    HBAR_SI * C_SI.powi(4) / (15360.0 * PI * G_SI.powi(2))
}
pub fn hawking_lifetime_s(mass_kg: f64) -> Option<f64> {
    if !mass_kg.is_finite() || mass_kg <= 0.0 { return None; }
    Some(mass_kg.powi(3) / (3.0 * hawking_alpha_kg3_s()))
}
pub fn hawking_power_w(mass_kg: f64) -> Option<f64> {
    if !mass_kg.is_finite() || mass_kg <= 0.0 { return None; }
    Some(hawking_alpha_kg3_s() * C_SI * C_SI / mass_kg.powi(2))
}
/// Analytic mass evolution under the constant-alpha approximation.
/// Returns 0 at the assumed evaporation endpoint; this endpoint is a *model
/// convention*, not a quantum-gravity prediction.
pub fn hawking_remaining_mass_kg(initial_mass_kg: f64, elapsed_s: f64) -> Option<f64> {
    if !initial_mass_kg.is_finite() || initial_mass_kg <= 0.0 ||
       !elapsed_s.is_finite() || elapsed_s < 0.0 { return None; }
    let fraction = elapsed_s / hawking_lifetime_s(initial_mass_kg)?;
    if fraction >= 1.0 { Some(0.0) }
    else { Some(initial_mass_kg * (1.0 - fraction).cbrt()) }
}
pub fn hawking_radiated_energy_j(initial_mass_kg: f64, elapsed_s: f64) -> Option<f64> {
    let remaining = hawking_remaining_mass_kg(initial_mass_kg, elapsed_s)?;
    Some((initial_mass_kg - remaining) * C_SI * C_SI)
}

/// Difference in Newtonian acceleration between a body's near-side surface and
/// its center along the line to an external attractor, divided by its own
/// surface gravitational acceleration. Defined only for non-overlapping bodies.
pub fn near_side_tidal_ratio_si(attractor: &SiBody, target: &SiBody) -> Option<f64> {
    if target.radius_m <= 0.0 || target.mass_kg <= 0.0 { return None; }
    let d = (target.pos_m - attractor.pos_m).len();
    if d <= target.radius_m + attractor.radius_m { return None; }
    let near_r = d - target.radius_m;
    let tidal = G_SI * attractor.mass_kg *
        (1.0 / (near_r * near_r) - 1.0 / (d * d));
    let self_gravity = G_SI * target.mass_kg / target.radius_m.powi(2);
    Some(tidal / self_gravity)
}

/// Approximate first-order vacuum deflection of light at distant impact b.
/// Reject near-horizon trajectories where the weak-field approximation fails.
pub fn weak_lensing_deflection_rad(mass_kg: f64, impact_parameter_m: f64) -> Option<f64> {
    let rs = schwarzschild_radius_m(mass_kg)?;
    if !impact_parameter_m.is_finite() || impact_parameter_m < 20.0 * rs { return None; }
    Some(4.0 * G_SI * mass_kg / (C_SI.powi(2) * impact_parameter_m))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(actual: f64, expected: f64, rel: f64) {
        let tolerance = rel * expected.abs().max(1e-280);
        assert!((actual - expected).abs() <= tolerance,
            "actual={actual:.13e} expected={expected:.13e} rel_err={:.6e}",
            (actual - expected).abs() / expected.abs().max(1e-280));
    }

    #[test]
    fn unit_constants_match_codatas_2022() {
        assert_eq!(C_SI, 299_792_458.0);
        assert_eq!(AU_M, 149_597_870_700.0);
        near(G_SI, 6.6743e-11, 1e-15);
    }

    #[test]
    fn solar_schwarzschild_radius_matches_2_95_km() {
        near(schwarzschild_radius_m(SOLAR_MASS_KG).unwrap(), 2953.34, 1e-5);
    }
    #[test]
    fn solar_hawking_temperature_matches_62_nanokelvin() {
        near(hawking_temperature_k(SOLAR_MASS_KG).unwrap(), 6.17e-8, 1e-3);
    }
    #[test]
    fn lifetime_scales_as_mass_cubed() {
        let a = hawking_lifetime_s(1.0e6).unwrap();
        let b = hawking_lifetime_s(2.0e6).unwrap();
        near(b / a, 8.0, 1e-14);
    }
    #[test]
    fn hawking_power_times_mass_squared_constant() {
        let p1 = hawking_power_w(1e9).unwrap();
        let p2 = hawking_power_w(2e9).unwrap();
        near(p1 / p2, 4.0, 1e-14);
        near(p1, 3.5616e14, 2e-5);
    }
    #[test]
    fn evaporation_analytic_partition_matches_lifetime() {
        let m = 1.0e5;
        let lifetime = hawking_lifetime_s(m).unwrap();
        near(hawking_remaining_mass_kg(m, lifetime * 0.875).unwrap(), m * 0.5, 1e-12);
        assert_eq!(hawking_remaining_mass_kg(m, lifetime * 1.01), Some(0.0));
        near(hawking_radiated_energy_j(m, lifetime).unwrap(), m * C_SI * C_SI, 1e-14);
    }
    #[test]
    fn lorentz_and_proper_time_at_point_eight_c() {
        let v = 0.8 * C_SI;
        near(lorentz_gamma(v).unwrap(), 5.0 / 3.0, 1e-14);
        near(proper_time_constant_speed_s(100.0, v).unwrap(), 60.0, 1e-14);
    }
    #[test]
    fn kinetic_energy_stays_nonzero_at_low_beta() {
        let v = 0.001 * C_SI;
        let calculated = kinetic_energy_relativistic_j(1.0, v).unwrap();
        near(calculated / (0.5 * v * v), 1.0, 1e-6);
    }
    #[test]
    fn relativistic_kinetic_energy_half_c() {
        let joules = kinetic_energy_relativistic_j(1.0, C_SI * 0.5).unwrap();
        near(joules, 1.39037910021726e16, 1e-12);
        assert!(kinetic_energy_relativistic_j(1.0, C_SI).is_none());
    }
    #[test]
    fn einstein_composition_stays_below_c() {
        let velocity = add_collinear_velocities_m_s(0.7 * C_SI, 0.7 * C_SI).unwrap();
        near(velocity / C_SI, 1.4 / 1.49, 1e-14);
        assert!(velocity < C_SI);
    }
    #[test]
    fn light_year_takes_julian_year() {
        near(light_travel_time_s(LIGHT_YEAR_M).unwrap(), JULIAN_YEAR_S, 1e-15);
    }
    #[test]
    fn retarded_stationary_emitter_delay_one_light_second() {
        let emitter = V2::new(C_SI, 0.0);
        let retarded = retarded_emission_si(emitter, V2::default(), V2::default(), 50.0).unwrap();
        near(retarded.delay_s, 1.0, 1e-14);
        near(retarded.emission_time_s, 49.0, 1e-14);
    }
    #[test]
    fn receding_emitter_retarded_delay_less_than_current_range_over_c() {
        let emitter = V2::new(C_SI, 0.0);
        let observed = retarded_emission_si(emitter, V2::new(0.5 * C_SI, 0.0),
            V2::default(), 0.0).unwrap();
        near(observed.delay_s, 2.0 / 3.0, 1e-14);
        near(observed.emission_pos_m.x, C_SI * (2.0 / 3.0), 1e-14);
    }
    #[test]
    fn near_side_tidal_acceleration_increases_at_smaller_distance() {
        let attractor = SiBody::new(SOLAR_MASS_KG, 6.96e8,
            V2::default(), V2::default());
        let target = SiBody::new(EARTH_MASS_KG, 6.371e6,
            V2::new(AU_M, 0.0), V2::default());
        let near_target = SiBody { pos_m: V2::new(AU_M / 2.0, 0.0), ..target };
        let far = near_side_tidal_ratio_si(&attractor, &target).unwrap();
        let close = near_side_tidal_ratio_si(&attractor, &near_target).unwrap();
        assert!(far > 0.0 && close > 7.9 * far);
    }
    #[test]
    fn lensing_weak_field_scales_inverse_impact() {
        let rs = schwarzschild_radius_m(SOLAR_MASS_KG).unwrap();
        let a = weak_lensing_deflection_rad(SOLAR_MASS_KG, 200.0 * rs).unwrap();
        let b = weak_lensing_deflection_rad(SOLAR_MASS_KG, 400.0 * rs).unwrap();
        near(a / b, 2.0, 1e-12);
        assert!(weak_lensing_deflection_rad(SOLAR_MASS_KG, 2.0 * rs).is_none());
    }
    #[test]
    fn nbody_newtons_third_law_pairwise() {
        let bodies = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
        let a = nbody_accelerations_si(&bodies, 0.0);
        let resultant = a[0] * bodies[0].mass_kg + a[1] * bodies[1].mass_kg;
        let fmag = a[0].len() * bodies[0].mass_kg;
        assert!(resultant.len() / fmag < 1e-12);
    }
    #[test]
    fn sun_earth_one_year_energy_and_orbit_are_stable() {
        let mut bodies = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
        let e0 = newtonian_energy_j(&bodies, 0.0);
        let p0 = newtonian_momentum_kg_m_s(&bodies);
        let start = bodies[1].pos_m;
        let period = kepler_period_s(SOLAR_MASS_KG + EARTH_MASS_KG, AU_M);
        evolve_nbody_si(&mut bodies, period, 3600.0 * 6.0, 0.0);
        let e1 = newtonian_energy_j(&bodies, 0.0);
        let p1 = newtonian_momentum_kg_m_s(&bodies);
        assert!((e1 - e0).abs() / e0.abs() < 1e-7);
        assert!((bodies[1].pos_m - start).len() / AU_M < 5e-4);
        let pscale = bodies[1].mass_kg * bodies[1].vel_m_s.len();
        assert!((p1 - p0).len() / pscale < 1e-10);
    }
    #[test]
    fn nbody_step_zero_preserves_state() {
        let mut bodies = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
        let start = bodies[1];
        velocity_verlet_si(&mut bodies, 0.0, 0.0);
        assert_eq!(start.pos_m, bodies[1].pos_m);
        assert_eq!(start.vel_m_s, bodies[1].vel_m_s);
    }
    #[test]
    fn invalid_physical_inputs_are_rejected() {
        assert!(schwarzschild_radius_m(0.0).is_none());
        assert!(hawking_power_w(-1.0).is_none());
        assert!(proper_time_constant_speed_s(10.0, C_SI).is_none());
        assert!(light_travel_time_s(-1.0).is_none());
    }
}
