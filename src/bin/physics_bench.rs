//! Run with: cargo run --release --bin physics_bench
//! Prints SI reference numbers and a deterministic one-year two-body orbital
//! regression; no graphics context or multiplayer/server is required.
use blighted_galaxy::physics::*;
use blighted_galaxy::V2;

fn main() {
    println!("Blighted Galaxy • SI physics benchmark");
    println!("NIST CODATA-2022: c = {:.0} m/s, G = {:.7e}", C_SI, G_SI);
    let rs = schwarzschild_radius_m(SOLAR_MASS_KG).unwrap();
    let th = hawking_temperature_k(SOLAR_MASS_KG).unwrap();
    println!("1 solar-mass Schwarzschild: radius {:.3} km, Hawking T {:.5e} K", rs/1000.0, th);

    let mass = 100_000.0;
    let t = hawking_lifetime_s(mass).unwrap();
    let energy = mass * C_SI * C_SI;
    println!("Idealized 100-ton micro-hole: lifetime {:.7} s, total rest-energy {:.4e} J", t, energy);
    println!("  Caution: constant species / semiclassical evaporation, no endpoint explosion claim");

    let impact = kinetic_energy_relativistic_j(2.0, 0.9993 * C_SI).unwrap();
    println!("2-kg RKM at 0.9993 c: kinetic energy {:.4e} J ({:.3e} Mt TNT)", impact, impact / 4.184e15);
    let tau = proper_time_constant_speed_s(1000.0 * JULIAN_YEAR_S, 0.7 * C_SI).unwrap();
    println!("1000 yr at 0.7c: traveler proper time {:.3} yr", tau / JULIAN_YEAR_S);
    let seen = retarded_emission_si(V2::new(LIGHT_YEAR_M, 0.0),
        V2::new(0.5*C_SI, 0.0), V2::default(), 0.0).unwrap();
    println!("Emitter receding at 0.5c from 1 ly today: retarded delay {:.6} yr",
        seen.delay_s / JULIAN_YEAR_S);

    let mut system = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
    let before = newtonian_energy_j(&system, 0.0);
    let period = kepler_period_s(SOLAR_MASS_KG + EARTH_MASS_KG, AU_M);
    let start = system[1].pos_m;
    evolve_nbody_si(&mut system, period, 6.0 * 3600.0, 0.0);
    let drift = (newtonian_energy_j(&system, 0.0) - before).abs() / before.abs();
    let closure = (system[1].pos_m - start).len() / AU_M;
    println!("Two-body 1-orbit test: energy drift {:.3e}, phase closure {:.3e} AU", drift, closure);
}
