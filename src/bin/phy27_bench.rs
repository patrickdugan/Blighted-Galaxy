//! BG-PHY-27: radially graded unequal-mass white-dwarf equilibrium experiment.
//!
//! Tests whether power-law radial mass grading `n_alpha(r) ~ rho(r)^(1-alpha)`
//! (`alpha = 0` equal mass through `alpha = 1` equal volume) resolves the
//! surface smoothing/tail contraction deficit under unmodified conservative
//! Price-Monaghan variable-h SPH + Poisson-matched gravity over a full
//! 3-second damped relaxation and 7.88-second (3.0 t_dyn) undamped hold.
//!
//! cargo run --release --bin phy27_bench -- --n 256 --seeds 0 --alphas 0.0,0.25,0.5,0.75,1.0 \
//!   --kernel both --relax-steps 300 --hold-steps 788 --dt 0.01 --csv bg-phy27-n256.csv
use blighted_galaxy::physics::{G_SI, SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self, Massive3, Sph3Config, V3};
use blighted_galaxy::sph_adaptive::{self, AdaptiveConfig, AdaptiveState};
use blighted_galaxy::sph_graded_seed;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::stellar_balance;
use blighted_galaxy::stellar_relax::{self, DampingSchedule, RelaxParameters};
use blighted_galaxy::white_dwarf::{profile_for_mass, DegenerateElectrons, WhiteDwarfProfile};
use std::fmt::Write as _;
use std::time::Instant;

fn arg<'a>(a: &'a [String], flag: &str) -> Option<&'a str> {
    a.windows(2).find(|x| x[0] == flag).map(|x| x[1].as_str())
}
fn opt<T: std::str::FromStr>(a: &[String], flag: &str, default: T) -> T {
    arg(a, flag).and_then(|v| v.parse().ok()).unwrap_or(default)
}
fn parse_list<T: std::str::FromStr>(raw: &str) -> Vec<T>
where
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    raw.split(',')
        .map(|s| s.trim().parse::<T>().expect("invalid list item"))
        .collect()
}
fn num(x: f64) -> String {
    format!("{x:.10e}")
}

fn kernel_config(star: &WhiteDwarfProfile, kernel: &str, eta_override: f64) -> AdaptiveConfig {
    let mut cfg = AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity()
        .with_newton_density();
    // Allow smaller surface smoothing lengths for graded-mass outer shells.
    cfg.min_h_m = 0.015 * star.radius_m;
    if kernel == "c4" {
        cfg = cfg.with_wendland_c4();
    }
    cfg.eta = if eta_override > 0.0 {
        eta_override
    } else if kernel == "m4" {
        1.2
    } else {
        1.7
    };
    cfg
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: usize = opt(&a, "--n", 256);
    let seeds: Vec<u64> = parse_list(arg(&a, "--seeds").unwrap_or("0"));
    let alphas: Vec<f64> = parse_list(arg(&a, "--alphas").unwrap_or("0.0,0.25,0.50,0.75,1.0"));
    let kernel = arg(&a, "--kernel").unwrap_or("both");
    let relax_steps: usize = opt(&a, "--relax-steps", 300);
    let hold_steps: usize = opt(&a, "--hold-steps", 788);
    let dt: f64 = opt(&a, "--dt", 0.01);
    let tau: f64 = opt(&a, "--damp-tau", 0.7);
    let glass_sweeps: usize = opt(&a, "--glass-sweeps", 20);
    let eta_override: f64 = opt(&a, "--eta", 0.0);
    assert!(n >= 32 && n % 8 == 0);
    assert!(["m4", "c4", "both"].contains(&kernel));

    let star = profile_for_mass(
        0.65 * SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),
        10_000.0,
    )
    .unwrap();
    let physics = Sph3Config::cold_white_dwarf(2.0).unwrap();
    let t_dyn = (star.radius_m.powi(3) / (G_SI * star.mass_kg)).sqrt();
    let kernels = if kernel == "both" {
        vec!["m4", "c4"]
    } else {
        vec![kernel]
    };

    let mut csv = String::from(
        "N,seed,kernel,eta,alpha,mass_ratio_max_min,max_adj_mass_ratio,init_rmax_bias_pct,init_last_p_over_g,init_cdf_rms,relax_s,post_rmax_km,post_rmax_bias_pct,post_r99_km,post_cdf_rms,post_virial,post_h_cap_frac,post_min_omega,hold_s,hold_tdyn,hold_max_rmax_excursion_pct,hold_max_r99_excursion_pct,hold_max_cdf_rms,hold_final_rmax_bias_pct,hold_energy_rel_drift,hold_momentum_rel_drift,issue8_radius_pass,issue8_hold_pass,issue8_cdf_pass,issue8_all_pass,wall_s\n",
    );

    let empty = Massive3 {
        pos_m: V3::default(),
        vel_m_s: V3::default(),
        mass_kg: 0.0,
        radius_m: 0.0,
    };

    for &seed in &seeds {
        let center = V3::new(1e8, 0.0, 0.0);
        for &alpha in &alphas {
            let (base, seed_rep) = sph_graded_seed::seed_graded_shells_oriented(
                &star,
                center,
                V3::default(),
                n,
                alpha,
                seed,
            );
            let p_init = if glass_sweeps > 0 {
                sph_graded_seed::relax_antipodal_angles_graded(
                    &base,
                    &star,
                    center,
                    1.2,
                    glass_sweeps,
                )
            } else {
                base
            };
            let init_diag = sph_massprofile::profile_diagnostics(&p_init, &star);
            for &k in &kernels {
                let started = Instant::now();
                let cfg = kernel_config(&star, k, eta_override);
                let mut p = p_init.clone();
                let mut state = AdaptiveState::new(&p, star.radius_m, cfg);
                let init_budget = stellar_balance::radial_balance(
                    &p,
                    &mut state,
                    physics,
                    cfg,
                    &star,
                    8,
                );
                let last_sh = init_budget.shells.last().unwrap();
                let init_last_pg =
                    last_sh.pressure_out_m_s2 / last_sh.hydrostatic_gravity_in_m_s2.max(1e-30);

                let relax_opts = RelaxParameters {
                    iterations: relax_steps,
                    dt_s: dt,
                    damping_time_s: tau,
                };
                let relax_rep = stellar_relax::relax_isolated_adaptive_with_schedule(
                    &mut p,
                    physics,
                    &star,
                    cfg,
                    &mut state,
                    relax_opts,
                    DampingSchedule::Constant { tau_s: tau },
                )
                .expect("graded relaxation succeeded");

                let post_diag = sph_massprofile::profile_diagnostics(&p, &star);
                let post_budget = stellar_balance::radial_balance(
                    &p,
                    &mut state,
                    physics,
                    cfg,
                    &star,
                    8,
                );

                // Conservative undamped hold, sampling every 10 steps for max excursion.
                let energy0 = sph_adaptive::adaptive_energy_j(&empty, &p, &mut state, physics, cfg);
                let mom0 = sph3::total_momentum(&empty, &p);
                let rmax0 = post_diag.radius_max_m;
                let r99_0 = post_diag.r99_m;
                let mut max_rmax_exc = 0.0_f64;
                let mut max_r99_exc = 0.0_f64;
                let mut max_cdf_rms = post_diag.radial_quantile_rms;
                let mut max_e_drift = 0.0_f64;

                let mut central = empty;
                for step in 1..=hold_steps {
                    sph_adaptive::advance(&mut central, &mut p, &mut state, physics, cfg, dt);
                    if step % 10 == 0 || step == hold_steps {
                        let d = sph_massprofile::profile_diagnostics(&p, &star);
                        max_rmax_exc = max_rmax_exc.max(((d.radius_max_m / rmax0) - 1.0).abs());
                        max_r99_exc = max_r99_exc.max(((d.r99_m / r99_0) - 1.0).abs());
                        max_cdf_rms = max_cdf_rms.max(d.radial_quantile_rms);
                        let e = sph_adaptive::energy_j_cached(&empty, &p, &state, physics, cfg);
                        max_e_drift = max_e_drift.max((e - energy0).abs() / energy0.abs().max(1.0));
                    }
                }
                let final_diag = sph_massprofile::profile_diagnostics(&p, &star);
                let mom1 = sph3::total_momentum(&empty, &p);
                let mom_scale = p.iter().map(|x| x.mass_kg * x.vel_m_s.len()).sum::<f64>().max(1.0);
                let mom_drift = (mom1 - mom0).len() / mom_scale;
                let hold_s = hold_steps as f64 * dt;
                let hold_tdyn = hold_s / t_dyn;
                let post_bias_pct = 100.0 * (post_diag.radius_max_m / star.radius_m - 1.0);
                let final_bias_pct = 100.0 * (final_diag.radius_max_m / star.radius_m - 1.0);

                let radius_pass = post_bias_pct.abs() <= 5.0;
                let hold_pass = hold_tdyn >= 3.0 && 100.0 * max_rmax_exc <= 5.0 && 100.0 * max_r99_exc <= 5.0;
                let cdf_pass = max_cdf_rms <= 0.05;
                let all_pass = radius_pass && hold_pass && cdf_pass && max_e_drift < 1e-3 && mom_drift < 1e-9;

                let row = [
                    n.to_string(),
                    seed.to_string(),
                    k.into(),
                    num(cfg.eta),
                    num(alpha),
                    num(seed_rep.mass_ratio_max_to_min),
                    num(seed_rep.max_adjacent_pair_mass_ratio),
                    num(seed_rep.initial_rmax_bias_pct),
                    num(init_last_pg),
                    num(init_diag.radial_quantile_rms),
                    num(relax_steps as f64 * dt),
                    num(post_diag.radius_max_m / 1000.0),
                    num(post_bias_pct),
                    num(post_diag.r99_m / 1000.0),
                    num(post_diag.radial_quantile_rms),
                    num(post_budget.pressure_to_gravity_virial),
                    num(relax_rep.capped_fraction),
                    num(relax_rep.minimum_omega),
                    num(hold_s),
                    num(hold_tdyn),
                    num(100.0 * max_rmax_exc),
                    num(100.0 * max_r99_exc),
                    num(max_cdf_rms),
                    num(final_bias_pct),
                    num(max_e_drift),
                    num(mom_drift),
                    radius_pass.to_string(),
                    hold_pass.to_string(),
                    cdf_pass.to_string(),
                    all_pass.to_string(),
                    num(started.elapsed().as_secs_f64()),
                ];
                let _ = writeln!(csv, "{}", row.join(","));
                println!(
                    "BG-PHY-27 N={n} s={seed} {k:2} eta={:.2} alpha={:.2} (mRatio={:.1}) initR={:+.2}% postR={:+.2}% holdExc(Rmax/R99)={:.2}%/{:.2}% maxCDF={:.4} dE={:.2e} GATE(R={},H={},C={},ALL={})",
                    cfg.eta,
                    alpha,
                    seed_rep.mass_ratio_max_to_min,
                    seed_rep.initial_rmax_bias_pct,
                    post_bias_pct,
                    100.0 * max_rmax_exc,
                    100.0 * max_r99_exc,
                    max_cdf_rms,
                    max_e_drift,
                    radius_pass,
                    hold_pass,
                    cdf_pass,
                    all_pass,
                );
            }
        }
    }

    if let Some(path) = arg(&a, "--csv") {
        std::fs::write(path, csv).expect("write phy27 csv");
        println!("wrote {path}");
    }
}
