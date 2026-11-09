//! BG-PHY-10: compare matched M4 gravitation and approximately cubic shells.
 //!
 //! cargo run --release --bin phy10_bench -- --csv bg-phy10-comparison.csv
 //! Each particle count has a common damping interval across variants.
 //! This is a measured A/B numerical experiment; the reference-star
 //! radius and free-evolution stability gates are not waived.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Sph3Config,V3};
use blighted_galaxy::sph_sampling::{
    self,StarSampling,seed_approximately_cubic_shells,
};
use blighted_galaxy::sph_adaptive::{
    AdaptiveConfig,AdaptiveState,GravityModel,
};
use blighted_galaxy::stellar_relax::{
    relax_isolated_adaptive,adaptive_undamped_hold,RelaxParameters,
};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::time::Instant;

fn name(method:StarSampling)->&'static str {
    match method {
        StarSampling::LegacyAntipodalFibonacci=>"legacy",
        StarSampling::ShellBalancedAntipodal=>"radius_shells",
        StarSampling::ApproximatelyCubicShells=>"nearly_cubic_shells",
    }
}
fn gravity_name(m:GravityModel)->&'static str {
    match m {
        GravityModel::FixedPlummer=>"fixed_plummer",
        GravityModel::MatchedAdaptivePlummer=>"variable_plummer",
        GravityModel::MatchedCubicSpline=>"poisson_cubic",
    }
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let path=args.windows(2).find(|a|a[0]=="--csv").map(|a|a[1].as_str());
    let quick=args.iter().any(|arg|arg=="--quick");
    let extended=args.iter().any(|arg|arg=="--extended-hold");
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    println!("BG-PHY-10 M4 POISSON GRAVITY / NEAR-CUBIC HYDROSTATIC SHELLS");
    println!("Reference star: 0.65 Msun, R={:.3} km; t_dyn={:.5} s",
        star.radius_m/1e3,t_dyn);
    println!("NOTE: all variants at each N share a damping interval; intervals vary by N.");
    let mut csv=String::from("N,method,gravity,relax_s,shells,median_shell_aspect,force_ratio,logrho_initial,logrho_final,radius_km,radius_bias_pct,omega_min,cap_fraction,wall_s,hold_s,hold_radius_change_pct,hold_energy_drift,hold_momentum_drift\n");
    let ns:Vec<usize>=if quick {vec![64,128]} else {vec![64,128,256,512,1024]};
    for n in ns {
        let iter=if n<=128 {300} else if n<=512 {120} else {60};
        let variants:Vec<(StarSampling,GravityModel)>=if n<=512 {
            vec![
                (StarSampling::LegacyAntipodalFibonacci,GravityModel::MatchedAdaptivePlummer),
                (StarSampling::ShellBalancedAntipodal,GravityModel::MatchedAdaptivePlummer),
                (StarSampling::ApproximatelyCubicShells,GravityModel::MatchedAdaptivePlummer),
                (StarSampling::ApproximatelyCubicShells,GravityModel::MatchedCubicSpline),
            ]
        } else {
            vec![
                (StarSampling::LegacyAntipodalFibonacci,GravityModel::MatchedAdaptivePlummer),
                (StarSampling::ApproximatelyCubicShells,GravityModel::MatchedAdaptivePlummer),
                (StarSampling::ApproximatelyCubicShells,GravityModel::MatchedCubicSpline),
            ]
        };
        let options=RelaxParameters {
            iterations:iter,dt_s:0.01,damping_time_s:0.7,
        };
        for (sampling,gravity) in variants {
            let center=V3::new(1e8,0.0,0.0);
            let mut particles=sph_sampling::seeded_star(
                &star,center,V3::default(),n,sampling);
            let report=if let StarSampling::ApproximatelyCubicShells=sampling {
                Some(seed_approximately_cubic_shells(&star,
                    center,V3::default(),n).1)
            } else {None};
            let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m);
            cfg=match gravity {
                GravityModel::FixedPlummer=>cfg,
                GravityModel::MatchedAdaptivePlummer=>cfg.with_adaptive_gravity(),
                GravityModel::MatchedCubicSpline=>cfg.with_cubic_gravity(),
            };
            let mut state=AdaptiveState::new(
                &particles,star.radius_m,cfg);
            let start=Instant::now();
            let relaxed=relax_isolated_adaptive(
                &mut particles,physics,&star,cfg,&mut state,options)
                .expect("relaxation must remain finite");
            let wall=start.elapsed().as_secs_f64();
            let radius_km=relaxed.report.max_radius_m/1000.0;
            let bias=100.0*relaxed.radius_bias;
            let shells=report.as_ref().map(|r|r.shells.len()).unwrap_or(0);
            let aspect=report.as_ref().map(|r|r.median_tile_aspect)
                .unwrap_or(f64::NAN);
            println!("N={n:4} {:20} {:16} R={radius_km:8.1} km bias={bias:+7.2}% force={:.4} rhoLog={:.4} wall={wall:.2}s shells={shells}",
                name(sampling),gravity_name(gravity),
                relaxed.report.force_reduction(),
                relaxed.report.final_density_log_rms);

            let mut hold_s=f64::NAN;
            let mut hold_radius=f64::NAN;
            let mut hold_energy=f64::NAN;
            let mut hold_momentum=f64::NAN;
            let run_hold=gravity==GravityModel::MatchedCubicSpline &&
                (n==128||n==256||(extended&&n==512));
            if run_hold {
                let step=0.01;
                let steps=(3.0*t_dyn/step).ceil() as usize;
                let hold=adaptive_undamped_hold(
                    &mut particles,physics,cfg,&mut state,step,steps)
                    .expect("undamped energy and momentum finite");
                hold_s=hold.elapsed_s;
                hold_radius=100.0*hold.radius_relative_change;
                hold_energy=hold.energy_relative_drift;
                hold_momentum=hold.momentum_relative_drift;
                println!("  undamped hold 3 t_dyn={hold_s:.3}s: dR={hold_radius:+.3}%, dE/E={hold_energy:.3e}, dP={hold_momentum:.3e}");
                assert!(hold.mass_relative_drift<1e-12);
                assert!(hold.momentum_relative_drift<1e-8);
            }
            assert!(relaxed.report.mass_relative_error<1e-12);
            assert!(relaxed.report.center_shift_m/star.radius_m<1e-8);
            assert!(relaxed.report.final_density_log_rms.is_finite());
            assert!(relaxed.minimum_omega>0.0);
            csv.push_str(&format!("{n},{},{},{:.3},{shells},{:.8},{:.8},{:.8},{:.8},{radius_km:.5},{bias:.6},{:.8},{:.7},{wall:.5},{hold_s:.5},{hold_radius:.8},{hold_energy:.9},{hold_momentum:.9}\n",
                name(sampling),gravity_name(gravity),
                iter as f64*options.dt_s,aspect,
                relaxed.report.force_reduction(),
                relaxed.report.initial_density_log_rms,
                relaxed.report.final_density_log_rms,
                relaxed.minimum_omega,relaxed.capped_fraction,
            ));
        }
    }
    if let Some(p)=path {
        std::fs::write(p,csv).expect("write BG-PHY-10 comparison");
        println!("wrote {p}");
    }
}
