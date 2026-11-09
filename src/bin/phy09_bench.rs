//! BG-PHY-09 reference: isolated cold white dwarf with three A/B variants.
//! cargo run --release --bin phy09_bench -- --csv bg-phy09-comparison.csv
//! This benchmark reports scientific acceptance gates without forcing a pass
//! when the star fails to hold the hydrostatic mass-radius profile.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_sampling::{
    self,StarSampling,
};
use blighted_galaxy::sph_adaptive::{
    AdaptiveConfig,AdaptiveState,GravityModel,
};
use blighted_galaxy::stellar_relax::{
    relax_isolated_adaptive,adaptive_undamped_hold,RelaxParameters,
};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::time::Instant;
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let csv_path=args.windows(2).find(|p|p[0]=="--csv")
        .map(|p|p[1].as_str());
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let reference=star.radius_m;
    println!("BG-PHY-09: variational adaptive Plummer gravity and radial-shell seeds");
    println!("hydrostatic 0.65 Msun radius {:.3} km, t_dyn {:.3} s",
        reference/1e3,t_dyn);
    let mut csv=String::from(
        "particles,seed,gravity,damped_seconds,force_ratio,density_log_rms,bias_percent,radius_km,minimum_omega,capped_fraction,elapsed_wall_s,hold_seconds,hold_radius_percent,hold_energy_drift,hold_momentum_drift\n"
    );
    for n in [64_usize,128,256] {
        let iterations=if n<=128 {300} else {120};
        let options=RelaxParameters{
            iterations,dt_s:0.01,damping_time_s:0.7,
        };
        for (seed,gravity) in [
            (StarSampling::LegacyAntipodalFibonacci,GravityModel::FixedPlummer),
            (StarSampling::ShellBalancedAntipodal,GravityModel::FixedPlummer),
            (StarSampling::ShellBalancedAntipodal,GravityModel::MatchedAdaptivePlummer),
        ] {
            let seed_label=match seed {
                StarSampling::LegacyAntipodalFibonacci=>"legacy",
                StarSampling::ShellBalancedAntipodal=>"shell",
            };
            let gravity_label=match gravity {
                GravityModel::FixedPlummer=>"fixed",
                GravityModel::MatchedAdaptivePlummer=>"variational",
            };
            let mut particles=sph_sampling::seeded_star(&star,
                V3::new(1e8,0.0,0.0),V3::default(),n,seed);
            let cfg=if gravity==GravityModel::MatchedAdaptivePlummer {
                AdaptiveConfig::from_star_radius(reference)
                    .with_adaptive_gravity()
            } else {AdaptiveConfig::from_star_radius(reference)};
            let mut state=AdaptiveState::new(&particles,reference,cfg);
            let timer=Instant::now();
            let result=relax_isolated_adaptive(
                &mut particles,physics,&star,cfg,&mut state,options)
                .expect("finite cold white-dwarf relaxation");
            let elapsed=timer.elapsed().as_secs_f64();
            println!("N={n:3} {seed_label:6} {gravity_label:11}: R={:.2} km ({:+.2}%), force={:.4}, log-rho={:.4}, minOmega={:.4}, elapsed={:.3}s",
                result.report.max_radius_m/1e3,
                result.radius_bias*100.0,
                result.report.force_reduction(),
                result.report.final_density_log_rms,
                result.minimum_omega,elapsed);
            let mut hold_s=f64::NAN;
            let mut hold_radius=f64::NAN;
            let mut hold_energy=f64::NAN;
            let mut hold_momentum=f64::NAN;
            if n<=128 && gravity==GravityModel::MatchedAdaptivePlummer {
                let steps=(3.0*t_dyn/0.01).ceil() as usize;
                let hold=adaptive_undamped_hold(&mut particles,physics,
                    cfg,&mut state,0.01,steps).unwrap();
                hold_s=hold.elapsed_s;
                hold_radius=100.0*hold.radius_relative_change;
                hold_energy=hold.energy_relative_drift;
                hold_momentum=hold.momentum_relative_drift;
                println!("    undamped 3t_dyn={:.3}s, ΔR={:+.3}%, ΔE/E={:.3e}, momentum drift={:.3e}",
                    hold_s,hold_radius,hold_energy,hold_momentum);
                assert!(hold.mass_relative_drift<1e-12);
                assert!(hold.momentum_relative_drift<1e-8);
            }
            assert!(result.report.mass_relative_error<1e-12);
            assert!(result.report.center_shift_m/star.radius_m<1e-8);
            assert!(result.report.final_accel_rms_m_s2.is_finite());
            assert!(result.report.final_density_log_rms.is_finite());
            csv.push_str(&format!(
                "{n},{seed_label},{gravity_label},{:.4},{:.8},{:.8},{:.7},{:.5},{:.8},{:.7},{:.5},{:.5},{:.7},{:.8},{:.8}\n",
                options.dt_s*options.iterations as f64,
                result.report.force_reduction(),result.report.final_density_log_rms,
                result.radius_bias*100.0,result.report.max_radius_m/1e3,
                result.minimum_omega,result.capped_fraction,elapsed,
                hold_s,hold_radius,hold_energy,hold_momentum));
        }
    }
    if let Some(path)=csv_path {
        std::fs::write(path,csv).expect("write experiment CSV");
        println!("wrote {path}");
    }
}
