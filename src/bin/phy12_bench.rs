//! BG-PHY-12: matched-duration white-dwarf force balance and damping tests.
 //!
 //! cargo run --release --bin phy12_bench -- --n 256 --relax-steps 600 \
 //!   --hold-steps 788 --csv bg-phy12-timeseries.csv \
 //!   --shell-csv bg-phy12-shells.csv
 //!
 //! A 1D cold-electron hydrostatic star is the independent reference.
 //! The original constant-tau schedule and geometric taper are compared
 //! with identical particle seeds, force laws, dt and physical durations.
 //! Physical acceptance is reported separately from numerical invariant
 //! assertions. Neither schedule is assumed to produce equilibrium.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::{self,StarSampling};
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::stellar_balance::{
    self,StellarBalance,
};
use blighted_galaxy::stellar_relax::{
    relax_isolated_adaptive_with_schedule,RelaxParameters,DampingSchedule,
};
use blighted_galaxy::white_dwarf::{
    DegenerateElectrons,WhiteDwarfProfile,profile_for_mass,
};
use std::fmt::Write as _;
use std::time::Instant;

fn arg_value<'a>(args:&'a [String],name:&str)->Option<&'a str> {
    args.windows(2).find(|x|x[0]==name).map(|x|x[1].as_str())
}
fn arg<T:std::str::FromStr>(args:&[String],name:&str,default:T)->T {
    arg_value(args,name).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn make_row(
    data:&mut String,label:&str,phase:&str,time:f64,
    profile:&MassProfileDiagnostic,
    balance:&StellarBalance,
    hydro:&WhiteDwarfProfile,
    energy_drift:f64,momentum_drift:f64,
) {
    let _=writeln!(data,
        "{label},{phase},{time:.6},{:.5},{:.5},{:.6},{:.8},{:.8},{:.8},{:.8},{:.6},{energy_drift:.9e},{momentum_drift:.9e}",
        profile.radius_max_m/1e3,profile.r99_m/1e3,
        100.0*(profile.radius_max_m/hydro.radius_m-1.0),
        profile.radial_quantile_rms,
        balance.pressure_to_gravity_virial,
        balance.residual_to_gravity_rms,
        balance.maximum_radial_shell_residual,
        balance.net_force_rms_m_s2,
    );
}
fn make_shell_rows(
    data:&mut String,label:&str,phase:&str,
    bal:&StellarBalance,
) {
    for (i,s) in bal.shells.iter().enumerate() {
        let _=writeln!(data,
            "{label},{phase},{i},{:.6},{:.6},{},{:.5},{:.7e},{:.7e},{:.7e},{:.7e},{:.7e},{:.7e}",
            s.mass_fraction_lo,s.mass_fraction_hi,s.particle_count,
            s.mean_radius_m/1e3,s.pressure_out_m_s2,
            s.gravity_in_m_s2,s.hydrostatic_gravity_in_m_s2,
            s.net_radial_m_s2,s.net_radial_rms_m_s2,
            s.tangential_net_rms_m_s2);
    }
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=arg(&args,"--n",128);
    let relax_steps:usize=arg(&args,"--relax-steps",300);
    let hold_steps:usize=arg(&args,"--hold-steps",262);
    let dt:f64=arg(&args,"--dt",0.01);
    let damping_tau:f64=arg(&args,"--tau",0.7);
    let output=arg_value(&args,"--csv");
    let shell_output=arg_value(&args,"--shell-csv");
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!((1..=50_000).contains(&relax_steps));
    assert!((1..=50_000).contains(&hold_steps));
    assert!(dt.is_finite()&&dt>0.0&&dt<=0.1);
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let seed=sph_sampling::seeded_star(
        &star,V3::new(1.0e8,0.0,0.0),V3::default(),n,
        StarSampling::ApproximatelyCubicShells);
    let options=RelaxParameters{
        iterations:relax_steps,dt_s:dt,damping_time_s:damping_tau,
    };
    let schedules=[
        ("constant",DampingSchedule::Constant{tau_s:damping_tau}),
        ("geometric_taper",DampingSchedule::TaperToFree{
            tau_initial_s:0.15,tau_final_s:1.5,
            taper_start:0.40,taper_end:0.90,
        }),
    ];
    let mut history=String::from(
        "schedule,phase,time_s,rmax_km,r99_km,rmax_bias_pct,radial_cdf_rms,pressure_gravity_virial,net_gravity_rms,max_shell_residual,force_rms_m_s2,energy_drift,momentum_drift\n"
    );
    let mut shells=String::from(
        "schedule,phase,shell,mass_fraction_lo,mass_fraction_hi,n,mean_radius_km,pressure_out_m_s2,gravity_in_m_s2,hydrostatic_in_m_s2,net_radial_m_s2,net_radial_rms_m_s2,tangential_rms_m_s2\n"
    );
    println!("BG-PHY-12 radial force balance / damping comparison");
    println!("  N={n}, R_ref={:.3} km, t_dyn={t_dyn:.6}s, damp={:.3}s, hold={:.3}s ({:.2} t_dyn)",
        star.radius_m/1e3,relax_steps as f64*dt,
        hold_steps as f64*dt,hold_steps as f64*dt/t_dyn);
    let mut global_results=Vec::new();
    for (label,schedule) in schedules {
        let mut particles=seed.clone();
        let mut state=AdaptiveState::new(
            &particles,star.radius_m,cfg);
        let before_profile=sph_massprofile::profile_diagnostics(
            &particles,&star);
        let before_balance=stellar_balance::radial_balance(
            &particles,&mut state,physics,cfg,&star,8);
        make_row(&mut history,label,"initial",0.0,
            &before_profile,&before_balance,&star,0.0,0.0);
        make_shell_rows(&mut shells,label,"initial",&before_balance);
        let started=Instant::now();
        let relaxed=relax_isolated_adaptive_with_schedule(
            &mut particles,physics,&star,cfg,&mut state,
            options,schedule).expect("finite SPH relaxation");
        let wall=started.elapsed().as_secs_f64();
        let after_profile=sph_massprofile::profile_diagnostics(
            &particles,&star);
        let after_balance=stellar_balance::radial_balance(
            &particles,&mut state,physics,cfg,&star,8);
        println!("  {label:18}: after damping Rmax={:.2} km ({:+.3}%), R99={:.2} km, CDF rms={:.5}, force ratio={:.4}, virial(P/|G|)={:.5}, shell imbalance={:.3}, wall={wall:.2}s",
            after_profile.radius_max_m/1e3,
            100.0*(after_profile.radius_max_m/star.radius_m-1.0),
            after_profile.r99_m/1e3,
            after_profile.radial_quantile_rms,
            relaxed.report.force_reduction(),
            after_balance.pressure_to_gravity_virial,
            after_balance.maximum_radial_shell_residual);
        make_row(&mut history,label,"post_damp",0.0,
            &after_profile,&after_balance,&star,0.0,0.0);
        make_shell_rows(&mut shells,label,"post_damp",&after_balance);
        let mut central=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let mass=sph3::total_mass_kg(&central,&particles);
        let momentum0=sph3::total_momentum(&central,&particles);
        let energy0=sph_adaptive::adaptive_energy_j(
            &central,&particles,&mut state,physics,cfg);
        let mut max_rmax_drift=0.0_f64;
        let mut max_r99_drift=0.0_f64;
        let mut max_cdf=after_profile.radial_quantile_rms;
        let sample_interval=(hold_steps/20).max(1);
        for step in 1..=hold_steps {
            sph_adaptive::advance(&mut central,&mut particles,
                &mut state,physics,cfg,dt);
            if step%sample_interval==0||step==hold_steps {
                let observed=sph_massprofile::profile_diagnostics(
                    &particles,&star);
                max_rmax_drift=max_rmax_drift.max((
                    observed.radius_max_m/after_profile.radius_max_m-1.0).abs());
                max_r99_drift=max_r99_drift.max((
                    observed.r99_m/after_profile.r99_m-1.0).abs());
                max_cdf=max_cdf.max(observed.radial_quantile_rms);
                // Recompute force and local EOS density consistently at
                // every recorded time, while maintaining undamped dynamics.
                let balance=stellar_balance::radial_balance(
                    &particles,&mut state,physics,cfg,&star,8);
                let e=sph_adaptive::energy_j_cached(
                    &central,&particles,&state,physics,cfg);
                let pe=sph3::total_momentum(&central,&particles);
                let err=(e-energy0).abs()/energy0.abs().max(1.0);
                let dp=(pe-momentum0).len()/(mass*3e5);
                assert!(err.is_finite()&&err<1e-3,
                    "energy integration error {err}");
                assert!(dp<1e-8&&dp.is_finite(),
                    "momentum integration error {dp}");
                assert!((sph3::total_mass_kg(&central,&particles)-mass).abs()/
                    mass<1e-12);
                let time=step as f64*dt;
                make_row(&mut history,label,"undamped",time,
                    &observed,&balance,&star,err,dp);
                if step==hold_steps {
                    make_shell_rows(&mut shells,label,"final",&balance);
                    println!("    undamped {:6.3}s: Rmax change={:+.3}%, R99 change={:+.3}%, force balance={:.5}, virial={:.5}, CDF={:.5}, relative dE={:.3e}",
                        time,100.0*(observed.radius_max_m/
                            after_profile.radius_max_m-1.0),
                        100.0*(observed.r99_m/after_profile.r99_m-1.0),
                        balance.residual_to_gravity_rms,
                        balance.pressure_to_gravity_virial,
                        observed.radial_quantile_rms,err);
                }
            }
        }
        let radius_bias=(after_profile.radius_max_m/star.radius_m-1.0).abs();
        let stable=max_rmax_drift<0.05&&max_r99_drift<0.05;
        let profile_good=max_cdf<0.05;
        let passed=radius_bias<0.05&&stable&&profile_good;
        println!("    MAX excursion: Rmax {:.3}% R99 {:.3}% CDF rms {:.5}; gates radius={}, stable={}, CDF={}, overall={}",
            max_rmax_drift*100.0,max_r99_drift*100.0,max_cdf,
            radius_bias<0.05,stable,profile_good,passed);
        global_results.push((label,radius_bias,max_rmax_drift,max_r99_drift,
            max_cdf,passed));
    }
    println!("  COMMON-DURATION A/B:");
    for (label,radius,outer,r99,cdf,pass) in global_results {
        println!("    {label}: radius error={:.3}%, max Rmax drift={:.3}%, max R99 drift={:.3}%, CDF={cdf:.5}, pass={pass}",
            radius*100.0,outer*100.0,r99*100.0);
    }
    if let Some(path)=output {
        std::fs::write(path,history).expect("write force time series");
        println!("  wrote {path}");
    }
    if let Some(path)=shell_output {
        std::fs::write(path,shells).expect("write radial shells");
        println!("  wrote {path}");
    }
}
