//! BG-PHY-08: moving-observer Kerr light cones and variable-h cold SPH.
//! Run: cargo run --release --bin phy08_bench -- --csv phy08.csv
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::KerrCamera;
use blighted_galaxy::kerr_moving_observer::{
    KerrMovingObserver,ObserverMotion,observe_recorded_source,
    source_aligned_to_observer_ray,
};
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::stellar_relax::{
    adaptive_undamped_hold,relax_isolated_adaptive,
    relax_isolated_white_dwarf,RelaxParameters,
};
use blighted_galaxy::white_dwarf::{
    DegenerateElectrons,profile_for_mass,
};
use std::f64::consts::FRAC_PI_2;
use std::time::Instant;

fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let csv_path=args.windows(2).find(|w|w[0]=="--csv").map(|w|w[1].as_str());
    let mass=10.0*SOLAR_MASS_KG;
    println!("BG-PHY-08 MOVING OBSERVER / CURVED NULL PHOTONS");
    for spin in [0.0,0.75,-0.75] {
        let geo=KerrGeometry::new(mass,spin).unwrap();
        let cam=KerrCamera::new(geo,85.0,FRAC_PI_2,36.0).unwrap();
        let observer=KerrMovingObserver::new(cam,
            ObserverMotion::CircularGeodesic{direction:1},0.25).unwrap();
        let now=0.06;
        let source=source_aligned_to_observer_ray(
            &observer,0.0,20.0,now,1,1.0).unwrap();
        let history=source.record(&cam,now-0.015,now+0.005,0.00008).unwrap();
        let acq=observe_recorded_source(&observer,&history,
            now,1.0,4.0,65,17.0).expect("moving Kerr optical match");
        let gap=(acq.apparent_pos_m-acq.current_pos_m).len()/geo.r_g_m;
        println!("  spin {spin:+.2} : beta_ZAMO={:+.6}, dτ/dt={:.6}, delay={:.9}s, target shift={:.3} rg, alpha={:+.4} rg",
            observer.beta_zamo,observer.dproper_dcoordinate(),
            acq.delay_s,gap,acq.alpha_camera_rg);
        assert!(acq.delay_s>0.0&&acq.delay_s<0.015);
        assert!(gap>1.0);
        assert!(acq.ray.max_radial_residual<1e-4);
        assert!(acq.ray.max_polar_residual<1e-4);
    }

    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let cfg=AdaptiveConfig::from_star_radius(star.radius_m);
    let mut csv=String::from("N,method,damped_seconds,force_ratio,log_rho_before,log_rho_after,radius_m,reference_radius_m,relative_radius_bias,omega_min,h_min_m,h_max_m,h_cap_fraction,cpu_seconds,hold_seconds,hold_radius_change,hold_energy_drift,hold_momentum_drift\n");

    println!("BG-PHY-08 VARIABLE-h COLD WHITE DWARF RELAXATION");
    println!("  hydrostatic reference: M=0.65 Msun R={:.2} km t_dyn={:.3} s",
        star.radius_m/1e3,(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt());
    for n in [64_usize,128,256] {
        let center=V3::new(1.0e8,0.0,0.0);
        let seed=sph3::seed_from_hydrostatic(&star,center,V3::default(),n);
        let options=RelaxParameters{
            iterations:if n>=256 {250} else {500},
            dt_s:0.01,damping_time_s:0.7,
        };
        let seconds=options.dt_s*options.iterations as f64;
        let mut fixed=seed.clone();
        let start=Instant::now();
        let old=relax_isolated_white_dwarf(
            &mut fixed,physics,&star,options).unwrap();
        let fixed_elapsed=start.elapsed().as_secs_f64();
        let fixed_bias=old.max_radius_m/star.radius_m-1.0;
        println!("  N={n} fixed h, {seconds:.1}s: radius={:.2} km ({:+.1}%), force ratio {:.5}, log density {:.4}",
            old.max_radius_m/1e3,fixed_bias*100.0,old.force_reduction(),
            old.final_density_log_rms);
        csv.push_str(&format!(
            "{n},fixed,{seconds:.4},{:.8},{:.8},{:.8},{:.5},{:.5},{:.8},NaN,NaN,NaN,NaN,{:.5},NaN,NaN,NaN,NaN\n",
            old.force_reduction(),old.initial_density_log_rms,
            old.final_density_log_rms,old.max_radius_m,star.radius_m,
            fixed_bias,fixed_elapsed));

        let mut p=seed;
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let start=Instant::now();
        let adaptive=relax_isolated_adaptive(
            &mut p,physics,&star,cfg,&mut state,options).unwrap();
        let elapsed=start.elapsed().as_secs_f64();
        println!("  N={n} adaptive, {seconds:.1}s: radius={:.2} km ({:+.1}%), force ratio {:.5}, log density {:.4}, Ωmin {:.4}, h=[{:.2},{:.2}] Mm, capped {:.1}%",
            adaptive.report.max_radius_m/1e3,
            adaptive.radius_bias*100.0,
            adaptive.report.force_reduction(),
            adaptive.report.final_density_log_rms,
            adaptive.minimum_omega,
            adaptive.minimum_h_m/1e6,adaptive.maximum_h_m/1e6,
            adaptive.capped_fraction*100.0);
        assert!(adaptive.report.mass_relative_error<1e-12);
        assert!(adaptive.report.center_shift_m/star.radius_m<1e-8);
        assert!(adaptive.minimum_omega>0.0);
        let mut hold_seconds=f64::NAN;
        let mut hold_radius=f64::NAN;
        let mut hold_energy=f64::NAN;
        let mut hold_p=f64::NAN;
        if n==64 {
            let dyn_time=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
            let dt=0.01;
            let steps=(3.0*dyn_time/dt).ceil() as usize;
            let hold=adaptive_undamped_hold(&mut p,physics,cfg,
                &mut state,dt,steps).unwrap();
            hold_seconds=hold.elapsed_s;
            hold_radius=hold.radius_relative_change;
            hold_energy=hold.energy_relative_drift;
            hold_p=hold.momentum_relative_drift;
            println!("    N=64 undamped hold {:.2}s (3 t_dyn): ΔR/R={:+.3}%, ΔE/E={:.3e}, ΔP/scale={:.3e}",
                hold_seconds,100.0*hold_radius,hold_energy,hold_p);
            assert!(hold.mass_relative_drift<1e-12);
            assert!(hold.momentum_relative_drift<1e-9);
        }
        csv.push_str(&format!(
            "{n},adaptive,{seconds:.4},{:.8},{:.8},{:.8},{:.5},{:.5},{:.8},{:.8},{:.5},{:.5},{:.8},{:.5},{:.5},{:.8},{:.8},{:.8}\n",
            adaptive.report.force_reduction(),
            adaptive.report.initial_density_log_rms,
            adaptive.report.final_density_log_rms,
            adaptive.report.max_radius_m,star.radius_m,
            adaptive.radius_bias,adaptive.minimum_omega,
            adaptive.minimum_h_m,adaptive.maximum_h_m,
            adaptive.capped_fraction,elapsed,
            hold_seconds,hold_radius,hold_energy,hold_p));
    }
    if let Some(path)=csv_path {
        std::fs::write(path,csv).expect("write adaptive SPH comparison CSV");
        println!("wrote {path}");
    }
}
