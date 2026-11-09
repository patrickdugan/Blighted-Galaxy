//! Headless BG-PHY-07 Kerr optical worldline and white-dwarf relaxation bench.
 //! cargo run --release --bin optical_relax_bench
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::KerrCamera;
use blighted_galaxy::kerr_optics::{
    source_aligned_to_past_ray,observe_combatant,
};
use blighted_galaxy::physics::{C_SI,SOLAR_MASS_KG};
use blighted_galaxy::stellar_relax::{
    relax_isolated_white_dwarf,RelaxParameters,
};
use blighted_galaxy::white_dwarf::{
    profile_for_mass,DegenerateElectrons,
};
use blighted_galaxy::sph3::{self,Sph3Config,V3};
use std::f64::consts::FRAC_PI_2;
use std::time::Instant;

fn main() {
    let argv=std::env::args().skip(1).collect::<Vec<_>>();
    let flag=|name:&str|argv.windows(2)
        .find(|x|x[0]==name).map(|x|x[1].as_str());
    let convergence_csv=flag("--csv");
    let optical_csv=flag("--optical-csv");
    let mut csv=String::from(
        "particles,scheme,physical_s,force_initial_m_s2,force_final_m_s2,force_ratio,log_density_initial,log_density_final,max_radius_km,reference_radius_km,max_speed_km_s,run_seconds\n");
    let mass=10.0*SOLAR_MASS_KG;
    println!("BG-PHY-07: recorded curved Kerr light-cone acquisition");
    for spin in [0.0,0.8,-0.8] {
        let geo=KerrGeometry::new(mass,spin).unwrap();
        let camera=KerrCamera::new(geo,85.0,FRAC_PI_2,36.0).unwrap();
        let now=0.06;
        let orbit=source_aligned_to_past_ray(
            &camera,0.0,20.0,now,1,1.0).unwrap();
        let history=orbit.record(&camera,0.0,0.09,0.0001).unwrap();
        let observed=observe_combatant(
            &camera,&history,now,1.0,5.0,81,17.0)
            .expect("curved photon must intersect source worldline");
        let euclidean=(
            camera.geometry.r_g_m*85.0-observed.observed_pos_m.len())/C_SI;
        println!(
            "  Kerr spin {spin:+.2} : delay {:.9}s, radial Euclidean {:.9}s, alpha {:+.3} rg, observed/current separation {:.3} rg",
            observed.light_delay_s,euclidean,observed.alpha_rg,
            (observed.observed_pos_m-observed.present_pos_m).len()/geo.r_g_m
        );
        println!(
            "    radial constraint {:.3e}; polar {:.3e}, {} sampled photons, historical validity {}",
            observed.radial_error,observed.polar_error,
            observed.ray.samples.len(),history.validate_segments()
        );
        if spin==0.8 {
            if let Some(path)=optical_csv {
                let mut data=String::from(
                    "photon_t_offset_s,photon_x_m,photon_y_m,source_x_m,source_y_m\n");
                for v in &observed.ray.samples {
                    if let Some(sample)=history.sample_at(now+v.t_offset_s) {
                        data.push_str(&format!("{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                            v.t_offset_s,v.chart_xy_m.x,v.chart_xy_m.y,
                            sample.pos_m.x,sample.pos_m.y));
                    }
                }
                std::fs::write(path,data).expect("write Kerr optical path CSV");
            }
        }
        assert!(observed.light_delay_s>0.0);
        assert!(observed.radial_error<1e-4);
        assert!(observed.polar_error<1e-4);
        assert!(history.validate_segments());
        if spin==0.0 {assert!(observed.light_delay_s>euclidean);}
    }

    println!("BG-PHY-07: isolated SPH relaxation with cold degenerate EOS");
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let profile=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let cfg=Sph3Config::cold_white_dwarf(2.0).unwrap();
    for (count,options,label) in [
        (64_usize,RelaxParameters::baseline(),"baseline"),
        (128_usize,RelaxParameters::baseline(),"baseline"),
        (128_usize,RelaxParameters {
            iterations:1500,dt_s:0.01,damping_time_s:0.7,
        },"extended"),
        (256_usize,RelaxParameters {
            iterations:1500,dt_s:0.01,damping_time_s:0.7,
        },"extended"),
        (512_usize,RelaxParameters {
            iterations:1500,dt_s:0.01,damping_time_s:0.7,
        },"extended"),
    ] {
        let mut particles=sph3::seed_from_hydrostatic(
            &profile,V3::new(1.0e8,0.0,0.0),V3::default(),count);
        let t=Instant::now();
        let report=relax_isolated_white_dwarf(
            &mut particles,cfg,&profile,options)
            .expect("finite relaxed particles");
        println!(
            "  N={} {} (physical {:.1}s) runtime {:.3}s, RMS accel {:e}->{:e}, force ratio {:.5}, log-density mismatch {:.5}->{:.5}",
            count,label,options.iterations as f64*options.dt_s,t.elapsed().as_secs_f64(),
            report.initial_accel_rms_m_s2,
            report.final_accel_rms_m_s2,
            report.force_reduction(),
            report.initial_density_log_rms,
            report.final_density_log_rms
        );
        println!(
            "    max radius {:.3} km, relative CM drift {:.3e}, max residual speed {:.3} km/s",
            report.max_radius_m/1000.0,
            report.center_shift_m/profile.radius_m,
            report.final_max_speed_relative_m_s/1000.0
        );
        csv.push_str(&format!(
            "{},{},{:.5},{:.12e},{:.12e},{:.8},{:.8},{:.8},{:.4},{:.4},{:.4},{:.4}\n",
            count,label,options.iterations as f64*options.dt_s,
            report.initial_accel_rms_m_s2,report.final_accel_rms_m_s2,
            report.force_reduction(),report.initial_density_log_rms,
            report.final_density_log_rms,report.max_radius_m/1000.0,
            profile.radius_m/1000.0,
            report.final_max_speed_relative_m_s/1000.0,
            t.elapsed().as_secs_f64()
        ));
        assert!(report.mass_relative_error<1e-12);
        assert!(report.center_shift_m/profile.radius_m<1e-8);
        assert!(report.force_reduction().is_finite());
        assert!(report.max_radius_m<3.0*profile.radius_m);
    }
    if let Some(path)=convergence_csv {
        std::fs::write(path,csv).expect("write stellar convergence CSV");
    }
}
