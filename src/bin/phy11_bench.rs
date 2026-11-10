//! BG-PHY-11 high-resolution isolated cold white-dwarf validation.
 //!
 //! cargo run --release --bin phy11_bench -- --n 1024 \
 //!   --relax-steps 300 --hold-steps 788 --csv bg-phy11-highres.csv
 //!
 //! Scientific acceptance is reported, never faked by weakening asserts.
 //! A 3-dynamical-time undamped hold must keep the star's radius stable.
 //! This remains an O(N²) exact symmetric gravity backend: cell neighbors
 //! accelerate only compact-support density and zeta terms.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_sampling::{self,StarSampling};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::sph_neighbors::CellNeighbors;
use blighted_galaxy::stellar_relax::{
    relax_isolated_adaptive,RelaxParameters,
};
use blighted_galaxy::white_dwarf::{
    DegenerateElectrons,WhiteDwarfProfile,profile_for_mass,
};
use std::time::Instant;
use std::fmt::Write as _;

fn value<'a>(args:&'a [String],key:&str)->Option<&'a str> {
    args.windows(2).find(|w|w[0]==key).map(|w|w[1].as_str())
}
fn parse<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T {
    value(args,key).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn mass_profile_row(
    csv:&mut String,phase:&str,sim_s:f64,profile:MassProfileDiagnostic,
    reference:&WhiteDwarfProfile,energy_err:f64,momentum_err:f64,wall_s:f64,
) {
    let radius_bias=(profile.radius_max_m/reference.radius_m-1.0)*100.0;
    let ref99=reference.radius_at_mass_fraction(0.99);
    let r99_bias=(profile.r99_m/ref99-1.0)*100.0;
    let _=writeln!(csv,
        "{phase},{sim_s:.5},{:.7},{:.6},{:.6},{:.6},{:.6},{:.6},{radius_bias:.6},{r99_bias:.6},NaN,{energy_err:.9e},{momentum_err:.9e},{wall_s:.4}",
        profile.radius_max_m/1e3,profile.r99_m/1e3,
        profile.r95_m/1e3,profile.r90_m/1e3,profile.r50_m/1e3,
        profile.radial_quantile_rms);
}
fn density_neighbor_audit(particles:&[Fluid3],state:&AdaptiveState)->(f64,f64,f64) {
    let max_h=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let neighbors=CellNeighbors::new(particles,max_h);
    let lists=neighbors.all_candidates();
    let mean=lists.iter().map(|l|l.len()).sum::<usize>() as f64/
        particles.len() as f64;
    let naive_start=Instant::now();
    let brute=(0..particles.len()).map(|i|{
        particles.iter().map(|p|
            p.mass_kg*sph3::cubic_spline_3d(
                (particles[i].pos_m-p.pos_m).len(),state.smoothing_m[i])
        ).sum::<f64>()
    }).collect::<Vec<_>>();
    let brute_sec=naive_start.elapsed().as_secs_f64();
    let indexed_start=Instant::now();
    let candidate=(0..particles.len()).map(|i|{
        lists[i].iter().map(|&j|
            particles[j].mass_kg*sph3::cubic_spline_3d(
                (particles[i].pos_m-particles[j].pos_m).len(),
                state.smoothing_m[i])
        ).sum::<f64>()
    }).collect::<Vec<_>>();
    let indexed_sec=indexed_start.elapsed().as_secs_f64();
    let max_relative=brute.iter().zip(candidate.iter()).map(|(a,b)|
        (a-b).abs()/a.abs().max(1e-50)).fold(0.0_f64,f64::max);
    assert!(max_relative<1e-12,"neighbor density audit {max_relative}");
    println!("  density neighbor audit: n={} mean candidates {:.1} ({:.1}% of all), max relative error {:.2e}",
        particles.len(),mean,100.0*mean/particles.len() as f64,max_relative);
    println!("  single-pass reference O(N^2) {:.5}s; cell-list {:.5}s",
        brute_sec,indexed_sec);
    (mean,max_relative,indexed_sec/brute_sec.max(1e-12))
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=parse(&args,"--n",128);
    let relax_steps:usize=parse(&args,"--relax-steps",150);
    let hold_steps:usize=parse(&args,"--hold-steps",262);
    let dt:f64=parse(&args,"--dt",0.01);
    let output=value(&args,"--csv");
    assert!(n>=32&&n<=8192&&n%2==0);
    assert!(relax_steps>=1&&relax_steps<=100_000);
    assert!(hold_steps>=1&&hold_steps<=100_000);
    assert!(dt>0.0&&dt<=0.1&&dt.is_finite());
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let cfg=AdaptiveConfig::from_star_radius(star.radius_m).with_cubic_gravity();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    println!("BG-PHY-11 / cold white-dwarf equilibrium and neighbor audit");
    println!("  M=0.65 solar, hydrostatic R={:.5} km, t_dyn={t_dyn:.5}s",
        star.radius_m/1e3);
    println!("  N={n} relax steps={relax_steps} ({:.3}s), undamped steps={hold_steps} ({:.3}s, {:.3} dyn times)",
        relax_steps as f64*dt,hold_steps as f64*dt,hold_steps as f64*dt/t_dyn);
    let mut particles=sph_sampling::seeded_star(
        &star,V3::new(1.0e8,0.0,0.0),V3::default(),n,
        StarSampling::ApproximatelyCubicShells);
    let mut state=AdaptiveState::new(&particles,star.radius_m,cfg);
    let seed=sph_massprofile::profile_diagnostics(&particles,&star);
    println!("  initial mass-CDF rms {:.6}, R99 {:.2} km, Rmax {:.2} km",
        seed.radial_quantile_rms,seed.r99_m/1e3,seed.radius_max_m/1e3);
    let (mean_neighbors,density_error,neighbor_time_ratio)=
        density_neighbor_audit(&particles,&state);
    let timer=Instant::now();
    let report=relax_isolated_adaptive(
        &mut particles,physics,&star,cfg,&mut state,
        RelaxParameters {
            iterations:relax_steps,dt_s:dt,damping_time_s:0.7,
        }).unwrap();
    let relax_wall=timer.elapsed().as_secs_f64();
    let relaxed=sph_massprofile::profile_diagnostics(&particles,&star);
    println!("  damped: Rmax {:.2} km ({:+.3}%); R99 {:.2} km; CDF rms {:.5}; force ratio {:.5}; elapsed {:.2}s",
        relaxed.radius_max_m/1e3,
        100.0*(relaxed.radius_max_m/star.radius_m-1.0),
        relaxed.r99_m/1e3,relaxed.radial_quantile_rms,
        report.report.force_reduction(),relax_wall);
    let mut csv=String::from(
        "phase,sim_s,rmax_km,r99_km,r95_km,r90_km,r50_km,radial_cdf_rms,radius_bias_percent,r99_bias_percent,rho_log_rms_placeholder,energy_drift,momentum_drift,wall_s\n"
    );
    let mut central=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    let e0=sph_adaptive::adaptive_energy_j(
        &central,&particles,&mut state,physics,cfg);
    let p0=sph3::total_momentum(&central,&particles);
    let m0=sph3::total_mass_kg(&central,&particles);
    mass_profile_row(&mut csv,"post_damping",0.0,relaxed,&star,
        0.0,0.0,relax_wall);
    let mut max_rmax_drift=0.0_f64;
    let mut max_r99_drift=0.0_f64;
    let mut max_radial_cdf=relaxed.radial_quantile_rms;
    let hold_time=Instant::now();
    let interval=(hold_steps/40).max(1);
    for step in 1..=hold_steps {
        sph_adaptive::advance(&mut central,&mut particles,
            &mut state,physics,cfg,dt);
        if step%interval==0 || step==hold_steps {
            let profile=sph_massprofile::profile_diagnostics(&particles,&star);
            let rmax_drift=profile.radius_max_m/relaxed.radius_max_m-1.0;
            let r99_drift=profile.r99_m/relaxed.r99_m-1.0;
            max_rmax_drift=max_rmax_drift.max(rmax_drift.abs());
            max_r99_drift=max_r99_drift.max(r99_drift.abs());
            max_radial_cdf=max_radial_cdf.max(profile.radial_quantile_rms);
            {
                let energy=sph_adaptive::adaptive_energy_j(
                    &central,&particles,&mut state,physics,cfg);
                let momentum=sph3::total_momentum(&central,&particles);
                let mass=sph3::total_mass_kg(&central,&particles);
                let e_drift=(energy-e0).abs()/e0.abs().max(1.0);
                let p_drift=(momentum-p0).len()/(m0*3e5);
                assert!((mass-m0).abs()/m0<1e-12);
                assert!(p_drift<1e-8);
                assert!(e_drift.is_finite()&&e_drift<1e-3,
                    "energy drift {e_drift}");
                let elapsed=step as f64*dt;
                println!("  undamped t={elapsed:.3}s ({:.2} t_dyn): Rmax {:.2} km ({:+.2}%), R99 {:.2} km ({:+.2}%), CDF rms {:.5}, dE {:.3e}",
                    elapsed/t_dyn,profile.radius_max_m/1e3,
                    rmax_drift*100.0,profile.r99_m/1e3,
                    r99_drift*100.0,profile.radial_quantile_rms,e_drift);
                mass_profile_row(&mut csv,"undamped",elapsed,
                    profile,&star,e_drift,p_drift,hold_time.elapsed().as_secs_f64());
            }
        }
    }
    let last=sph_massprofile::profile_diagnostics(&particles,&star);
    let post_relax_radius_error=(relaxed.radius_max_m/star.radius_m-1.0).abs();
    let stable=max_rmax_drift<0.05 && max_r99_drift<0.05;
    let profile_ok=max_radial_cdf<0.05;
    let passed=post_relax_radius_error<0.05&&stable&&profile_ok;
    println!("  worst hold variations: outer R {:.2}%, R99 {:.2}%; max CDF rms {:.5}",
        max_rmax_drift*100.0,max_r99_drift*100.0,max_radial_cdf);
    println!("  5%-radius / 5%-stability / 5%-CDF physics gates: radius={}, stability={}, CDF={}, overall={}",
        post_relax_radius_error<0.05,stable,profile_ok,passed);
    println!("  final undamped Rmax {:.2} km, R99 {:.2} km",
        last.radius_max_m/1000.0,last.r99_m/1000.0);
    let _=writeln!(csv,
        "#meta,N={n},t_dyn={t_dyn:.6},mean_neighbors={mean_neighbors:.3},density_max_relative_error={density_error:.3e},neighbor_time_ratio={neighbor_time_ratio:.4},max_outer_drift={max_rmax_drift:.6},max_r99_drift={max_r99_drift:.6},cdf_max={max_radial_cdf:.6},physics_gate_pass={passed}");
    if let Some(path)=output {
        std::fs::write(path,csv).expect("write BG-PHY-11 diagnostics");
        println!("  wrote {path}");
    }
}
