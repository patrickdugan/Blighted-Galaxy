//! BG-PHY-17: matched full-duration physical hold for selected eta target.
//!
//! Retains the 0.65 solar mass mu_e=2 cold EOS reference, original
//! damping tau=0.7s, 3-second preparation and 3*t_dyn undamped hold.
//! In contrast to the frozen neighbor sweep, this integrates the dynamics.
//! Existing physical acceptance thresholds are unchanged.
//!
//! cargo run --release --bin phy17_hold -- --n 1024 --seed 0
//!   --layout glass --kernel c4 --target-neighbors 150
//!   --relax-steps 300 --hold-steps 788 --dt 0.01
//!   --csv bg-phy17-hold.csv --trace-csv bg-phy17-hold-trace.csv
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_neighbor_calibration;
use blighted_galaxy::sph_force_audit;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::stellar_relax::{relax_isolated_adaptive,RelaxParameters};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
use std::time::Instant;

fn value<'a>(args:&'a [String],key:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==key).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T{
    value(args,key).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn empty_center()->Massive3 {
    Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    }
}
fn append_row(text:&mut String,fields:&[String]){
    let _=writeln!(text,"{}",fields.join(","));
}
fn f(x:f64)->String{format!("{x:.9e}")}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",1024);
    let seed:u64=opt(&args,"--seed",0);
    let layout=value(&args,"--layout").unwrap_or("glass");
    let kernel=value(&args,"--kernel").unwrap_or("c4");
    let target:f64=opt(&args,"--target-neighbors",150.0);
    let fixed_eta:f64=opt(&args,"--fixed-eta",1.2);
    let relax_steps:usize=opt(&args,"--relax-steps",300);
    let hold_steps:usize=opt(&args,"--hold-steps",788);
    let dt:f64=opt(&args,"--dt",0.01);
    let calibration_steps:usize=opt(&args,"--calibration-steps",12);
    let glass_sweeps:usize=opt(&args,"--glass-sweeps",20);
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(layout=="glass"||layout=="shell");
    assert!(kernel=="m4"||kernel=="c4");
    assert!(target==0.0||(target>=2.0&&target<=n as f64));
    assert!(fixed_eta>0.5&&fixed_eta<3.0);
    assert!(relax_steps>=1&&relax_steps<=30000);
    assert!(hold_steps>=1&&hold_steps<=30000);
    assert!(dt.is_finite()&&dt>0.0&&dt<=0.1);
    assert!(glass_sweeps>=1&&glass_sweeps<=100);
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let original=seed_cubic_shells_oriented(
        &star,center,V3::default(),n,seed).0;
    let mut p=if layout=="glass" {
        sph_glass::relax_antipodal_angles(
            &original,&star,center,1.2,glass_sweeps).0
    } else {original};
    let sampled_radius=star.radius_at_mass_fraction(1.0-1.0/n as f64);
    let initial_profile=sph_massprofile::profile_diagnostics(&p,&star);
    assert!((initial_profile.radius_max_m-sampled_radius).abs()/star.radius_m<1e-10);
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity();
    if kernel=="c4" {cfg=cfg.with_wendland_c4();}
    let calibration=if target>0.0 {
        let result=sph_neighbor_calibration::calibrate(
            &p,star.radius_m,cfg,target,0.75,2.85,calibration_steps);
        assert!(result.bracketed&&result.abs_error<=target*0.05,
            "selected target cannot be calibrated within 5% realized neighbors");
        cfg.eta=result.eta;
        Some(result)
    } else {
        cfg.eta=fixed_eta;
        None
    };
    let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
    let initial_forces=sph_force_audit::frozen_force_audit(
        &p,&mut state,physics,cfg,&star,8);
    let initial_neighbors=sph_neighbor_calibration::measure(&p,&state,cfg).mean;
    println!("BG-PHY-17 hold: N={n} seed={seed} {layout} {kernel} target={target:.1} eta={:.6} realized={initial_neighbors:.2}, reference R={:.2} km, sampled R={:.2} km",
        cfg.eta,star.radius_m/1000.0,sampled_radius/1000.0);
    let start=Instant::now();
    let relaxed=relax_isolated_adaptive(
        &mut p,physics,&star,cfg,&mut state,
        RelaxParameters {
            iterations:relax_steps,dt_s:dt,damping_time_s:0.7,
        }).expect("finite stellar damping");
    assert!(relaxed.report.mass_relative_error<1e-12);
    let post=sph_massprofile::profile_diagnostics(&p,&star);
    let post_forces=sph_force_audit::frozen_force_audit(
        &p,&mut state,physics,cfg,&star,8);
    let post_neighbors=sph_neighbor_calibration::measure(&p,&state,cfg).mean;
    let post_exact={
        let exact_cfg=cfg.with_newton_density();
        let mut exact_state=state.clone();
        sph_force_audit::frozen_force_audit(
            &p,&mut exact_state,physics,exact_cfg,&star,8)
    };
    let mut central=empty_center();
    let mass0=sph3::total_mass_kg(&central,&p);
    let momentum0=sph3::total_momentum(&central,&p);
    let energy0=sph_adaptive::adaptive_energy_j(
        &central,&p,&mut state,physics,cfg);
    let mut worst_rmax=0.0_f64;
    let mut worst_r99=0.0_f64;
    let mut worst_cdf=post.radial_quantile_rms;
    let mut max_energy_drift=0.0_f64;
    let mut max_momentum_drift=0.0_f64;
    let mut trace=String::from(
        "N,seed,layout,kernel,target_neighbors,eta,t_s,rmax_km,r99_km,cdf_rms,rmax_drift_pct,r99_drift_pct,energy_drift,momentum_drift\n");
    let stride=(hold_steps/30).max(1);
    for step in 1..=hold_steps {
        sph_adaptive::advance(&mut central,&mut p,&mut state,physics,cfg,dt);
        if step%stride==0||step==hold_steps {
            let q=sph_massprofile::profile_diagnostics(&p,&star);
            let dr=(q.radius_max_m/post.radius_max_m-1.0).abs();
            let d99=(q.r99_m/post.r99_m-1.0).abs();
            worst_rmax=worst_rmax.max(dr);
            worst_r99=worst_r99.max(d99);
            worst_cdf=worst_cdf.max(q.radial_quantile_rms);
            let energy=sph_adaptive::adaptive_energy_j(
                &central,&p,&mut state,physics,cfg);
            let de=(energy-energy0).abs()/energy0.abs().max(1.0);
            let dm=(sph3::total_momentum(&central,&p)-momentum0).len()/
                (mass0*3.0e5);
            max_energy_drift=max_energy_drift.max(de);
            max_momentum_drift=max_momentum_drift.max(dm);
            assert!(de.is_finite()&&de<1e-3);
            assert!(dm.is_finite()&&dm<1e-8);
            assert!((sph3::total_mass_kg(&central,&p)-mass0).abs()/mass0<1e-12);
            append_row(&mut trace,&[
                n.to_string(),seed.to_string(),layout.to_string(),kernel.to_string(),
                f(target),f(cfg.eta),f(step as f64*dt),
                f(q.radius_max_m/1000.0),f(q.r99_m/1000.0),
                f(q.radial_quantile_rms),f(100.0*dr),f(100.0*d99),f(de),f(dm),
            ]);
        }
    }
    let elapsed=start.elapsed().as_secs_f64();
    let hold_duration=hold_steps as f64*dt;
    let dyn_time=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let full_hold=hold_duration>=3.0*dyn_time;
    let bias=100.0*(post.radius_max_m/star.radius_m-1.0);
    let physical_gate=full_hold&&bias.abs()<5.0&&
        worst_rmax<0.05&&worst_r99<0.05&&worst_cdf<0.05;
    let mut csv=String::from(
        "N,seed,layout,kernel,target_neighbors,eta,calibration_error,calibration_bracketed,initial_neighbors,post_neighbors,initial_sampling_bias_pct,post_rmax_km,post_r99_km,post_rmax_bias_pct,contraction_from_seed_pct,post_cdf_rms,initial_pressure_ref_rms,post_pressure_ref_rms,post_gravity_ref_rms,post_net_ref_rms,post_exact_root_pressure_rms,hold_s,t_dyn_s,full_hold,worst_rmax_drift_pct,worst_r99_drift_pct,max_cdf_rms,max_energy_drift,max_momentum_drift,physical_gate_pass,wall_s\n");
    append_row(&mut csv,&[
        n.to_string(),seed.to_string(),layout.to_string(),kernel.to_string(),
        f(target),f(cfg.eta),
        f(calibration.map(|x|x.abs_error).unwrap_or(0.0)),
        calibration.map(|x|x.bracketed).unwrap_or(true).to_string(),
        f(initial_neighbors),f(post_neighbors),
        f(100.0*(sampled_radius/star.radius_m-1.0)),
        f(post.radius_max_m/1000.0),f(post.r99_m/1000.0),
        f(bias),f(100.0*(post.radius_max_m/sampled_radius-1.0)),
        f(post.radial_quantile_rms),
        f(initial_forces.global.pressure_ref_rms),
        f(post_forces.global.pressure_ref_rms),
        f(post_forces.global.gravity_ref_rms),
        f(post_forces.global.net_ref_rms),
        f(post_exact.global.pressure_ref_rms),
        f(hold_duration),f(dyn_time),full_hold.to_string(),
        f(100.0*worst_rmax),f(100.0*worst_r99),
        f(worst_cdf),f(max_energy_drift),f(max_momentum_drift),
        physical_gate.to_string(),f(elapsed),
    ]);
    println!("post Rmax bias={bias:+.3}%, contraction={:+.3}%, P/G/net={:.4}/{:.4}/{:.4}, Rmax/R99 drift {:.2}%/{:.2}%, CDFmax={worst_cdf:.5}, gate={physical_gate}",
        100.0*(post.radius_max_m/sampled_radius-1.0),
        post_forces.global.pressure_ref_rms,
        post_forces.global.gravity_ref_rms,
        post_forces.global.net_ref_rms,
        100.0*worst_rmax,100.0*worst_r99);
    if let Some(path)=value(&args,"--csv") {
        std::fs::write(path,csv).expect("write physical hold summary");
        println!("wrote {path}");
    }
    if let Some(path)=value(&args,"--trace-csv") {
        std::fs::write(path,trace).expect("write physical hold trace");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn zero_mass_central_placeholder_is_finite() {
        let x=empty_center();
        assert_eq!(x.mass_kg,0.0);
        assert_eq!(x.vel_m_s,V3::default());
    }
}
