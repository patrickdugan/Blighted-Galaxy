//! BG-PHY-15: frozen hydrostatic force decomposition, with seed-radius bias.
//!
//! cargo run --release --bin phy15_bench -- --n 256 --seeds 0,17
//!   --relax-steps 300 --dt 0.01 --csv phy15.csv --shell-csv phy15-shells.csv
//!
//! Each law starts from IDENTICAL hydrostatic positions at a given seed.
//! The initial state and the post-damping state are audited without moving
//! particles during force evaluation. The finite-sample last-particle radius
//! is compared with the exact reference inverse mass CDF at 1-1/N.
//! Snapshot force error metrics have their own reference-gravity RMS scale.
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{self,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::sph_force_audit::{self,FrozenForceAudit};
use blighted_galaxy::stellar_relax::{relax_isolated_adaptive,RelaxParameters};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;

fn value<'a>(args:&'a [String],flag:&str)->Option<&'a str> {
    args.windows(2).find(|w|w[0]==flag).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],flag:&str,default:T)->T {
    value(args,flag).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn seeds(args:&[String])->Vec<u64> {
    let list=value(args,"--seeds").unwrap_or("0,17")
        .split(',').map(|x|x.parse::<u64>().expect("invalid seed"))
        .collect::<Vec<_>>();
    assert!(!list.is_empty()&&list.len()<=8);
    assert!(list.windows(2).all(|w|w[0]<w[1]));
    list
}
fn growth_eta(n:usize)->f64 {
    1.05*(n as f64/256.0).powf(1.0/6.0)
}
fn append_snapshot(
    summary:&mut String,shell_csv:&mut String,
    n:usize,seed:u64,law:&str,eta:f64,stage:&str,seconds:f64,
    ref_radius:f64,sampling_bias_pct:f64,
    contraction_pct:f64,p:&MassProfileDiagnostic,a:&FrozenForceAudit,
) {
    let g=a.global;
    let rmax_bias=100.0*(p.radius_max_m/ref_radius-1.0);
    let _=writeln!(summary,
        "{n},{seed},{law},{eta:.7},{stage},{seconds:.4},{:.6},{rmax_bias:.7},{sampling_bias_pct:.7},{contraction_pct:.7},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8}",
        p.radius_max_m/1000.0,
        p.radial_quantile_rms,
        g.reference_gravity_rms_m_s2,
        g.pressure_ref_rms,g.gravity_ref_rms,g.net_ref_rms,
        g.pressure_radial_ratio,g.gravity_radial_ratio,
        g.pressure_tangent_rms,g.gravity_tangent_rms,g.net_tangent_rms,
        ref_radius/1000.0,
    );
    for (index,s) in a.shells.iter().enumerate() {
        let m=s.metrics;
        let _=writeln!(shell_csv,
            "{n},{seed},{law},{eta:.7},{stage},{index},{:.6},{:.6},{},{:.6},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8},{:.8}",
            s.mass_fraction_lo,s.mass_fraction_hi,m.particle_count,
            m.mean_radius_m/1000.0,
            m.reference_gravity_rms_m_s2,
            m.pressure_ref_rms,m.gravity_ref_rms,m.net_ref_rms,
            m.pressure_radial_ratio,m.gravity_radial_ratio,
            m.pressure_tangent_rms,m.gravity_tangent_rms,m.net_tangent_rms,
        );
    }
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",128);
    let relax_steps:usize=opt(&args,"--relax-steps",300);
    let dt:f64=opt(&args,"--dt",0.01);
    let shell_count:usize=opt(&args,"--shells",8);
    let fixed_eta:f64=opt(&args,"--fixed-eta",1.2);
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(relax_steps>0&&relax_steps<=30000);
    assert!(dt.is_finite()&&dt>0.0&&dt<=0.1);
    assert!(shell_count>=2&&shell_count<=n/4);
    assert!(fixed_eta>0.5&&fixed_eta<3.0);
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let expected_seed_rmax=star.radius_at_mass_fraction(1.0-1.0/n as f64);
    let expected_sampling_bias_pct=100.0*(expected_seed_rmax/star.radius_m-1.0);
    let mut summary=String::from(
        "N,seed,law,eta,stage,elapsed_s,rmax_km,rmax_bias_pct,initial_sampling_bias_pct,contraction_from_seed_pct,cdf_rms,reference_gravity_rms_m_s2,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms,net_tangent_rms,reference_radius_km\n"
    );
    let mut shell_csv=String::from(
        "N,seed,law,eta,stage,shell,mass_fraction_lo,mass_fraction_hi,particle_count,mean_radius_km,reference_gravity_rms_m_s2,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms,net_tangent_rms\n"
    );
    println!("BG-PHY-15 N={n}, Rref={:.3} km, sampled last-particle Rref={:.3} km (bias {expected_sampling_bias_pct:+.3}%)",
        star.radius_m/1000.0,expected_seed_rmax/1000.0);
    let eta_laws=[("fixed",fixed_eta),("growth",growth_eta(n))];
    for seed in seeds(&args) {
        let (initial,_)=seed_cubic_shells_oriented(
            &star,V3::new(1e8,0.0,0.0),V3::default(),n,seed);
        let seed_profile=sph_massprofile::profile_diagnostics(&initial,&star);
        assert!((seed_profile.radius_max_m-expected_seed_rmax).abs()/star.radius_m<1e-10,
            "orientation changed the reference radial mass quantiles");
        for (law,eta) in eta_laws {
            let mut p=initial.clone();
            let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
                .with_cubic_gravity();
            cfg.eta=eta;
            let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
            let a0=sph_force_audit::frozen_force_audit(
                &p,&mut state,physics,cfg,&star,shell_count);
            append_snapshot(&mut summary,&mut shell_csv,
                n,seed,law,eta,"initial",0.0,star.radius_m,
                expected_sampling_bias_pct,0.0,&seed_profile,&a0);
            let relaxed=relax_isolated_adaptive(
                &mut p,physics,&star,cfg,&mut state,
                RelaxParameters {
                    iterations:relax_steps,dt_s:dt,damping_time_s:0.7,
                }).expect("finite relaxation");
            assert!(relaxed.report.mass_relative_error<1e-12);
            let post_profile=sph_massprofile::profile_diagnostics(&p,&star);
            let a1=sph_force_audit::frozen_force_audit(
                &p,&mut state,physics,cfg,&star,shell_count);
            let contraction_pct=100.0*(
                post_profile.radius_max_m/seed_profile.radius_max_m-1.0);
            append_snapshot(&mut summary,&mut shell_csv,
                n,seed,law,eta,"post_damping",relax_steps as f64*dt,
                star.radius_m,expected_sampling_bias_pct,contraction_pct,
                &post_profile,&a1);
            println!(
                "seed={seed}, {law}, eta={eta:.4}, initial bias={expected_sampling_bias_pct:+.2}%, post bias={:+.2}%, contraction={contraction_pct:+.2}%, post P/gref={:.3}, G/gref={:.3}, P-err={:.3}, G-err={:.3}, net={:.3}",
                100.0*(post_profile.radius_max_m/star.radius_m-1.0),
                a1.global.pressure_radial_ratio,
                a1.global.gravity_radial_ratio,
                a1.global.pressure_ref_rms,
                a1.global.gravity_ref_rms,
                a1.global.net_ref_rms,
            );
        }
    }
    if let Some(path)=value(&args,"--csv") {
        std::fs::write(path,summary).expect("write summary CSV");
        println!("wrote {path}");
    }
    if let Some(path)=value(&args,"--shell-csv") {
        std::fs::write(path,shell_csv).expect("write shell CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn growth_eta_matches_phy14_exponent() {
        assert!((growth_eta(256)-1.05).abs()<1e-12);
        assert!((growth_eta(1024)/growth_eta(256)-2.0_f64.powf(1.0/3.0)).abs()<1e-12);
    }
}
