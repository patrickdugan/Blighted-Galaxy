//! BG-PHY-16: controlled 2x2 seed/kernel frozen-force experiment.
//!
//! N and orientation seed held fixed; compare legacy shell directions
//! against angular glass, and cubic M4 against Poisson-matched Wendland C4.
//! Each eta also compares the fixed and growing-neighbor laws.
//! By default ONLY frozen snapshots are measured; optional --relax-steps
//! tests matched dynamical behavior without weakening Issue #8 acceptance.
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass::{self,AngularGlassReport};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::sph_force_audit::{self,FrozenForceAudit};
use blighted_galaxy::sph_neighbors::CellNeighbors;
use blighted_galaxy::stellar_relax::{relax_isolated_adaptive,RelaxParameters};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use blighted_galaxy::sph3::Fluid3;
use std::fmt::Write as _;

fn arg<'a>(args:&'a [String],key:&str)->Option<&'a str> {
    args.windows(2).find(|x|x[0]==key).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T {
    arg(args,key).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn seeds(args:&[String])->Vec<u64> {
    let values=arg(args,"--seeds").unwrap_or("0,17");
    let list=values.split(',').map(|x|x.parse::<u64>()
        .expect("invalid seed")).collect::<Vec<_>>();
    assert!(!list.is_empty()&&list.len()<=8);
    assert!(list.windows(2).all(|w|w[0]<w[1]));
    list
}
fn growth_eta(n:usize)->f64 {1.05*(n as f64/256.0).powf(1.0/6.0)}
fn f(x:f64)->String {format!("{x:.9e}")}
fn write(csv:&mut String,cols:&[String]) {
    let _=writeln!(csv,"{}",cols.join(","));
}
fn actual_neighbors(p:&[Fluid3],state:&AdaptiveState)->f64 {
    let hmax=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let c=CellNeighbors::new(p,hmax).all_candidates();
    let count=c.iter().enumerate().map(|(i,js)|{
        js.iter().filter(|&&j|
            (p[i].pos_m-p[j].pos_m).len_sq()<
            (2.0*state.smoothing_m[i]).powi(2)).count()
    }).sum::<usize>();
    count as f64/p.len() as f64
}
fn emit(
    summary:&mut String,shells:&mut String,
    n:usize,seed:u64,law:&str,eta:f64,
    layout:&str,kernel:DensityKernel,stage:&str,elapsed_s:f64,
    reference_radius:f64,seed_radius:f64,
    optimization:AngularGlassReport,
    p:&[Fluid3],state:&AdaptiveState,
    profile:&MassProfileDiagnostic,a:&FrozenForceAudit,
) {
    let m=a.global;
    let contraction_pct=100.0*(profile.radius_max_m/seed_radius-1.0);
    let bias_pct=100.0*(profile.radius_max_m/reference_radius-1.0);
    let initial_sampling_pct=100.0*(seed_radius/reference_radius-1.0);
    let mean_h=state.smoothing_m.iter().sum::<f64>()/n as f64;
    write(summary,&[
        n.to_string(),seed.to_string(),law.to_string(),f(eta),
        layout.to_string(),kernel.label().to_string(),stage.to_string(),
        f(elapsed_s),f(profile.radius_max_m/1000.0),
        f(bias_pct),f(initial_sampling_pct),f(contraction_pct),
        f(profile.radial_quantile_rms),f(mean_h/1000.0),
        f(actual_neighbors(p,state)),
        f(optimization.starting_objective),f(optimization.ending_objective),
        optimization.accepted_sweeps.to_string(),
        f(optimization.maximum_rotation_degrees),
        f(m.pressure_ref_rms),f(m.gravity_ref_rms),f(m.net_ref_rms),
        f(m.pressure_radial_ratio),f(m.gravity_radial_ratio),
        f(m.pressure_tangent_rms),f(m.gravity_tangent_rms),f(m.net_tangent_rms),
        f(a.shells[0].metrics.pressure_ref_rms),
        f(a.shells[0].metrics.pressure_tangent_rms),
        f(a.shells.last().unwrap().metrics.pressure_radial_ratio),
        f(a.shells.last().unwrap().metrics.gravity_radial_ratio),
    ]);
    for (index,s) in a.shells.iter().enumerate() {
        let m=s.metrics;
        write(shells,&[
            n.to_string(),seed.to_string(),law.to_string(),f(eta),
            layout.to_string(),kernel.label().to_string(),
            stage.to_string(),index.to_string(),
            f(s.mass_fraction_lo),f(s.mass_fraction_hi),
            m.particle_count.to_string(),f(m.mean_radius_m/1000.0),
            f(m.pressure_ref_rms),f(m.gravity_ref_rms),f(m.net_ref_rms),
            f(m.pressure_radial_ratio),f(m.gravity_radial_ratio),
            f(m.pressure_tangent_rms),f(m.gravity_tangent_rms),
            f(m.net_tangent_rms),
        ]);
    }
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",128);
    let sweeps:usize=opt(&args,"--glass-sweeps",20);
    let relax_steps:usize=opt(&args,"--relax-steps",0);
    let dt:f64=opt(&args,"--dt",0.01);
    let shells:usize=opt(&args,"--shells",8);
    let fixed_eta:f64=opt(&args,"--fixed-eta",1.2);
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(sweeps>0&&sweeps<=200);
    assert!(relax_steps<=30000);
    assert!(dt>0.0&&dt<=0.1&&dt.is_finite());
    assert!(shells>=2&&shells<=n/4);
    assert!(fixed_eta>0.5&&fixed_eta<3.0);
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let seed_radius=star.radius_at_mass_fraction(1.0-1.0/n as f64);
    let mut summary=String::from("N,seed,eta_law,eta,layout,kernel,stage,elapsed_s,rmax_km,rmax_bias_pct,initial_sampling_pct,contraction_from_seed_pct,cdf_rms,mean_h_km,actual_neighbors,glass_start_objective,glass_end_objective,accepted_glass_sweeps,max_glass_rotation_degrees,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms,net_tangent_rms,inner_shell_pressure_err,inner_shell_pressure_tangent_rms,outer_shell_pressure_ratio,outer_shell_gravity_ratio\n");
    let mut shell_csv=String::from("N,seed,eta_law,eta,layout,kernel,stage,shell,mass_fraction_lo,mass_fraction_hi,particle_count,mean_radius_km,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms,net_tangent_rms\n");
    println!("BG-PHY-16 N={n}, R_ref={:.3} km, initial sampled radius {:.3} km, glass sweeps={sweeps}, dynamic steps={relax_steps}",
        star.radius_m/1000.0,seed_radius/1000.0);
    for seed in seeds(&args) {
        let center=V3::new(1e8,0.0,0.0);
        let original=seed_cubic_shells_oriented(&star,center,
            V3::default(),n,seed).0;
        let original_radius=sph_massprofile::profile_diagnostics(&original,&star).radius_max_m;
        assert!((original_radius-seed_radius).abs()/star.radius_m<1e-10);
        for (law,eta) in [("fixed",fixed_eta),("growth",growth_eta(n))] {
            let (glass,report)=sph_glass::relax_antipodal_angles(
                &original,&star,center,eta,sweeps);
            let original_report=AngularGlassReport {
                starting_objective:report.starting_objective,
                ending_objective:report.starting_objective,
                accepted_sweeps:0,maximum_rotation_degrees:0.0,
            };
            for (layout,p0,optim) in [
                ("shell",original.as_slice(),original_report),
                ("glass",glass.as_slice(),report),
            ] {
                let radius0=sph_massprofile::profile_diagnostics(p0,&star).radius_max_m;
                assert!((radius0-seed_radius).abs()/star.radius_m<1e-10,
                    "angular glass changed radial mass CDF");
                for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4] {
                    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
                        .with_cubic_gravity();
                    cfg.eta=eta;
                    if kernel==DensityKernel::WendlandC4 {
                        cfg=cfg.with_wendland_c4();
                    }
                    let mut p=p0.to_vec();
                    let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
                    let profile0=sph_massprofile::profile_diagnostics(&p,&star);
                    let a0=sph_force_audit::frozen_force_audit(
                        &p,&mut state,physics,cfg,&star,shells);
                    emit(&mut summary,&mut shell_csv,n,seed,law,eta,layout,
                        kernel,"initial",0.0,star.radius_m,seed_radius,
                        optim,&p,&state,&profile0,&a0);
                    println!("seed={seed} {law:6} {layout:5} {:12} P-err={:.4} G-err={:.4} inner-P={:.4} outer-P={:.4} glass objective {:.4}->{:.4}",
                        kernel.label(),a0.global.pressure_ref_rms,
                        a0.global.gravity_ref_rms,
                        a0.shells[0].metrics.pressure_ref_rms,
                        a0.shells.last().unwrap().metrics.pressure_radial_ratio,
                        optim.starting_objective,optim.ending_objective);
                    if relax_steps>0 {
                        let relaxed=relax_isolated_adaptive(
                            &mut p,physics,&star,cfg,&mut state,
                            RelaxParameters {
                                iterations:relax_steps,dt_s:dt,
                                damping_time_s:0.7,
                            }).expect("finite seeded stellar relaxation");
                        assert!(relaxed.report.mass_relative_error<1e-12);
                        let post=sph_massprofile::profile_diagnostics(&p,&star);
                        let a1=sph_force_audit::frozen_force_audit(
                            &p,&mut state,physics,cfg,&star,shells);
                        emit(&mut summary,&mut shell_csv,n,seed,law,eta,layout,
                            kernel,"post_damping",relax_steps as f64*dt,
                            star.radius_m,seed_radius,optim,
                            &p,&state,&post,&a1);
                        println!("    post contraction={:+.2}% net={:.4} P={:.4} G={:.4}",
                            100.0*(post.radius_max_m/seed_radius-1.0),
                            a1.global.net_ref_rms,
                            a1.global.pressure_ref_rms,
                            a1.global.gravity_ref_rms);
                    }
                }
            }
        }
    }
    if let Some(path)=arg(&args,"--csv"){
        std::fs::write(path,summary).expect("write BG-PHY-16 summary");
        println!("wrote {path}");
    }
    if let Some(path)=arg(&args,"--shell-csv"){
        std::fs::write(path,shell_csv).expect("write BG-PHY-16 shells");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn compact_experiment_has_expected_comparison_counts() {
        let variants=2*2*2;
        assert_eq!(variants,8);
        assert!((growth_eta(256)-1.05).abs()<1e-12);
    }
}
