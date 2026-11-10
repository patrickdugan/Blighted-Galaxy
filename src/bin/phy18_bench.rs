//! BG-PHY-18: non-dynamical Bonet-Lok pressure gradient audit.
//!
//! Compare direct SPH difference gradients with per-particle moment-inverse
//! linear-gradient estimators. Each method is applied to the numerical
//! density-derived pressure field and the independently sampled hydrostatic
//! pressure field. The latter separates particle-geometry and boundary
//! interpolation errors from density/EOS sampling bias.
//!
//! The original pairwise energy/momentum-conserving force model is unchanged.
//!
//! cargo run --release --bin phy18_bench -- --n 1024 --seeds 0,17
//!   --targets 100,150,200 --csv bg-phy18-n1024.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Sph3Config,V3};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_neighbor_calibration;
use blighted_galaxy::sph_consistency;
use blighted_galaxy::sph_force_audit;
use blighted_galaxy::sph_gradient_audit::{self,GradientMetrics};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;

fn value<'a>(args:&'a [String],key:&str)->Option<&'a str> {
    args.windows(2).find(|x|x[0]==key).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T {
    value(args,key).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn integers(args:&[String],key:&str,default:&str)->Vec<u64> {
    let numbers=value(args,key).unwrap_or(default).split(',')
        .map(|x|x.parse::<u64>().expect("invalid seed"))
        .collect::<Vec<_>>();
    assert!(!numbers.is_empty()&&numbers.len()<=8);
    assert!(numbers.windows(2).all(|x|x[0]<x[1]));
    numbers
}
fn targets(args:&[String])->Vec<f64> {
    let values=value(args,"--targets").unwrap_or("100,150,200");
    let nums=values.split(',').map(|x|x.parse::<f64>()
        .expect("invalid neighbor target")).collect::<Vec<_>>();
    assert!(!nums.is_empty()&&nums.len()<=6);
    assert!(nums.windows(2).all(|x|x[0]<x[1]));
    nums
}
fn f(x:f64)->String {format!("{x:.10e}")}
fn emit(csv:&mut String,n:usize,seed:u64,layout:&str,
    kernel:DensityKernel,series:&str,target:f64,
    eta:f64,reachable:bool,cal_error:f64,
    support:&sph_neighbor_calibration::NeighborStats,
    force:&sph_force_audit::FrozenForceAudit,
    e0:&sph_consistency::ConsistencyMetrics,
    m:&GradientMetrics,inner:&GradientMetrics,outer:&GradientMetrics) {
    let cols=[
        n.to_string(),seed.to_string(),layout.to_string(),
        kernel.label().to_string(),series.to_string(),f(target),
        f(eta),reachable.to_string(),f(cal_error),
        f(support.mean),f(support.inner_mean),f(support.outer_mean),
        f(support.mean_h_m/1000.0),
        f(support.min_omega),
        f(force.global.pressure_ref_rms),f(force.global.gravity_ref_rms),
        f(force.global.net_ref_rms),
        f(force.shells.last().unwrap().metrics.pressure_radial_ratio),
        f(e0.e0_rms),f(e0.e0_outer_rms),
        f(m.usable_fraction),f(m.moment_condition_max),
        f(m.moment_trace_mean),f(m.density_ln_ratio_mean),
        f(m.density_ln_ratio_rms),
        f(m.direct_num_error),f(m.corrected_num_error),
        f(m.direct_ref_error),f(m.corrected_ref_error),
        f(m.direct_num_radial_ratio),f(m.corrected_num_radial_ratio),
        f(m.direct_ref_radial_ratio),f(m.corrected_ref_radial_ratio),
        f(inner.direct_num_error),f(inner.corrected_num_error),
        f(inner.direct_ref_error),f(inner.corrected_ref_error),
        f(inner.density_ln_ratio_mean),f(inner.moment_trace_mean),
        f(outer.usable_fraction),f(outer.direct_num_error),
        f(outer.corrected_num_error),f(outer.direct_ref_error),
        f(outer.corrected_ref_error),
        f(outer.direct_num_radial_ratio),
        f(outer.corrected_num_radial_ratio),
        f(outer.direct_ref_radial_ratio),
        f(outer.corrected_ref_radial_ratio),
        f(outer.density_ln_ratio_mean),f(outer.density_ln_ratio_rms),
        f(outer.moment_trace_mean),f(outer.moment_condition_max),
        f(outer.corrected_num_tangent_rms),
        f(outer.corrected_ref_tangent_rms),
    ];
    let _=writeln!(csv,"{}",cols.join(","));
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",1024);
    let seeds=integers(&args,"--seeds","0,17");
    let target_list=targets(&args);
    let sweeps:usize=opt(&args,"--glass-sweeps",20);
    let steps:usize=opt(&args,"--calibration-steps",12);
    let c4_only=args.iter().any(|a|a=="--c4-only");
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(sweeps>=1&&sweeps<=100);
    assert!(steps>=1&&steps<=20);
    assert!(target_list.iter().all(|x|*x>=2.0&&*x<=n as f64));
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut csv=String::from("N,seed,layout,kernel,series,target_neighbors,eta,calibration_reachable,calibration_error,realized_neighbors,inner_neighbors,outer_neighbors,mean_h_km,min_omega,standard_pressure_vector_error,standard_gravity_vector_error,standard_net_vector_error,standard_outer_pressure_ratio,e0_rms,e0_outer_rms,usable_fraction,moment_condition_max,moment_trace_mean,density_ln_ratio_mean,density_ln_ratio_rms,direct_numeric_pressure_error,corrected_numeric_pressure_error,direct_oracle_pressure_error,corrected_oracle_pressure_error,direct_numeric_pressure_ratio,corrected_numeric_pressure_ratio,direct_oracle_pressure_ratio,corrected_oracle_pressure_ratio,inner_direct_numeric_error,inner_corrected_numeric_error,inner_direct_oracle_error,inner_corrected_oracle_error,inner_density_ln_ratio_mean,inner_moment_trace_mean,outer_usable_fraction,outer_direct_numeric_error,outer_corrected_numeric_error,outer_direct_oracle_error,outer_corrected_oracle_error,outer_direct_numeric_pressure_ratio,outer_corrected_numeric_pressure_ratio,outer_direct_oracle_pressure_ratio,outer_corrected_oracle_pressure_ratio,outer_density_ln_ratio_mean,outer_density_ln_ratio_rms,outer_moment_trace_mean,outer_moment_condition_max,outer_corrected_numeric_tangent_rms,outer_corrected_oracle_tangent_rms\n");
    println!("BG-PHY-18 pressure-gradient consistency N={n}, seeds={seeds:?}, targets={target_list:?}");
    for seed in seeds {
        let center=V3::new(1e8,0.0,0.0);
        let original=seed_cubic_shells_oriented(
            &star,center,V3::default(),n,seed).0;
        let glass=sph_glass::relax_antipodal_angles(
            &original,&star,center,1.2,sweeps).0;
        for (layout,p) in [
            ("shell",original.as_slice()),
            ("glass",glass.as_slice()),
        ] {
            for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4] {
                if c4_only&&kernel==DensityKernel::CubicM4 {continue;}
                let mut base=AdaptiveConfig::from_star_radius(star.radius_m)
                    .with_cubic_gravity();
                if kernel==DensityKernel::WendlandC4 {
                    base=base.with_wendland_c4();
                }
                let mut scenarios=vec![("eta120",0.0,1.2,true,0.0)];
                if kernel==DensityKernel::WendlandC4 {
                    for target in &target_list {
                        let cal=sph_neighbor_calibration::calibrate(
                            p,star.radius_m,base,*target,0.75,2.85,steps);
                        scenarios.push(("target",*target,cal.eta,
                            cal.bracketed,cal.abs_error));
                    }
                }
                for (series,target,eta,reachable,err) in scenarios {
                    let mut cfg=base.with_newton_density();
                    cfg.eta=eta;
                    let mut state=AdaptiveState::new(p,star.radius_m,cfg);
                    let forces=sph_force_audit::frozen_force_audit(
                        p,&mut state,physics,cfg,&star,8);
                    let ngb=sph_neighbor_calibration::measure(p,&state,cfg);
                    let e0=sph_consistency::audit(p,&state,cfg.kernel);
                    let gradient=sph_gradient_audit::audit(
                        p,&state,cfg.kernel,&star);
                    emit(&mut csv,n,seed,layout,kernel,series,target,eta,
                        reachable,err,&ngb,&forces,&e0,
                        &gradient.global,&gradient.inner,&gradient.outer);
                    println!("seed={seed:2} {layout:5} {:11} {series:6} neighbors={:.1} eta={eta:.4} Pstd={:.3} Pdiff={:.3} Pcorr={:.3} Poracle={:.3}->{:.3} outer correction {:.3}->{:.3} rho-ln={:+.3}, moment-trace={:.3} valid={:.0}%",
                        kernel.label(),ngb.mean,
                        forces.global.pressure_ref_rms,
                        gradient.global.direct_num_error,
                        gradient.global.corrected_num_error,
                        gradient.global.direct_ref_error,
                        gradient.global.corrected_ref_error,
                        gradient.outer.direct_ref_error,
                        gradient.outer.corrected_ref_error,
                        gradient.outer.density_ln_ratio_mean,
                        gradient.outer.moment_trace_mean,
                        100.0*gradient.global.usable_fraction);
                }
            }
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write BG-PHY-18 CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gradient_case_defaults_are_valid() {
        assert_eq!(targets(&[]).len(),3);
        assert_eq!(integers(&[],"--seeds","0,17"),vec![0,17]);
    }
}
