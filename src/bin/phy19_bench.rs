//! BG-PHY-19: fitted pressure and enthalpy gradients near a free surface.
//!
//! Frozen N256/1024/2048 factorial using C4 100/150/200 physical neighbors
//! and M4 eta1.20 baseline. Compares weighted linear and quadratic local
//! fits on simulated SPH vs independent hydrostatic pressure/enthalpy fields.
//! Only diagnostics; dynamics and physical acceptance remain unchanged.
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Sph3Config,V3};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_neighbor_calibration;
use blighted_galaxy::sph_force_audit;
use blighted_galaxy::sph_gradient_audit;
use blighted_galaxy::sph_quadratic_audit::{self,FIELD_LABELS,PolynomialMetrics};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
fn value<'a>(args:&'a [String],flag:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==flag).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],flag:&str,default:T)->T{
    value(args,flag).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn list<T:std::str::FromStr+PartialOrd>(
    args:&[String],flag:&str,default:&str)->Vec<T>{
    let out=value(args,flag).unwrap_or(default)
        .split(',').map(|x|x.parse().ok().expect("invalid value"))
        .collect::<Vec<T>>();
    assert!(!out.is_empty()&&out.len()<=8);
    assert!(out.windows(2).all(|x|x[0]<x[1]));
    out
}
fn f(v:f64)->String{format!("{v:.10e}")}
fn row(csv:&mut String,cols:&[String]) {
    let _=writeln!(csv,"{}",cols.join(","));
}
fn push_metrics(cols:&mut Vec<String>,m:&PolynomialMetrics){
    cols.push(f(m.usable_fraction));
    cols.push(f(m.max_linear_pivot_spread));
    cols.push(f(m.max_quadratic_pivot_spread));
    cols.push(f(m.density_ln_bias_mean));
    for x in m.vector_error {cols.push(f(x));}
    for x in m.radial_ratio {cols.push(f(x));}
}
fn header()->String{
    let mut parts=vec![
        "N".into(),"seed".into(),"layout".into(),"kernel".into(),
        "series".into(),"target_neighbors".into(),"eta".into(),
        "calibration_reachable".into(),"calibration_error".into(),
        "realized_neighbors".into(),"mean_h_km".into(),
        "original_conservative_pressure_error".into(),
        "original_outer_pressure_ratio".into(),
        "original_pressure_gradient_linear_error".into(),
        "original_pressure_gradient_moment_corrected_error".into(),
    ];
    for scope in ["global","inner","outer"] {
        for key in ["usable","max_linear_pivot_spread",
                    "max_quadratic_pivot_spread","mean_density_ln_bias"] {
            parts.push(format!("{scope}_{key}"));
        }
        for kind in FIELD_LABELS{
            parts.push(format!("{scope}_{kind}_vector_error"));
        }
        for kind in FIELD_LABELS{
            parts.push(format!("{scope}_{kind}_radial_ratio"));
        }
    }
    format!("{}\n",parts.join(","))
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",1024);
    let seeds=list::<u64>(&args,"--seeds","0,17");
    let targets=list::<f64>(&args,"--targets","100,150,200");
    let sweeps:usize=opt(&args,"--glass-sweeps",20);
    let calibration_steps:usize=opt(&args,"--calibration-steps",12);
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(sweeps>0&&sweeps<=100);
    assert!(calibration_steps>=1&&calibration_steps<=20);
    assert!(targets.iter().all(|t|*t>=2.0&&*t<=n as f64));
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut csv=header();
    for seed in seeds {
        let center=V3::new(1e8,0.0,0.0);
        let source=seed_cubic_shells_oriented(
            &star,center,V3::default(),n,seed).0;
        let glass=sph_glass::relax_antipodal_angles(
            &source,&star,center,1.2,sweeps).0;
        for (layout,p) in [
            ("shell",source.as_slice()),("glass",glass.as_slice())
        ]{
            for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
                let mut base=AdaptiveConfig::from_star_radius(star.radius_m)
                    .with_cubic_gravity();
                if kernel==DensityKernel::WendlandC4 {
                    base=base.with_wendland_c4();
                }
                let mut cases=vec![("eta120",0.0,1.2,true,0.0)];
                if kernel==DensityKernel::WendlandC4{
                    for &target in &targets{
                        let cal=sph_neighbor_calibration::calibrate(
                            p,star.radius_m,base,target,0.75,2.85,
                            calibration_steps);
                        cases.push(("target",target,cal.eta,
                            cal.bracketed,cal.abs_error));
                    }
                }
                for (series,target,eta,reachable,cal_error) in cases{
                    let mut cfg=base.with_newton_density();
                    cfg.eta=eta;
                    let mut state=AdaptiveState::new(p,star.radius_m,cfg);
                    let forces=sph_force_audit::frozen_force_audit(
                        p,&mut state,physics,cfg,&star,8);
                    let existing=sph_gradient_audit::audit(
                        p,&state,cfg.kernel,&star);
                    let polynomial=sph_quadratic_audit::audit(
                        p,&state,cfg.kernel,&star);
                    let ngb=sph_neighbor_calibration::measure(p,&state,cfg);
                    let mut cols=vec![
                        n.to_string(),seed.to_string(),layout.into(),
                        kernel.label().into(),series.into(),f(target),f(eta),
                        reachable.to_string(),f(cal_error),
                        f(ngb.mean),f(ngb.mean_h_m/1000.0),
                        f(forces.global.pressure_ref_rms),
                        f(forces.shells.last().unwrap().metrics.pressure_radial_ratio),
                        f(existing.global.direct_num_error),
                        f(existing.global.corrected_num_error),
                    ];
                    push_metrics(&mut cols,&polynomial.global);
                    push_metrics(&mut cols,&polynomial.inner);
                    push_metrics(&mut cols,&polynomial.outer);
                    row(&mut csv,&cols);
                    println!("N={n} seed={seed} {layout:5} {:11} {series:6} ngb={:.2} Pstd={:.4} quadP={:.4} quadPref={:.4} quadH={:.4} quadHref={:.4} outerHref ratio={:.3} usable={:.1}%",
                        kernel.label(),ngb.mean,
                        forces.global.pressure_ref_rms,
                        polynomial.global.vector_error[1],
                        polynomial.global.vector_error[3],
                        polynomial.global.vector_error[5],
                        polynomial.global.vector_error[7],
                        polynomial.outer.radial_ratio[7],
                        100.0*polynomial.global.usable_fraction);
                }
            }
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write quadratic gradient CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn csv_header_has_all_four_field_pairs(){
        let h=header();
        assert!(h.contains("outer_enthalpy_ref_quadratic_radial_ratio"));
        assert!(h.contains("global_pressure_num_linear_vector_error"));
    }
}
