//! BG-PHY-20: frozen density reconstruction and enthalpy-gradient audit.
//!
//! Compare SPH density, Shepard, local affine/quadratic MLS, and a
//! spherical mass-shell secant, at matched particle positions and h_i.
//! The independent 1D rho is used ONLY for scoring. Entropy/EOS, smoothing,
//! gravity, Omega, live forces and pairwise conservation are unchanged.
//!
//! cargo run --release --bin phy20_bench -- --n 1024 --seeds 0,17
//!   --targets 100,150,200 --csv bg-phy20-1024.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_neighbor_calibration;
use blighted_galaxy::sph_density_audit;
use blighted_galaxy::sph_quadratic_audit;
use blighted_galaxy::sph_force_audit;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::stellar_relax::{relax_isolated_adaptive,RelaxParameters};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass,WhiteDwarfProfile};
use std::fmt::Write as _;
fn value<'a>(args:&'a [String],key:&str)->Option<&'a str> {
    args.windows(2).find(|x|x[0]==key).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T{
    value(args,key).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn numbers<T:std::str::FromStr+PartialOrd>(
    args:&[String],key:&str,default:&str
)->Vec<T>{
    let results=value(args,key).unwrap_or(default).split(',')
        .map(|x|x.parse().ok().expect("invalid numeric value"))
        .collect::<Vec<T>>();
    assert!(!results.is_empty()&&results.len()<=8);
    assert!(results.windows(2).all(|x|x[0]<x[1]));
    results
}
fn f(v:f64)->String{format!("{v:.10e}")}
fn row(csv:&mut String,fields:&[String]){
    let _=writeln!(csv,"{}",fields.join(","));
}
#[derive(Clone,Copy,Default)]
struct DensityScores {
    mean_log:f64,
    rms_log:f64,
    mean_ratio:f64,
    mass:f64,
}
impl DensityScores {
    fn score(p:&[Fluid3],indices:&[usize],rho:&[f64],
        rho_ref:&[f64])->Self{
        let mut result=Self::default();
        let mut sum_sq=0.0;
        for &i in indices {
            let mass=p[i].mass_kg;
            let log=(rho[i]/rho_ref[i]).ln();
            assert!(log.is_finite());
            result.mass+=mass;
            result.mean_log+=mass*log;
            result.mean_ratio+=mass*rho[i]/rho_ref[i];
            sum_sq+=mass*log*log;
        }
        result.mean_log/=result.mass;
        result.rms_log=(sum_sq/result.mass).sqrt();
        result.mean_ratio/=result.mass;
        result
    }
}
fn metrics_csv(
    csv:&mut String,n:usize,seed:u64,layout:&str,stage:&str,
    kernel:DensityKernel,series:&str,target:f64,
    eta:f64,cal_error:f64,reachable:bool,mean_neighbors:f64,
    mean_h:f64,cdf_rms:f64,
    method:&sph_density_audit::DensityField,
    unity_all:f64,unity_outer:f64,
    global:DensityScores,inner:DensityScores,outer:DensityScores,
    baseline_force:f64,
    enthalpy:&sph_quadratic_audit::PolynomialAudit,
) {
    let p=&enthalpy.global;
    let o=&enthalpy.outer;
    row(csv,&[
        n.to_string(),seed.to_string(),layout.into(),stage.into(),
        kernel.label().into(),series.into(),f(target),f(eta),
        reachable.to_string(),f(cal_error),
        f(mean_neighbors),f(mean_h/1000.0),f(cdf_rms),
        method.method.name().into(),f(method.fallback_fraction),
        f(method.maximum_pivot_spread),
        f(unity_all),f(unity_outer),
        f(global.mean_log),f(global.rms_log),f(global.mean_ratio),
        f(inner.mean_log),f(inner.rms_log),f(inner.mean_ratio),
        f(outer.mean_log),f(outer.rms_log),f(outer.mean_ratio),
        f(baseline_force),f(p.usable_fraction),
        f(p.vector_error[4]),f(p.vector_error[5]),
        f(p.vector_error[7]),
        f(p.radial_ratio[5]),
        f(o.usable_fraction),f(o.vector_error[4]),f(o.vector_error[5]),
        f(o.vector_error[7]),f(o.radial_ratio[5]),
        f(o.radial_ratio[7]),f(o.density_ln_bias_mean),
    ]);
}
fn audit_case(
    csv:&mut String,p:&[Fluid3],star:&WhiteDwarfProfile,
    physics:Sph3Config,n:usize,seed:u64,layout:&str,stage:&str,
    kernel:DensityKernel,series:&str,target:f64,eta:f64,
    cal_error:f64,reachable:bool,
) {
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity();
    if kernel==DensityKernel::WendlandC4{
        cfg=cfg.with_wendland_c4();
    }
    cfg.eta=eta;
    cfg=cfg.with_newton_density();
    let mut state=AdaptiveState::new(p,star.radius_m,cfg);
    let conservative=sph_force_audit::frozen_force_audit(
        p,&mut state,physics,cfg,star,8);
    let support=sph_neighbor_calibration::measure(p,&state,cfg);
    let batch=sph_density_audit::estimate(p,&state,cfg.kernel);
    let total=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),|a,x|a+x.pos_m*x.mass_kg)/total;
    let mut ordered=(0..n).collect::<Vec<_>>();
    ordered.sort_by(|&i,&j|
        (p[i].pos_m-center).len_sq()
        .total_cmp(&(p[j].pos_m-center).len_sq()).then(i.cmp(&j)));
    let inner=&ordered[..n/4];
    let outer=&ordered[n-n/4..];
    let reference=p.iter().map(|x|star.density_at((x.pos_m-center).len()))
        .collect::<Vec<_>>();
    assert!(reference.iter().all(|x|x.is_finite()&&*x>0.0),
        "only seed configurations inside the original hydrostatic surface");
    let partition_all=batch.partition.iter().sum::<f64>()/n as f64;
    let partition_outer=outer.iter().map(|&i|batch.partition[i])
        .sum::<f64>()/outer.len() as f64;
    let cdf=sph_massprofile::profile_diagnostics(p,star).radial_quantile_rms;
    for field in &batch.fields {
        let scores=DensityScores::score(p,&ordered,&field.rho,&reference);
        let iscores=DensityScores::score(p,inner,&field.rho,&reference);
        let oscores=DensityScores::score(p,outer,&field.rho,&reference);
        let polynomial=sph_quadratic_audit::audit_with_density_samples(
            p,&state,&field.rho,cfg.kernel,star);
        metrics_csv(
            csv,n,seed,layout,stage,kernel,series,target,eta,cal_error,
            reachable,support.mean,support.mean_h_m,cdf,
            field,partition_all,partition_outer,scores,iscores,oscores,
            conservative.global.pressure_ref_rms,&polynomial);
        println!("N={n} seed={seed} {layout:5} {stage:15} {:11} {series:7} target={target:5.0} neighbors={:.1} {:18} rhoLN outer={:+.4} rhoLN RMS={:.4} Hquad={:.5} Hquad-outer={:.5} fallback={:.1}%",
            kernel.label(),support.mean,field.method.name(),
            oscores.mean_log,oscores.rms_log,
            polynomial.global.vector_error[5],
            polynomial.outer.vector_error[5],
            100.0*field.fallback_fraction);
    }
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",1024);
    let seeds=numbers::<u64>(&args,"--seeds","0,17");
    let targets=numbers::<f64>(&args,"--targets","100,150,200");
    let sweeps:usize=opt(&args,"--glass-sweeps",20);
    let rounds:usize=opt(&args,"--calibration-steps",12);
    let relax_steps:usize=opt(&args,"--relax-steps",0);
    let dt:f64=opt(&args,"--dt",0.01);
    let c4_only=args.iter().any(|x|x=="--c4-only");
    assert!(n>=32&&n<=2048&&n%2==0);
    assert!(targets.iter().all(|x|*x>=2.0&&*x<=n as f64));
    assert!(sweeps>=1&&sweeps<=100);
    assert!(rounds>=1&&rounds<=20);
    assert!(relax_steps<=3000&&dt>0.0&&dt<=0.1);
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut csv=String::from(
        "N,seed,layout,stage,kernel,series,target_neighbors,eta,calibration_reachable,calibration_abs_error,realized_neighbors,mean_h_km,cdf_rms,estimator,fallback_fraction,max_pivot_spread,partition_unity_global,partition_unity_outer,density_global_log_bias,density_global_log_rms,density_global_mean_ratio,density_inner_log_bias,density_inner_log_rms,density_inner_mean_ratio,density_outer_log_bias,density_outer_log_rms,density_outer_mean_ratio,original_conservative_pressure_error,polynomial_usable_fraction,enthalpy_linear_error,enthalpy_quadratic_error,reference_enthalpy_quadratic_error,enthalpy_quadratic_radial_ratio,outer_polynomial_usable_fraction,outer_enthalpy_linear_error,outer_enthalpy_quadratic_error,outer_reference_enthalpy_quadratic_error,outer_enthalpy_quadratic_radial_ratio,outer_reference_enthalpy_quadratic_radial_ratio,outer_polynomial_density_ln_bias\n"
    );
    println!("BG-PHY-20 free-surface density experiments N={n}, seeds={seeds:?}, target neighbors={targets:?}");
    for seed in seeds {
        let center=V3::new(1e8,0.0,0.0);
        let base=seed_cubic_shells_oriented(
            &star,center,V3::default(),n,seed).0;
        let glass=sph_glass::relax_antipodal_angles(
            &base,&star,center,1.2,sweeps).0;
        for (layout,initial) in [
            ("shell",base.as_slice()),("glass",glass.as_slice()),
        ] {
            let mut stages=vec![("initial",initial.to_vec())];
            if relax_steps>0 {
                let mut p=initial.to_vec();
                let m4=AdaptiveConfig::from_star_radius(star.radius_m)
                    .with_cubic_gravity();
                let mut state=AdaptiveState::new(&p,star.radius_m,m4);
                let result=relax_isolated_adaptive(
                    &mut p,physics,&star,m4,&mut state,
                    RelaxParameters{
                        iterations:relax_steps,dt_s:dt,damping_time_s:0.7,
                    }).expect("finite M4 preparation");
                assert!(result.report.mass_relative_error<1e-12);
                stages.push(("post_m4_damping",p));
            }
            for (stage,p) in stages {
                let mut kernels=vec![DensityKernel::WendlandC4];
                if !c4_only {kernels.insert(0,DensityKernel::CubicM4);}
                for kernel in kernels{
                    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
                        .with_cubic_gravity();
                    if kernel==DensityKernel::WendlandC4 {
                        cfg=cfg.with_wendland_c4();
                    }
                    let mut cases=vec![("eta120",0.0,1.2,0.0,true)];
                    if kernel==DensityKernel::WendlandC4 {
                        for &target in &targets{
                            let result=sph_neighbor_calibration::calibrate(
                                &p,star.radius_m,cfg,target,0.75,2.85,rounds);
                            cases.push(("target",target,result.eta,
                                result.abs_error,result.bracketed));
                        }
                    }
                    for (series,target,eta,error,bracketed) in cases {
                        audit_case(
                            &mut csv,&p,&star,physics,n,seed,layout,stage,
                            kernel,series,target,eta,error,bracketed);
                    }
                }
            }
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write density comparison CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn smoke_arguments_are_sorted(){
        assert_eq!(numbers::<u64>(&[],"--seeds","0,17"),vec![0,17]);
        assert_eq!(numbers::<f64>(&[],"--targets","100,150,200").len(),3);
    }
}
