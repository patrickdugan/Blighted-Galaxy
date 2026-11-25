//! BG-PHY-17: high-neighbour Wendland C4 convergence / consistency sweep.
//!
//! The eta search CALIBRATES the mean realized physical support neighbour
//! count using independently safeguarded exact h/rho roots on each frozen
//! seed layout. It then compares underrelaxed and exact forces at the same
//! particle positions, recording E0/E1 consistency and outer-shell support.
//! The optional M4 control is calibrated to the same realized neighbour
//! targets. All snapshots are frozen; physical radius and hold acceptance
//! remain separate in Issue #8.
//!
//! cargo run --release --bin phy17_bench -- --n 1024 --seeds 0,17
//!   --targets 100,150,200 --calibration-steps 12
//!   --csv bg-phy17-n1024.csv --shell-csv bg-phy17-n1024-shells.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass::{self,AngularGlassReport};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::sph_force_audit::{self,FrozenForceAudit};
use blighted_galaxy::sph_neighbor_calibration::{self,NeighborStats,CalibratedEta};
use blighted_galaxy::sph_consistency::{self,ConsistencyMetrics};
use blighted_galaxy::sph_hroot;
use blighted_galaxy::sph_neighbors::CellNeighbors;
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass,WhiteDwarfProfile};
use std::time::Instant;
use std::fmt::Write as _;

fn value<'a>(args:&'a [String],flag:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==flag).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],flag:&str,default:T)->T{
    value(args,flag).and_then(|s|s.parse().ok()).unwrap_or(default)
}
fn numbers<T:std::str::FromStr+PartialOrd+Copy>(
    args:&[String],flag:&str,default:&str
)->Vec<T> {
    let values=value(args,flag).unwrap_or(default);
    let list=values.split(',').map(|x|x.parse::<T>()
        .ok().expect("invalid sweep number")).collect::<Vec<_>>();
    assert!(!list.is_empty()&&list.len()<=12);
    assert!(list.windows(2).all(|w|w[0]<w[1]),"must be increasing");
    list
}
fn growth_eta(n:usize)->f64{1.05*(n as f64/256.0).powf(1.0/6.0)}
fn num(v:f64)->String{format!("{v:.10e}")}
fn row(text:&mut String,cols:&[String]){
    let _=writeln!(text,"{}",cols.join(","));
}
fn max_uncapped_residual(p:&[Fluid3],state:&AdaptiveState,cfg:AdaptiveConfig)->f64{
    let lists=CellNeighbors::new(p,cfg.max_h_m).all_candidates();
    let mut max=0.0_f64;
    for i in 0..p.len(){
        let h=state.smoothing_m[i];
        if (h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8 ||
           (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8 {continue;}
        let residual=sph_hroot::closure_residual_with_kernel(
            i,p,&lists[i],cfg.eta,h,cfg.kernel);
        max=max.max(residual);
    }
    max
}
struct Experiment<'a>{
    n:usize,seed:u64,layout:&'a str,kernel:DensityKernel,
    label:&'a str,target:f64,calibration:Option<CalibratedEta>,
    eta:f64,glass:AngularGlassReport,
}
fn audit_case(
    summary:&mut String,shell_csv:&mut String,
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    case:&Experiment<'_>,
    n_shells:usize,
){
    let started=Instant::now();
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity();
    if case.kernel==DensityKernel::WendlandC4 {
        cfg=cfg.with_wendland_c4();
    }
    cfg.eta=case.eta;
    let initial=sph_massprofile::profile_diagnostics(p,star);
    let expected_radius=star.radius_at_mass_fraction(1.0-1.0/case.n as f64);
    assert!((initial.radius_max_m-expected_radius).abs()/star.radius_m<1e-10);
    let bias=100.0*(initial.radius_max_m/star.radius_m-1.0);
    for (mode,cfg_mode) in [
        ("underrelaxed",cfg),
        ("exact_root",cfg.with_newton_density()),
    ] {
        let mut state=AdaptiveState::new(p,star.radius_m,cfg_mode);
        let f=sph_force_audit::frozen_force_audit(
            p,&mut state,physics,cfg_mode,star,n_shells);
        let support=sph_neighbor_calibration::measure(p,&state,cfg_mode);
        let consistency=sph_consistency::audit(p,&state,cfg_mode.kernel);
        let closure=max_uncapped_residual(p,&state,cfg_mode);
        let elapsed=started.elapsed().as_secs_f64();
        emit(summary,shell_csv,case,mode,
            &support,&consistency,&f,bias,closure,elapsed);
        if mode=="exact_root" {
            let g=f.global;
            let outer=f.shells.last().unwrap().metrics;
            println!("N={} seed={} {:5} {:11} {:8} target={:5.1} eta={:.5} neighbors={:.2} inner/outer={:.1}/{:.1} P={:.5} G={:.5} Pouter={:.4} E0={:.5} closure={:.2e} bracketed={}",
                case.n,case.seed,case.layout,case.kernel.label(),
                case.label,case.target,case.eta,support.mean,
                support.inner_mean,support.outer_mean,g.pressure_ref_rms,
                g.gravity_ref_rms,outer.pressure_radial_ratio,
                consistency.e0_rms,closure,
                case.calibration.map(|c|c.bracketed).unwrap_or(true));
        }
    }
}
fn emit(
    summary:&mut String,shell_csv:&mut String,
    case:&Experiment<'_>,mode:&str,
    support:&NeighborStats,consistency:&ConsistencyMetrics,
    audit:&FrozenForceAudit,bias:f64,closure:f64,seconds:f64,
){
    let g=audit.global;
    let inner=audit.shells.first().unwrap().metrics;
    let outer=audit.shells.last().unwrap().metrics;
    let cal=case.calibration;
    let reachable=cal.map(|c|c.bracketed).unwrap_or(true);
    let error=cal.map(|c|c.abs_error).unwrap_or(0.0);
    let attempts=cal.map(|c|c.evaluations).unwrap_or(0);
    row(summary,&[
        case.n.to_string(),case.seed.to_string(),
        case.layout.to_string(),case.kernel.label().to_string(),
        case.label.to_string(),num(case.target),num(case.eta),
        mode.to_string(),reachable.to_string(),num(error),
        attempts.to_string(),
        num(support.mean),support.min.to_string(),support.max.to_string(),
        num(support.inner_mean),num(support.outer_mean),
        num(support.mean_h_m/1000.0),
        num(support.lower_h_caps),num(support.upper_h_caps),
        num(support.min_omega),num(closure),
        num(g.pressure_ref_rms),num(g.gravity_ref_rms),num(g.net_ref_rms),
        num(g.pressure_radial_ratio),num(g.gravity_radial_ratio),
        num(g.pressure_tangent_rms),num(g.gravity_tangent_rms),
        num(inner.pressure_ref_rms),num(inner.pressure_tangent_rms),
        num(outer.pressure_ref_rms),num(outer.pressure_radial_ratio),
        num(outer.gravity_radial_ratio),
        num(consistency.e0_rms),num(consistency.e0_inner_rms),
        num(consistency.e0_outer_rms),
        num(consistency.e0_outer_mean_radial),
        num(consistency.e0_inner_tangent_rms),
        num(consistency.e1_trace_mean),num(consistency.e1_trace_outer_mean),
        num(bias),num(case.glass.starting_objective),
        num(case.glass.ending_objective),
        case.glass.accepted_sweeps.to_string(),
        num(case.glass.maximum_rotation_degrees),
        num(seconds),
    ]);
    for (s,shell) in audit.shells.iter().enumerate() {
        let m=shell.metrics;
        row(shell_csv,&[
            case.n.to_string(),case.seed.to_string(),
            case.layout.to_string(),case.kernel.label().to_string(),
            case.label.to_string(),num(case.target),num(case.eta),
            mode.to_string(),s.to_string(),
            num(shell.mass_fraction_lo),num(shell.mass_fraction_hi),
            m.particle_count.to_string(),num(m.mean_radius_m/1000.0),
            num(m.pressure_ref_rms),num(m.gravity_ref_rms),num(m.net_ref_rms),
            num(m.pressure_radial_ratio),num(m.gravity_radial_ratio),
            num(m.pressure_tangent_rms),num(m.gravity_tangent_rms),
        ]);
    }
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",1024);
    let seeds=numbers::<u64>(&args,"--seeds","0,17");
    let targets=numbers::<f64>(&args,"--targets","100,150,200");
    let calibration_steps:usize=opt(&args,"--calibration-steps",12);
    let glass_sweeps:usize=opt(&args,"--glass-sweeps",20);
    let eta_min:f64=opt(&args,"--eta-min",0.75);
    let eta_max:f64=opt(&args,"--eta-max",2.85);
    let n_shells:usize=opt(&args,"--shells",8);
    let c4_only=args.iter().any(|s|s=="--c4-only");
    assert!(n>=32&&n<=4096&&n%2==0);
    assert!(n_shells>=2&&n_shells<=n/4);
    assert!(glass_sweeps>0&&glass_sweeps<=100);
    assert!(calibration_steps>=1&&calibration_steps<=20);
    assert!(eta_min>0.5&&eta_min<eta_max&&eta_max<3.0);
    assert!(targets.iter().all(|x|*x>=2.0&&*x<=n as f64));
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut summary=String::from("N,seed,layout,kernel,series,target_neighbors,eta,closure_mode,calibration_bracketed,calibration_abs_error,calibration_evaluations,actual_neighbors,neighbor_min,neighbor_max,inner_neighbors,outer_neighbors,mean_h_km,lower_h_cap_fraction,upper_h_cap_fraction,min_omega,max_free_closure_residual,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms,inner_pressure_ref_rms,inner_pressure_tangent_rms,outer_pressure_ref_rms,outer_pressure_radial_ratio,outer_gravity_radial_ratio,e0_rms,e0_inner_rms,e0_outer_rms,e0_outer_mean_radial,e0_inner_tangent_rms,e1_trace_mean,e1_trace_outer_mean,initial_sampling_bias_pct,glass_start_objective,glass_end_objective,accepted_glass_sweeps,max_glass_rotation_degrees,elapsed_s\n");
    let mut shell_csv=String::from("N,seed,layout,kernel,series,target_neighbors,eta,closure_mode,shell,mass_fraction_lo,mass_fraction_hi,particle_count,mean_radius_km,pressure_ref_rms,gravity_ref_rms,net_ref_rms,pressure_radial_ratio,gravity_radial_ratio,pressure_tangent_rms,gravity_tangent_rms\n");
    println!("BG-PHY-17: N={n}, target counts={targets:?}, seeds={seeds:?}, glass sweeps={glass_sweeps}, eta range {eta_min}..{eta_max}");
    for seed in seeds {
        let center=V3::new(1e8,0.0,0.0);
        let original=seed_cubic_shells_oriented(&star,center,
            V3::default(),n,seed).0;
        // Fix geometry across eta and both kernels, avoiding a changing
        // angular glass objective as a confounder in the neighbour sweep.
        let (glass,optim)=sph_glass::relax_antipodal_angles(
            &original,&star,center,1.2,glass_sweeps);
        let original_optim=AngularGlassReport {
            starting_objective:optim.starting_objective,
            ending_objective:optim.starting_objective,
            accepted_sweeps:0,
            maximum_rotation_degrees:0.0,
        };
        for (layout,p,glass_report) in [
            ("shell",original.as_slice(),original_optim),
            ("glass",glass.as_slice(),optim),
        ] {
            let mut kernels=vec![DensityKernel::WendlandC4];
            if !c4_only {kernels.insert(0,DensityKernel::CubicM4);}
            for kernel in kernels {
                let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
                    .with_cubic_gravity();
                if kernel==DensityKernel::WendlandC4 {
                    cfg=cfg.with_wendland_c4();
                }
                for (series,eta) in [
                    ("fixed_120",1.2),
                    ("growth",growth_eta(n)),
                ] {
                    let case=Experiment {
                        n,seed,layout,kernel,label:series,
                        target:0.0,calibration:None,eta,glass:glass_report,
                    };
                    audit_case(&mut summary,&mut shell_csv,p,
                        &star,physics,&case,n_shells);
                }
                for target in targets.iter().copied() {
                    let cal=sph_neighbor_calibration::calibrate(
                        p,star.radius_m,cfg,target,
                        eta_min,eta_max,calibration_steps);
                    let case=Experiment {
                        n,seed,layout,kernel,label:"target",
                        target,calibration:Some(cal),eta:cal.eta,
                        glass:glass_report,
                    };
                    audit_case(&mut summary,&mut shell_csv,p,
                        &star,physics,&case,n_shells);
                }
            }
        }
    }
    if let Some(path)=value(&args,"--csv") {
        std::fs::write(path,summary).expect("write BG-PHY-17 CSV");
        println!("wrote {path}");
    }
    if let Some(path)=value(&args,"--shell-csv") {
        std::fs::write(path,shell_csv).expect("write BG-PHY-17 shell CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn eta_growth_reference() {
        assert!((growth_eta(256)-1.05).abs()<1e-12);
        assert!((growth_eta(1024)-1.05*4.0_f64.powf(1.0/6.0)).abs()<1e-12);
    }
    #[test]
    #[should_panic]
    fn duplicate_neighbor_targets_are_rejected() {
        let args=vec!["--targets".to_string(),"100,100".to_string()];
        numbers::<f64>(&args,"--targets","100,150,200");
    }
}
