//! BG-PHY-25: frozen free-surface reconstruction screening.
//!
//! Every candidate uses particle geometry, SPH volumes and the unchanged
//! cold Fermi EOS. The independent 1D stellar profile is used only
//! to *score* density, enthalpy-gradient and force error.
//! New candidate forces are computed via central differences of the
//! entire re-solved scalar cold-EOS plus Poisson-gravity potential.
//! No candidate is installed in the live integrator.
//!
//! cargo run --release --bin phy25_bench -- --n 64 --seeds 0
//!   --layout glass --stage initial --kernel both
//!   --summary bg-phy25-summary.csv --shells-csv bg-phy25-shells.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::AdaptiveConfig;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_rk_compatible::DensityChoice;
use blighted_galaxy::sph_surface_gravity::HamiltonianPilot;
use blighted_galaxy::sph_boundary_potential::{self,Candidate};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass,WhiteDwarfProfile};
use std::fmt::Write as _;
use std::time::Instant;
fn arg<'a>(a:&'a [String],flag:&str)->Option<&'a str>{
    a.windows(2).find(|x|x[0]==flag).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(a:&[String],flag:&str,default:T)->T{
    arg(a,flag).and_then(|v|v.parse().ok()).unwrap_or(default)
}
fn ids(a:&[String])->Vec<u64>{
    let v=arg(a,"--seeds").unwrap_or("0,17").split(',')
        .map(|x|x.parse::<u64>().expect("invalid seed"))
        .collect::<Vec<_>>();
    assert!(!v.is_empty()&&v.len()<=8);
    assert!(v.windows(2).all(|x|x[0]<x[1]));
    v
}
fn num(x:f64)->String{format!("{x:.10e}")}
fn maybe(x:Option<f64>)->String{x.map(num).unwrap_or_default()}
fn csv_row(buffer:&mut String,cols:&[String]){
    let _=writeln!(buffer,"{}",cols.join(","));
}
fn kernel_config(star:&WhiteDwarfProfile,kernel:&str)->AdaptiveConfig {
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity().with_newton_density();
    if kernel=="c4" {cfg=cfg.with_wendland_c4();}
    cfg.eta=if kernel=="m4"{1.2}else{1.7};
    cfg
}
fn source_geometry(
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,source:&str,steps:usize,dt:f64,tau:f64,
)->Vec<Fluid3>{
    if source=="initial"{return p.to_vec();}
    let choice=if source=="m4_original"||source=="c4_original"{
        DensityChoice::Original
    }else{DensityChoice::SurfaceBlended};
    let mut pilot=HamiltonianPilot{
        particles:p.to_vec(),physics,cfg,
        reference_radius_m:star.radius_m,density:choice,
    };
    for _ in 0..steps {pilot.damped_advance(dt,tau);}
    pilot.particles
}
fn evaluate(
    summary:&mut String,shells:&mut String,
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    cfg:AdaptiveConfig,n:usize,seed:u64,layout:&str,
    stage:&str,source:&str,kernel:&str,choice:Candidate,
    intervals:usize,lower:f64,upper:f64,epsilon:f64,
) {
    let started=Instant::now();
    let batch=sph_boundary_potential::state(
        p,star,physics,cfg,choice,true);
    let optimum=sph_boundary_potential::homology_minimum(
        p,star,physics,cfg,choice,lower,upper,intervals,epsilon);
    let shell=sph_boundary_potential::shell_radial_forces(
        p,star,physics,cfg,choice,8,epsilon);
    let q_last=shell.last().unwrap();
    let q_penultimate=&shell[shell.len()-2];
    let centroid={
        let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
        p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/mass
    };
    let rmax=p.iter().map(|x|(x.pos_m-centroid).len())
        .fold(0.0_f64,f64::max);
    let old_force_check=if n<=64{
        sph_boundary_potential::validated_control_error(
            p,star,physics,cfg,choice,epsilon)
    }else{None};
    if let Some(error)=old_force_check{
        assert!(error<3e-4,
            "old Hamiltonian force and scalar shell derivative differ: {error}");
    }
    let cols=vec![
        n.to_string(),seed.to_string(),layout.into(),stage.into(),
        source.into(),kernel.into(),choice.label().into(),
        num(cfg.eta),num(rmax/1000.0),
        num(100.0*(rmax/star.radius_m-1.0)),
        optimum.bracketed.to_string(),
        maybe(optimum.s_star),maybe(optimum.radius_bias_pct),
        optimum.curvature_positive.to_string(),
        num(q_last.pressure_ratio),num(q_last.gravity_inward_ratio),
        num(q_last.net_ratio),num(q_penultimate.pressure_ratio),
        num(q_penultimate.gravity_inward_ratio),
        num(q_penultimate.net_ratio),
        num(batch.last_shell_partition),
        num(batch.last_shell_normal_alignment),
        num(batch.last_shell_moment_ratio),
        num(batch.last_shell_partition_depth_h),
        num(batch.last_shell_moment_depth_h),
        num(batch.last_shell_partition_below_half_fraction),
        num(batch.last_shell_true_reference_depth_h),
        num(batch.last_shell_partition_depth_abs_error_h),
        num(batch.last_shell_moment_depth_abs_error_h),
        num(batch.density_global_log_rms),
        num(batch.density_outer_log_rms),
        num(batch.density_last_log_rms),
        num(batch.enthalpy_gradient_rms),
        num(batch.enthalpy_outer_gradient_rms),
        num(batch.fallback_fraction),
        num(batch.h_cap_fraction),num(batch.minimum_omega),
        num(batch.scalar_internal_j),num(batch.matched_gravity_j),
        num(batch.total_j),
        maybe(old_force_check),
        num(started.elapsed().as_secs_f64()),
    ];
    csv_row(summary,&cols);
    for row in &shell{
        csv_row(shells,&[
            n.to_string(),seed.to_string(),layout.into(),stage.into(),
            source.into(),kernel.into(),choice.label().into(),
            num(cfg.eta),row.shell.to_string(),
            num(row.mean_radius_km),
            num(row.pressure_ratio),
            num(row.gravity_inward_ratio),
            num(row.net_ratio),
            num(row.pressure_virial_j),
            num(row.gravity_virial_j),
        ]);
    }
    println!("BG-PHY-25 N={n} seed={seed} {layout:5} {stage:7} {source:12} {kernel:2} {:23} last P/G={:.3}/{:.3} penult P/G={:.3}/{:.3} s*={} R*={} rhoLast={:.3} Hgrad={:.3} fallback={:.1}%",
        choice.label(),q_last.pressure_ratio,q_last.gravity_inward_ratio,
        q_penultimate.pressure_ratio,q_penultimate.gravity_inward_ratio,
        optimum.s_star.map(|x|format!("{x:.5}"))
            .unwrap_or("unbracketed".into()),
        optimum.radius_bias_pct.map(|x|format!("{x:+.2}%"))
            .unwrap_or("unbracketed".into()),
        batch.density_last_log_rms,batch.enthalpy_gradient_rms,
        100.0*batch.fallback_fraction);
}
fn main(){
    let a=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&a,"--n",256);
    let seeds=ids(&a);
    let layout=arg(&a,"--layout").unwrap_or("glass");
    let stage=arg(&a,"--stage").unwrap_or("initial");
    let kernel=arg(&a,"--kernel").unwrap_or("both");
    let intervals:usize=opt(&a,"--intervals",8);
    let lower:f64=opt(&a,"--lower",0.75);
    let upper:f64=opt(&a,"--upper",1.30);
    let epsilon:f64=opt(&a,"--eps",1e-4);
    let damp_steps:usize=opt(&a,"--damp-steps",300);
    let dt:f64=opt(&a,"--dt",0.01);
    let tau:f64=opt(&a,"--damp-tau",0.7);
    assert!(n>=32&&n<=256&&n%8==0);
    assert!(["shell","glass","both"].contains(&layout));
    assert!(["initial","post","both"].contains(&stage));
    assert!(["m4","c4","both"].contains(&kernel));
    assert!(intervals>=3&&intervals<=24&&lower>=0.4&&upper<=1.5&&
        lower<1.0&&upper>1.0);
    assert!(epsilon>1e-7&&epsilon<0.01);
    assert!(damp_steps<=3000&&dt>0.0&&dt<=0.02&&tau>0.0);
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut summary=String::from(
        "N,seed,layout,stage,source,kernel,candidate,eta,Rmax_km,Rmax_bias_pct,homology_bracketed,stationary_scale,stationary_Rmax_bias_pct,stationary_curvature_positive,last_pressure_ratio,last_gravity_ratio,last_net_ratio,penultimate_pressure_ratio,penultimate_gravity_ratio,penultimate_net_ratio,last_partition,last_normal_radial_alignment,last_Q_over_S,last_partition_inferred_depth_h,last_moment_ratio_inferred_depth_h,last_partition_below_half_fraction,last_reference_depth_h,last_partition_depth_abs_error_h,last_moment_ratio_depth_abs_error_h,density_global_log_rms,density_outer_log_rms,density_last_log_rms,enthalpy_gradient_error,enthalpy_outer_gradient_error,fallback_fraction,h_cap_fraction,min_omega,cold_internal_energy_j,matched_gravity_energy_j,total_potential_j,validated_old_force_error,wall_s\n");
    let mut shell_csv=String::from(
        "N,seed,layout,stage,source,kernel,candidate,eta,shell,mean_radius_km,pressure_ratio,gravity_inward_ratio,net_ratio,pressure_virial_j,gravity_virial_j\n");
    let kernels=if kernel=="both"{vec!["m4","c4"]}else{vec![kernel]};
    let layouts=if layout=="both"{vec!["shell","glass"]}else{vec![layout]};
    let mut cases=0;
    for seed in seeds{
        let c=V3::new(1e8,0.0,0.0);
        let base=seed_cubic_shells_oriented(
            &star,c,V3::default(),n,seed).0;
        for geo in &layouts{
            let p=if *geo=="glass"{
                sph_glass::relax_antipodal_angles(
                    &base,&star,c,1.2,20).0
            }else{base.clone()};
            for k in &kernels{
                let cfg=kernel_config(&star,k);
                if stage=="initial"||stage=="both"{
                    for choice in Candidate::all(){
                        evaluate(&mut summary,&mut shell_csv,&p,
                            &star,physics,cfg,n,seed,geo,"initial",
                            "common_seed",k,choice,intervals,lower,upper,epsilon);
                        cases+=1;
                    }
                }
                if stage=="post"||stage=="both"{
                    for (source,density_mode) in [
                        ("original",DensityChoice::Original),
                        ("surface_mls",DensityChoice::SurfaceBlended),
                    ]{
                        let label=if *k=="m4" {
                            if density_mode==DensityChoice::Original{
                                "m4_original"
                            }else{"m4_blend"}
                        }else if density_mode==DensityChoice::Original{
                            "c4_original"
                        }else{"c4_blend"};
                        let relaxed=source_geometry(&p,&star,physics,cfg,
                            label,damp_steps,dt,tau);
                        for choice in Candidate::all(){
                            evaluate(&mut summary,&mut shell_csv,&relaxed,
                                &star,physics,cfg,n,seed,geo,
                                "post",source,k,choice,
                                intervals,lower,upper,epsilon);
                            cases+=1;
                        }
                    }
                }
            }
        }
    }
    println!("BG-PHY-25 completed {cases} frozen boundary comparisons");
    if let Some(path)=arg(&a,"--summary"){
        std::fs::write(path,summary).expect("write summary");
        println!("wrote {path}");
    }
    if let Some(path)=arg(&a,"--shells-csv"){
        std::fs::write(path,shell_csv).expect("write shell diagnostics");
        println!("wrote {path}");
    }
}
