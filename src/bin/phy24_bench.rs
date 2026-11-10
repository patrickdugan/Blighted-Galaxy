//! BG-PHY-24: homology energy, stationary radius and shell virial benchmark.
//!
//! Every case uses the independently integrated 0.65-Msun mu_e=2 star,
//! matched Poisson gravity, and the BG-PHY-23 conservative Hamiltonian.
//! Seed and post-3s-damping states are evaluated at fixed coordinates.
//! A bounded homology scan measures where the *discrete* energy prefers
//! the radius, with no EOS/force rescaling or dynamics intervention.
//!
//! cargo run --release --bin phy24_bench -- --n 256 --seeds 0,17
//!   --layout both --stage initial --kernel both --density both
//!   --csv bg-phy24-summary.csv --profile-csv bg-phy24-homology.csv
//!   --shell-csv bg-phy24-shells.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::AdaptiveConfig;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::sph_homology::{self,HomologyMinimum};
use blighted_galaxy::sph_rk_compatible::DensityChoice;
use blighted_galaxy::sph_surface_gravity::HamiltonianPilot;
use blighted_galaxy::white_dwarf::{DegenerateElectrons,WhiteDwarfProfile,profile_for_mass};
use std::fmt::Write as _;
use std::time::Instant;
fn arg<'a>(args:&'a [String],flag:&str)->Option<&'a str>{
    args.windows(2).find(|x|x[0]==flag).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],flag:&str,default:T)->T{
    arg(args,flag).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn flag(args:&[String],name:&str)->bool{args.iter().any(|x|x==name)}
fn seeds(args:&[String])->Vec<u64>{
    let list=arg(args,"--seeds").unwrap_or("0,17").split(',')
        .map(|v|v.parse().expect("seed must be u64"))
        .collect::<Vec<u64>>();
    assert!(!list.is_empty()&&list.len()<=8);
    assert!(list.windows(2).all(|w|w[0]<w[1]));
    list
}
fn number(x:f64)->String{format!("{x:.10e}")}
fn optional(x:Option<f64>)->String{x.map(number).unwrap_or_default()}
fn write(csv:&mut String,fields:&[String]){
    let _=writeln!(csv,"{}",fields.join(","));
}
fn cfg_for(star:&WhiteDwarfProfile,kernel:&str)->AdaptiveConfig{
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity().with_newton_density();
    if kernel=="c4"{cfg=cfg.with_wendland_c4();}
    cfg.eta=if kernel=="m4"{1.2}else{1.7};
    cfg
}
fn density_mode(name:&str)->DensityChoice{
    if name=="blend"{DensityChoice::SurfaceBlended}
    else {DensityChoice::Original}
}
#[derive(Clone,Copy)]
struct Source<'a>{
    n:usize,seed:u64,layout:&'a str,stage:&'a str,
    source_model:&'a str,kernel:&'a str,density:&'a str,
}
fn key(s:&Source<'_>)->Vec<String>{
    vec![s.n.to_string(),s.seed.to_string(),
        s.layout.into(),s.stage.into(),s.source_model.into(),
        s.kernel.into(),s.density.into()]
}
fn emit(
    summary:&mut String,profile:&mut String,shell_csv:&mut String,
    p:&[Fluid3],star:&WhiteDwarfProfile,physics:Sph3Config,
    source:&Source<'_>,cfg:AdaptiveConfig,
    lower:f64,upper:f64,step:f64,epsilon:f64,
) {
    let started=Instant::now();
    let mode=density_mode(source.density);
    let scan=sph_homology::scan(
        p,star,physics,cfg,mode,lower,upper,step,epsilon);
    let shell=sph_homology::shell_virial(p,star,physics,cfg,mode,8);
    let at1=sph_homology::point(
        p,star,physics,cfg,mode,1.0,epsilon);
    let radius=sph_massprofile::profile_diagnostics(p,star);
    let virial_p=shell.iter().map(|x|x.pressure_virial_j).sum::<f64>();
    let virial_g=shell.iter().map(|x|x.gravity_virial_j).sum::<f64>();
    let outer_pressure=0.5*(shell[6].pressure_support_ratio+
        shell[7].pressure_support_ratio);
    let outer_gravity=0.5*(shell[6].gravity_inward_ratio+
        shell[7].gravity_inward_ratio);
    let outer_partition=0.5*(shell[6].partition_mean+
        shell[7].partition_mean);
    let max_fd=scan.samples.iter().filter(|x|x.active_set_stable)
        .map(|x|x.derivative_relative_error)
        .fold(at1.derivative_relative_error,f64::max);
    let unstable=scan.samples.iter().filter(|x|!x.active_set_stable).count();
    let mut columns=key(source);
    columns.extend([
        number(cfg.eta),number(star.radius_m/1000.0),
        number(radius.radius_max_m/1000.0),
        number(100.0*(radius.radius_max_m/star.radius_m-1.0)),
        number(radius.radial_quantile_rms),number(lower),number(upper),
        number(step),scan.bracketed.to_string(),
        optional(scan.stationary_scale),
        optional(scan.stationary_rmax_bias_pct),
        scan.stationary_curvature_positive.to_string(),
        optional(scan.stationary_total_j),
        number(at1.internal_j),number(at1.gravity_j),
        number(at1.total_j),number(at1.analytic_dv_ds_j),
        number(at1.pressure_dv_ds_j),number(at1.gravity_dv_ds_j),
        number(-virial_p),number(-virial_g),
        number(outer_pressure),number(outer_gravity),
        number(outer_partition),number(at1.cap_fraction),
        number(at1.minimum_omega),
        number(at1.mls_fallback_fraction),
        number(max_fd),unstable.to_string(),
        number(scan.lower_energy_j),number(scan.upper_energy_j),
        number(started.elapsed().as_secs_f64),
    ]);
    write(summary,&columns);
    for pt in &scan.samples{
        let mut v=key(source);
        v.extend([
            number(cfg.eta),number(pt.scale),number(pt.rmax_km),
            number(pt.rmax_bias_pct),number(pt.internal_j),
            number(pt.gravity_j),number(pt.total_j),
            number(pt.pressure_dv_ds_j),number(pt.gravity_dv_ds_j),
            number(pt.analytic_dv_ds_j),number(pt.fd_dv_ds_j),
            number(pt.derivative_relative_error),number(pt.curvature_fd_j),
            pt.active_set_stable.to_string(),number(pt.cap_fraction),
            number(pt.minimum_omega),number(pt.mls_fallback_fraction),
            number(pt.net_force_ratio),number(pt.net_torque_ratio),
        ]);
        write(profile,&v);
    }
    for m in &shell{
        let mut v=key(source);
        v.extend([
            number(cfg.eta),m.index.to_string(),
            m.count.to_string(),number(m.mass_fraction_lo),
            number(m.mass_fraction_hi),number(m.mean_radius_km),
            number(m.pressure_support_ratio),
            number(m.gravity_inward_ratio),number(m.net_radial_ratio),
            number(m.pressure_tangent_rms),number(m.net_tangent_rms),
            number(m.pressure_virial_j),number(m.gravity_virial_j),
            number(m.total_virial_j),number(m.partition_mean),
            number(m.sp_density_log_bias_mean),
            number(m.h_capped_fraction),number(m.minimum_omega),
            number(m.mls_fallback_fraction),
        ]);
        write(shell_csv,&v);
    }
    println!("BG-PHY-24 N={} seed={} {:5} {:13} source={} {:2}/{:8} eta={:.2} R={:+.2}% s*={} R*(ref)={} virialP={:+.4e} virialG={:+.4e} outerP={:.3} outerG={:.3} outerS={:.3} maxFD={:.2e} active_changes={unstable}",
        source.n,source.seed,source.layout,source.stage,source.source_model,
        source.kernel,source.density,cfg.eta,at1.rmax_bias_pct,
        scan.stationary_scale.map(|x|format!("{x:.5}"))
            .unwrap_or("unbracketed".into()),
        scan.stationary_rmax_bias_pct.map(|x|format!("{x:+.2}%"))
            .unwrap_or("unbracketed".into()),
        virial_p,virial_g,outer_pressure,outer_gravity,
        outer_partition,max_fd);
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",256);
    let ids=seeds(&args);
    let layout=arg(&args,"--layout").unwrap_or("both");
    let stage=arg(&args,"--stage").unwrap_or("initial");
    let kernels=arg(&args,"--kernel").unwrap_or("both");
    let densities=arg(&args,"--density").unwrap_or("both");
    let lower:f64=opt(&args,"--s-lower",0.75);
    let upper:f64=opt(&args,"--s-upper",1.30);
    let step:f64=opt(&args,"--s-step",0.05);
    let epsilon:f64=opt(&args,"--fd-eps",1e-4);
    let damp_steps:usize=opt(&args,"--damp-steps",300);
    let dt:f64=opt(&args,"--dt",0.01);
    let tau:f64=opt(&args,"--damp-tau",0.7);
    let glass_sweeps:usize=opt(&args,"--glass-sweeps",20);
    assert!(n>=32&&n<=256&&n%8==0);
    assert!(["shell","glass","both"].contains(&layout));
    assert!(["initial","post","both"].contains(&stage));
    assert!(["m4","c4","both"].contains(&kernels));
    assert!(["original","blend","both"].contains(&densities));
    assert!((0.40..1.0).contains(&lower)&&upper>1.0&&upper<=1.5);
    assert!(step>0.0&&step<=0.2);
    assert!(epsilon>1e-7&&epsilon<1e-2);
    assert!(damp_steps<=30000&&dt>0.0&&dt<=0.02&&tau>0.0);
    assert!(glass_sweeps>0&&glass_sweeps<=100);
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let mut summary=String::from("N,seed,layout,stage,source_model,kernel,density,eta,Rref_km,Rmax_km,Rmax_bias_pct,mass_cdf_rms,scan_lower,scan_upper,scan_step,minimum_bracketed,stationary_homology_scale,stationary_Rmax_bias_pct,minimum_curvature_positive,stationary_total_energy_j,internal_energy_s1_j,gravity_energy_s1_j,total_energy_s1_j,total_dV_ds_s1_j,pressure_dV_ds_s1_j,gravity_dV_ds_s1_j,shell_pressure_dV_ds_s1_j,shell_gravity_dV_ds_s1_j,outer_two_pressure_support_ratio,outer_two_gravity_ratio,outer_two_partition_mean,h_capped_fraction,min_omega,mls_fallback_fraction,max_stable_analytic_fd_relative_error,active_set_transition_count,lower_energy_j,upper_energy_j,wall_s\n");
    let mut profile=String::from("N,seed,layout,stage,source_model,kernel,density,eta,homology_scale,Rmax_km,Rmax_bias_pct,cold_internal_energy_j,matched_gravity_energy_j,total_potential_j,pressure_dV_ds_j,gravity_dV_ds_j,total_analytic_dV_ds_j,independent_fd_dV_ds_j,relative_dV_ds_error,fd_energy_curvature_j,active_set_stable,h_cap_fraction,min_omega,mls_fallback_fraction,net_force_ratio,net_torque_ratio\n");
    let mut shells=String::from("N,seed,layout,stage,source_model,kernel,density,eta,shell,n_particles,mass_fraction_lo,mass_fraction_hi,mean_radius_km,pressure_support_ratio,gravity_inward_ratio,net_radial_ratio,pressure_tangent_rms,net_tangent_rms,pressure_virial_j,gravity_virial_j,total_virial_j,partition_mean,sp_density_log_bias_mean,h_cap_fraction,min_omega,mls_fallback_fraction\n");
    let mut cross=String::from("N,seed,layout,stage,source_model,source_kernel,source_density,target_kernel,target_density,eta,Rmax_km,total_homology_derivative_j,pressure_homology_derivative_j,gravity_homology_derivative_j,outer_two_pressure_support_ratio,outer_two_gravity_ratio,outer_two_net_radial_ratio,outer_partition_mean,outer_pressure_virial_j,outer_gravity_virial_j,active_set_stable\n");
    let mut total_cases=0;
    for seed in ids {
        let center=V3::new(1e8,0.0,0.0);
        let original=seed_cubic_shells_oriented(&star,
            center,V3::default(),n,seed).0;
        let geometries=if layout=="both" {
            vec!["shell","glass"]
        }else{vec![layout]};
        let chosen_kernels=if kernels=="both" {
            vec!["m4","c4"]
        }else{vec![kernels]};
        let chosen_densities=if densities=="both" {
            vec!["original","blend"]
        }else{vec![densities]};
        for geom in geometries {
            let seeded=if geom=="glass"{
                sph_glass::relax_antipodal_angles(
                    &original,&star,center,1.2,glass_sweeps).0
            } else {original.clone()};
            for &kernel in &chosen_kernels {
                for &density in &chosen_densities {
                    let cfg=cfg_for(&star,kernel);
                    let label=format!("{kernel}_{density}");
                    if stage=="initial"||stage=="both" {
                        let source=Source{
                            n,seed,layout:geom,stage:"initial",
                            source_model:"common_seed",kernel,density,
                        };
                        emit(&mut summary,&mut profile,&mut shells,
                            &seeded,&star,physics,&source,cfg,
                            lower,upper,step,epsilon);
                        total_cases+=1;
                    }
                    if stage=="post"||stage=="both" {
                        let mut pilot=HamiltonianPilot{
                            particles:seeded.clone(),physics,cfg,
                            reference_radius_m:star.radius_m,
                            density:density_mode(density),
                        };
                        for _ in 0..damp_steps{pilot.damped_advance(dt,tau);}
                        let source=Source{
                            n,seed,layout:geom,stage:"post_damping",
                            source_model:&label,kernel,density,
                        };
                        emit(&mut summary,&mut profile,&mut shells,
                            &pilot.particles,&star,physics,&source,cfg,
                            lower,upper,step,epsilon);
                        total_cases+=1;
                        if arg(&args,"--cross-csv").is_some() {
                            // Evaluate original and surface-gated density
                            // at identical post-damping coordinates.
                            for &target_density in &chosen_densities {
                                let mode=density_mode(target_density);
                                let q=sph_homology::point(
                                    &pilot.particles,&star,physics,cfg,
                                    mode,1.0,epsilon);
                                let shell=sph_homology::shell_virial(
                                    &pilot.particles,&star,physics,cfg,
                                    mode,8);
                                let outer=&shell[6..8];
                                let mean_pressure=outer.iter()
                                    .map(|x|x.pressure_support_ratio)
                                    .sum::<f64>()/2.0;
                                let mean_gravity=outer.iter()
                                    .map(|x|x.gravity_inward_ratio)
                                    .sum::<f64>()/2.0;
                                let mean_net=outer.iter()
                                    .map(|x|x.net_radial_ratio)
                                    .sum::<f64>()/2.0;
                                let mean_partition=outer.iter()
                                    .map(|x|x.partition_mean)
                                    .sum::<f64>()/2.0;
                                let cols=vec![
                                    n.to_string(),seed.to_string(),
                                    geom.into(),"post_damping".into(),
                                    label.clone(),kernel.into(),
                                    density.into(),kernel.into(),
                                    target_density.into(),number(cfg.eta),
                                    number(q.rmax_km),
                                    number(q.analytic_dv_ds_j),
                                    number(q.pressure_dv_ds_j),
                                    number(q.gravity_dv_ds_j),
                                    number(mean_pressure),
                                    number(mean_gravity),
                                    number(mean_net),
                                    number(mean_partition),
                                    number(outer.iter().map(|x|x.pressure_virial_j)
                                        .sum::<f64>()),
                                    number(outer.iter().map(|x|x.gravity_virial_j)
                                        .sum::<f64>()),
                                    q.active_set_stable.to_string(),
                                ];
                                write(&mut cross,&cols);
                            }
                        }
                    }
                }
            }
        }
    }
    println!("BG-PHY-24 finished {total_cases} homology studies, {} profile rows, {} shell rows",
        profile.lines().count()-1,shells.lines().count()-1);
    if let Some(path)=arg(&args,"--csv"){
        std::fs::write(path,summary).expect("write homology summary");
        println!("wrote {path}");
    }
    if let Some(path)=arg(&args,"--profile-csv"){
        std::fs::write(path,profile).expect("write homology scan CSV");
        println!("wrote {path}");
    }
    if let Some(path)=arg(&args,"--shell-csv"){
        std::fs::write(path,shells).expect("write shell force CSV");
        println!("wrote {path}");
    }
    if let Some(path)=arg(&args,"--cross-csv"){
        std::fs::write(path,cross).expect("write fixed-position density cross-audit");
        println!("wrote {path}");
    }
    // Read-only geometric comparisons, never a gameplay force switch.
    assert!(total_cases>0);
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn passed_defaults_select_all_controls(){
        assert_eq!(seeds(&[]),vec![0,17]);
        assert_eq!(number(0.5),"5.0000000000e-1");
        assert!(!flag(&[],"--unrelated"));
    }
}
