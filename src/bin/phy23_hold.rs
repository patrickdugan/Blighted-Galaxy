//! BG-PHY-23: opt-in self-gravitating cold white-dwarf physical pilot.
//!
//! The pilot uses exactly the independent 0.65 Msun, mu_e=2 reference,
//! three seconds of tau=0.7s constant-rate damping, then at least 3 t_dyn
//! undamped. Track radius/CDF and energy/momentum/angular-momentum drift.
//! All acceptance gates are reported independently of numeric CI success.
//! The shipped playable physics solver remains unchanged.
//!
//! cargo run --release --bin phy23_hold -- --n 256 --kernel c4 --density blend
//!   --relax-steps 300 --hold-steps 788 --dt 0.01 --csv bg-phy23-hold.csv
//!   --trace-csv bg-phy23-hold-trace.csv
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::AdaptiveConfig;
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_rk_compatible::DensityChoice;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::sph_surface_gravity::HamiltonianPilot;
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass,WhiteDwarfProfile};
use std::time::Instant;
use std::fmt::Write as _;
fn value<'a>(args:&'a [String],key:&str)->Option<&'a str> {
    args.windows(2).find(|v|v[0]==key).map(|v|v[1].as_str())
}
fn opt<T:std::str::FromStr>(a:&[String],key:&str,d:T)->T {
    value(a,key).and_then(|v|v.parse().ok()).unwrap_or(d)
}
fn f(x:f64)->String {format!("{x:.9e}")}
fn cross(a:V3,b:V3)->V3{
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn momentum(p:&[Fluid3])->V3 {
    p.iter().fold(V3::default(),|s,x|s+x.vel_m_s*x.mass_kg)
}
fn angular(p:&[Fluid3])->V3 {
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let c=p.iter().fold(V3::default(),|sum,x|sum+x.pos_m*x.mass_kg)/m;
    let v=momentum(p)*(1.0/m);
    p.iter().fold(V3::default(),|s,x|
        s+cross(x.pos_m-c,(x.vel_m_s-v)*x.mass_kg))
}
fn metric_radii(p:&[Fluid3],star:&WhiteDwarfProfile)->MassProfileDiagnostic{
    sph_massprofile::profile_diagnostics(p,star)
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",128);
    let seed:u64=opt(&args,"--seed",0);
    let kernel=value(&args,"--kernel").unwrap_or("c4");
    let density=value(&args,"--density").unwrap_or("blend");
    let layout=value(&args,"--layout").unwrap_or("glass");
    let eta:f64=opt(&args,"--eta",1.7);
    let relax_steps:usize=opt(&args,"--relax-steps",300);
    let hold_steps:usize=opt(&args,"--hold-steps",788);
    let dt:f64=opt(&args,"--dt",0.01);
    let tau:f64=opt(&args,"--damping-tau",0.7);
    let glass_sweeps:usize=opt(&args,"--glass-sweeps",20);
    assert!(n>=32&&n<=256&&n%2==0);
    assert!(kernel=="c4"||kernel=="m4");
    assert!(density=="blend"||density=="original");
    assert!(layout=="glass"||layout=="shell");
    assert!(relax_steps>=1&&relax_steps<=30000);
    assert!(hold_steps>=1&&hold_steps<=30000);
    assert!(dt>0.0&&dt<=0.02&&tau>0.0);
    assert!(eta>0.5&&eta<3.0);
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let original=seed_cubic_shells_oriented(&star,
        center,V3::default(),n,seed).0;
    let particles=if layout=="glass"{
        sph_glass::relax_antipodal_angles(
            &original,&star,center,1.2,glass_sweeps).0
    }else{original};
    let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
        .with_cubic_gravity().with_newton_density();
    if kernel=="c4" {cfg=cfg.with_wendland_c4();}
    cfg.eta=eta;
    let selected=if density=="blend"{
        DensityChoice::SurfaceBlended
    }else{DensityChoice::Original};
    let mut pilot=HamiltonianPilot {
        particles,physics,cfg,reference_radius_m:star.radius_m,
        density:selected,
    };
    let seed_profile=metric_radii(&pilot.particles,&star);
    let initial_center=seed_profile.center_m;
    let reference_mass=pilot.particles.iter().map(|x|x.mass_kg).sum::<f64>();
    let initial_momentum=momentum(&pilot.particles);
    let started=Instant::now();
    for _ in 0..relax_steps {
        pilot.damped_advance(dt,tau);
    }
    let post=metric_radii(&pilot.particles,&star);
    let post_force=pilot.gradient();
    let post_energy=pilot.energy_j();
    let post_momentum=momentum(&pilot.particles);
    let post_angular=angular(&pilot.particles);
    let dyn_time=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let full_hold=hold_steps as f64*dt>=3.0*dyn_time;
    let mut max_rmax=0.0_f64;
    let mut max_r99=0.0_f64;
    let mut max_cdf=post.radial_quantile_rms;
    let mut max_energy=0.0_f64;
    let mut max_momentum=0.0_f64;
    let mut max_angular=0.0_f64;
    let mut trace=String::from("N,seed,kernel,density,eta,layout,elapsed_s,rmax_km,r99_km,cdf_rms,rmax_drift_pct,r99_drift_pct,energy_drift,momentum_drift,angular_drift\n");
    let stride=(hold_steps/32).max(1);
    for step in 1..=hold_steps {
        pilot.advance(dt);
        if step%stride==0||step==hold_steps {
            let now=metric_radii(&pilot.particles,&star);
            let dr=(now.radius_max_m/post.radius_max_m-1.0).abs();
            let dr99=(now.r99_m/post.r99_m-1.0).abs();
            let energy=(pilot.energy_j()-post_energy).abs()/
                post_energy.abs().max(1.0);
            let momentum=(momentum(&pilot.particles)-post_momentum)
                .len()/(reference_mass*3e5);
            let angular=(angular(&pilot.particles)-post_angular)
                .len()/(reference_mass*star.radius_m*3e5);
            max_rmax=max_rmax.max(dr);
            max_r99=max_r99.max(dr99);
            max_cdf=max_cdf.max(now.radial_quantile_rms);
            max_energy=max_energy.max(energy);
            max_momentum=max_momentum.max(momentum);
            max_angular=max_angular.max(angular);
            assert!(now.radius_max_m.is_finite()&&
                now.r99_m.is_finite()&&now.radial_quantile_rms.is_finite());
            assert!(energy.is_finite()&&momentum.is_finite()&&
                angular.is_finite());
            let line=vec![
                n.to_string(),seed.to_string(),kernel.into(),density.into(),
                f(eta),layout.into(),f(step as f64*dt),
                f(now.radius_max_m/1000.0),f(now.r99_m/1000.0),
                f(now.radial_quantile_rms),
                f(100.0*dr),f(100.0*dr99),
                f(energy),f(momentum),f(angular),
            ];
            let _=writeln!(trace,"{}",line.join(","));
        }
    }
    let bias=100.0*(post.radius_max_m/star.radius_m-1.0);
    let contraction=100.0*(post.radius_max_m/seed_profile.radius_max_m-1.0);
    let mass_now=pilot.particles.iter().map(|x|x.mass_kg).sum::<f64>();
    let mass_drift=(mass_now-reference_mass).abs()/reference_mass;
    let center_shift=(post.center_m-initial_center).len()/star.radius_m;
    let initial_bulk=(post_momentum-initial_momentum).len()/(reference_mass*3e5);
    let physical_gate=full_hold&&bias.abs()<5.0&&
        max_rmax<0.05&&max_r99<0.05&&max_cdf<0.05&&
        max_energy<1e-3&&max_momentum<1e-8&&max_angular<1e-8;
    let elapsed=started.elapsed().as_secs_f64();
    let mut csv=String::from("N,seed,layout,kernel,density,eta,relax_s,hold_s,dynamical_time_s,full_hold,seed_rmax_km,post_rmax_km,post_r99_km,rmax_bias_pct,contraction_from_seed_pct,post_cdf_rms,post_energy_j,post_internal_j,post_gravity_j,post_net_force_ratio,post_net_torque_ratio,post_old_sph_force_discrepancy,post_capped_fraction,post_min_omega,worst_rmax_drift_pct,worst_r99_drift_pct,worst_mass_cdf_rms,max_energy_drift,max_momentum_drift,max_angular_momentum_drift,mass_drift,initial_bulk_momentum_drift,center_shift_ref_radius,physical_gate,wall_s\n");
    let line=vec![
        n.to_string(),seed.to_string(),layout.into(),
        kernel.into(),density.into(),f(eta),
        f(relax_steps as f64*dt),f(hold_steps as f64*dt),f(dyn_time),
        full_hold.to_string(),
        f(seed_profile.radius_max_m/1000.0),
        f(post.radius_max_m/1000.0),
        f(post.r99_m/1000.0),f(bias),f(contraction),
        f(post.radial_quantile_rms),f(post_energy),
        f(post_force.cold_internal_energy_j),
        f(post_force.gravity_potential_j),
        f(post_force.net_force_ratio),f(post_force.net_torque_ratio),
        f(post_force.old_sph_relative_force_error),
        f(post_force.capped_fraction),f(post_force.minimum_omega),
        f(100.0*max_rmax),f(100.0*max_r99),f(max_cdf),
        f(max_energy),f(max_momentum),f(max_angular),
        f(mass_drift),f(initial_bulk),f(center_shift),
        physical_gate.to_string(),f(elapsed),
    ];
    let _=writeln!(csv,"{}",line.join(","));
    println!("BG-PHY-23 hold N={n}, {kernel} {density} eta={eta}: radius bias={bias:+.3}%, contraction={contraction:+.3}%, Rmax drift={:.3}%, R99 drift={:.3}%, CDF={:.5}, E={max_energy:.2e}, P={max_momentum:.2e}, L={max_angular:.2e}, t_dyn={dyn_time:.4}s, gate={physical_gate}",
        100.0*max_rmax,100.0*max_r99,max_cdf);
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write BG-PHY-23 hold summary");
        println!("wrote {path}");
    }
    if let Some(path)=value(&args,"--trace-csv"){
        std::fs::write(path,trace).expect("write BG-PHY-23 hold trace");
        println!("wrote {path}");
    }
}
