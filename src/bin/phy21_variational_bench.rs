//! BG-PHY-21B total variational pressure-gradient oracle (N16-32).
//!
//! Recompute implicit h/rho, MLS density and cold EOS energy for every
//! coordinate perturbation; forces = -finite_diff(U). Compare with fast
//! antisymmetric CRK pair pressure and test short pressure-only Verlet.
//! No direct changes to the gameplay dynamics and no gravity here.
//!
//! cargo run --release --bin phy21_variational_bench --
//!   --n 16 --steps 2 --dt 0.002 --csv bg-phy21b-n16.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::AdaptiveConfig;
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_rk_compatible::DensityChoice;
use blighted_galaxy::sph_rk_variational::{self,VariationalPilot};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
fn value<'a>(a:&'a [String],key:&str)->Option<&'a str>{
    a.windows(2).find(|w|w[0]==key).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(a:&[String],key:&str,default:T)->T{
    value(a,key).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn fmt(x:f64)->String{format!("{x:.10e}")}
fn total_momentum(p:&[Fluid3])->V3{
    p.iter().fold(V3::default(),|s,x|s+x.vel_m_s*x.mass_kg)
}
fn report_force_step_comparison(
    a:&sph_rk_variational::VariationalGradient,
    b:&sph_rk_variational::VariationalGradient,
)->f64 {
    let mut numerator=0.0;
    let mut denominator=0.0;
    for (fa,fb) in a.forces_n.iter().zip(&b.forces_n){
        numerator+=(*fa-*fb).len_sq();
        denominator+=fb.len_sq();
    }
    (numerator/denominator.max(1.0)).sqrt()
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",16);
    let seed:u64=opt(&args,"--seed",17);
    let steps:usize=opt(&args,"--steps",2);
    let dt:f64=opt(&args,"--dt",0.002);
    let frac:f64=opt(&args,"--fd-fraction",1e-4);
    let eta:f64=opt(&args,"--eta",1.7);
    let layout=value(&args,"--layout").unwrap_or("shell");
    assert!(n>=16&&n<=32&&n%2==0);
    assert!(steps>=1&&steps<=12);
    assert!(dt>0.0&&dt<=0.02&&dt.is_finite());
    assert!((1e-7..=1e-2).contains(&frac));
    assert!(eta>0.5&&eta<3.0);
    assert!(layout=="shell"||layout=="glass");
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let original=seed_cubic_shells_oriented(
        &star,center,V3::default(),n,seed).0;
    let p=if layout=="glass"{
        sph_glass::relax_antipodal_angles(&original,&star,
            center,1.2,8).0
    }else{original};
    let mut csv=String::from(
        "N,seed,layout,kernel,density,eta,fd_fraction,fd_step_m,gradient_fd_halving_rms,variational_net_force_ratio,variational_net_torque_ratio,rk_to_variational_force_rms,relative_hamiltonian_drift,relative_momentum_drift,reverse_position_error,reverse_velocity_error,steps,dt_s\n");
    for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4] {
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==DensityKernel::WendlandC4 {cfg=cfg.with_wendland_c4();}
        cfg.eta=eta;
        for density in [DensityChoice::Original,DensityChoice::SurfaceBlended] {
            let mut pilot=VariationalPilot{
                particles:p.clone(),physics,cfg,
                reference_radius_m:star.radius_m,mode:density,
                difference_fraction:frac,
            };
            let original=pilot.clone();
            let g0=pilot.gradient();
            let finer=sph_rk_variational::total_gradient(
                &p,physics,cfg,star.radius_m,density,frac*0.5);
            let halving=report_force_step_comparison(&g0,&finer);
            let e0=pilot.energy();
            let p0=total_momentum(&pilot.particles);
            let mass=pilot.particles.iter().map(|x|x.mass_kg).sum::<f64>();
            for _ in 0..steps {pilot.advance(dt);}
            let e1=pilot.energy();
            let p1=total_momentum(&pilot.particles);
            let drift=(e1-e0).abs()/e0.abs();
            let mom=(p1-p0).len()/(mass*1e6);
            for _ in 0..steps {pilot.advance(-dt);}
            let dx=pilot.particles.iter().zip(&original.particles)
                .map(|(a,b)|(a.pos_m-b.pos_m).len())
                .fold(0.0_f64,f64::max)/star.radius_m;
            let dv=pilot.particles.iter().zip(&original.particles)
                .map(|(a,b)|(a.vel_m_s-b.vel_m_s).len())
                .fold(0.0_f64,f64::max)/1e6;
            let kind=match density {
                DensityChoice::Original=>"original",
                DensityChoice::SurfaceBlended=>"surface_blended_mls",
            };
            let row=vec![
                n.to_string(),seed.to_string(),layout.into(),
                kernel.label().into(),kind.into(),fmt(eta),fmt(frac),
                fmt(g0.finite_difference_step_m),fmt(halving),
                fmt(g0.net_force_ratio),fmt(g0.net_torque_ratio),
                fmt(g0.rk_discrepancy_rms),fmt(drift),fmt(mom),
                fmt(dx),fmt(dv),steps.to_string(),fmt(dt),
            ];
            let _=writeln!(csv,"{}",row.join(","));
            println!("BG-PHY-21B N={n} {} {kind:18} eta={eta:.3} fd-halving={halving:.2e} netF={:.2e} netTorque={:.2e} RK-vs-var={:.3} energy={drift:.2e} momentum={mom:.2e} reverse=({dx:.2e},{dv:.2e})",
                kernel.label(),g0.net_force_ratio,g0.net_torque_ratio,
                g0.rk_discrepancy_rms);
            assert!(g0.net_force_ratio<5e-3);
            assert!(g0.net_torque_ratio<5e-3);
            assert!(halving<0.05);
            assert!(drift<2e-3);
            assert!(mom<1e-5);
        }
    }
    if let Some(path)=value(&args,"--csv") {
        std::fs::write(path,csv).expect("write BG-PHY-21B CSV");
        println!("wrote {path}");
    }
}
