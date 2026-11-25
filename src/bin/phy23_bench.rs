//! BG-PHY-23: compare complete cold EOS + matched Poisson gravity force
//! with independent full scalar-energy coordinate finite differences.
//!
//! The pressure component is BG-PHY-22's analytic surface variation.
//! Gravity uses the existing Poisson kernel potential and zeta/Omega
//! grad-h terms. Same variable h and baseline SPH density in both.
//! No gameplay solver is changed.
//!
//! cargo run --release --bin phy23_bench -- --n 24 --steps 3
//!   --csv bg-phy23-n24.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::AdaptiveConfig;
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_rk_compatible::DensityChoice;
use blighted_galaxy::sph_surface_gravity::{self,HamiltonianPilot};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
use std::time::Instant;
fn arg<'a>(a:&'a [String],key:&str)->Option<&'a str> {
    a.windows(2).find(|x|x[0]==key).map(|x|x[1].as_str())
}
fn opt<T:std::str::FromStr>(a:&[String],key:&str,d:T)->T{
    arg(a,key).and_then(|x|x.parse().ok()).unwrap_or(d)
}
fn f(x:f64)->String {format!("{x:.10e}")}
fn momentum(p:&[Fluid3])->V3 {
    p.iter().fold(V3::default(),|s,x|s+x.vel_m_s*x.mass_kg)
}
fn cross(a:V3,b:V3)->V3 {
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn angular(p:&[Fluid3])->V3 {
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let c=p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/m;
    let v=momentum(p)*(1.0/m);
    p.iter().fold(V3::default(),|s,x|
        s+cross(x.pos_m-c,(x.vel_m_s-v)*x.mass_kg))
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",24);
    let seed:u64=opt(&args,"--seed",17);
    let layout=arg(&args,"--layout").unwrap_or("glass");
    let eta:f64=opt(&args,"--eta",1.7);
    let steps:usize=opt(&args,"--steps",4);
    let dt:f64=opt(&args,"--dt",0.002);
    let fd_fraction:f64=opt(&args,"--fd-fraction",1e-4);
    let no_oracle=args.iter().any(|x|x=="--no-oracle");
    let deform=args.iter().any(|x|x=="--deform");
    assert!((16..=128).contains(&n)&&n%2==0);
    assert!(layout=="shell"||layout=="glass");
    assert!(no_oracle||n<=64);
    assert!(dt>0.0&&dt<=0.02&&steps>0&&steps<=16);
    assert!(eta>0.5&&eta<3.0);
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let initial=seed_cubic_shells_oriented(&star,
        center,V3::default(),n,seed).0;
    let mut p=if layout=="glass" {
        sph_glass::relax_antipodal_angles(
            &initial,&star,center,1.2,8).0
    }else{initial};
    if deform{
        for (i,x) in p.iter_mut().enumerate(){
            let t=i as f64+0.3;
            x.pos_m+=V3::new((0.73*t).sin(),
                (1.19*t).cos(),(1.53*t).sin())*
                (0.013*star.radius_m);
            x.mass_kg*=1.0+0.055*(0.41*t).sin();
        }
    }
    let mut csv=String::from(
        "N,seed,layout,deformed,kernel,density,eta,fd_fraction,analytic_ms,oracle_ms,force_rms_error,scalar_energy_relative_error,old_sph_force_error,net_force_ratio,net_torque_ratio,pressure_energy_j,matched_gravity_j,minimum_omega,capped_fraction,MLS_fallback_fraction,relative_total_energy_drift,momentum_drift,angular_momentum_drift,reversed_position_error,reversed_velocity_error,steps,dt_s\n");
    for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4] {
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==DensityKernel::WendlandC4 {
            cfg=cfg.with_wendland_c4();
        }
        cfg.eta=eta;
        for density in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
            let started=Instant::now();
            let analytic=sph_surface_gravity::snapshot(
                &p,physics,cfg,star.radius_m,density);
            let analytic_ms=1000.0*started.elapsed().as_secs_f64();
            let (oracle_ms,error,scalar_error)=if no_oracle {
                (0.0,-1.0,0.0)
            }else{
                let t=Instant::now();
                let oracle=sph_surface_gravity::independent_total_gradient(
                    &p,physics,cfg,star.radius_m,density,fd_fraction);
                (1000.0*t.elapsed().as_secs_f64(),
                 oracle.relative_force_difference,
                 oracle.relative_potential_difference)
            };
            let mut pilot=HamiltonianPilot{
                particles:p.clone(),physics,cfg,
                reference_radius_m:star.radius_m,density,
            };
            let before=pilot.clone();
            let energy0=pilot.energy_j();
            let momentum0=momentum(&pilot.particles);
            let angular0=angular(&pilot.particles);
            for _ in 0..steps{pilot.advance(dt);}
            let energy1=pilot.energy_j();
            let mass=pilot.particles.iter()
                .map(|x|x.mass_kg).sum::<f64>();
            let momentum_drift=(momentum(&pilot.particles)-momentum0)
                .len()/(mass*1e6);
            let angular_drift=(angular(&pilot.particles)-angular0)
                .len()/(mass*star.radius_m*1e6);
            let drift=(energy1-energy0).abs()/
                energy0.abs().max(1.0);
            for _ in 0..steps{pilot.advance(-dt);}
            let dx=pilot.particles.iter().zip(&before.particles)
                .map(|(a,b)|(a.pos_m-b.pos_m).len())
                .fold(0.0_f64,f64::max)/star.radius_m;
            let dv=pilot.particles.iter().zip(&before.particles)
                .map(|(a,b)|(a.vel_m_s-b.vel_m_s).len())
                .fold(0.0_f64,f64::max)/1e6;
            assert!(no_oracle||error<2e-5,
                "total Hamiltonian analytic/FD error {error}");
            assert!(analytic.net_force_ratio<1e-9);
            assert!(analytic.net_torque_ratio<1e-9);
            assert!(drift<2e-5,"energy drift {drift}");
            assert!(momentum_drift<1e-9&&angular_drift<1e-9);
            assert!(dx<1e-7&&dv<1e-7);
            if density==DensityChoice::Original {
                assert!(analytic.old_sph_relative_force_error<1e-10,
                    "original SPH full-force control {}",
                    analytic.old_sph_relative_force_error);
            }
            let label=match density {
                DensityChoice::Original=>"original",
                DensityChoice::SurfaceBlended=>"surface_blended_mls",
            };
            let row=vec![
                n.to_string(),seed.to_string(),layout.into(),deform.to_string(),
                kernel.label().into(),label.into(),f(eta),f(fd_fraction),
                f(analytic_ms),f(oracle_ms),f(error),f(scalar_error),
                f(analytic.old_sph_relative_force_error),
                f(analytic.net_force_ratio),f(analytic.net_torque_ratio),
                f(analytic.cold_internal_energy_j),f(analytic.gravity_potential_j),
                f(analytic.minimum_omega),f(analytic.capped_fraction),
                f(analytic.pressure.mls_fallback_fraction),
                f(drift),f(momentum_drift),f(angular_drift),
                f(dx),f(dv),steps.to_string(),f(dt),
            ];
            let _=writeln!(csv,"{}",row.join(","));
            println!("BG-PHY-23 N={n} {} {label:18} FD={error:.3e} scalar={scalar_error:.2e} old-SPH={:.3e} caps={:.1}% E={drift:.2e} P={momentum_drift:.2e} L={angular_drift:.2e} reverse={dx:.2e}/{dv:.2e}",
                kernel.label(),analytic.old_sph_relative_force_error,
                100.0*analytic.capped_fraction);
        }
    }
    if let Some(path)=arg(&args,"--csv"){
        std::fs::write(path,csv).expect("write BG-PHY-23 CSV");
        println!("wrote {path}");
    }
}
