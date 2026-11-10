//! BG-PHY-22 analytic complete cold-barotropic density/MLS variation.
//!
//! The fast analytic reverse accumulation differentiates SPH rho, the
//! implicit h root, source volumes, S, affine MLS and smoothstep blend.
//! Compare each analytic force with BG-PHY-21B's independent full
//! finite-difference energy-gradient oracle. Cold EOS is unchanged.
//!
//! cargo run --release --bin phy22_bench -- --n 16 --steps 3 --dt 0.002
//!   --csv bg-phy22-n16.csv
//! cargo run --release --bin phy22_bench -- --n 64 --no-oracle --steps 4
//!   --csv bg-phy22-n64.csv
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_rk_compatible::{self,DensityChoice};
use blighted_galaxy::sph_rk_variational;
use blighted_galaxy::sph_surface_variation::{self,AnalyticPilot};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
use std::time::Instant;
fn value<'a>(args:&'a [String],key:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==key).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T{
    value(args,key).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn f(x:f64)->String{format!("{x:.10e}")}
fn rms_diff(a:&[V3],b:&[V3])->f64{
    assert_eq!(a.len(),b.len());
    let (mut num,mut denom)=(0.0,0.0);
    for (x,y) in a.iter().zip(b) {
        num+=(*x-*y).len_sq();denom+=y.len_sq();
    }
    (num/denom.max(1e-200)).sqrt()
}
fn momentum(p:&[Fluid3])->V3{
    p.iter().fold(V3::default(),|sum,x|sum+x.vel_m_s*x.mass_kg)
}
fn cross(a:V3,b:V3)->V3{
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn angular(p:&[Fluid3])->V3{
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let c=p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/m;
    let v=momentum(p)*(1.0/m);
    p.iter().fold(V3::default(),|s,x|
        s+cross(x.pos_m-c,(x.vel_m_s-v)*x.mass_kg))
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",16);
    let seed:u64=opt(&args,"--seed",17);
    let layout=value(&args,"--layout").unwrap_or("shell");
    let steps:usize=opt(&args,"--steps",3);
    let dt:f64=opt(&args,"--dt",0.002);
    let eta:f64=opt(&args,"--eta",1.7);
    let step_fraction:f64=opt(&args,"--fd-fraction",1e-4);
    let oracle=!args.iter().any(|x|x=="--no-oracle");
    assert!((16..=128).contains(&n)&&n%2==0);
    assert!(layout=="shell"||layout=="glass");
    assert!(steps>0&&steps<=30&&dt.is_finite()&&dt>0.0&&dt<=0.02);
    assert!(eta>0.5&&eta<3.0);
    assert!(!oracle||n<=64);
    assert!((1e-7..=1e-2).contains(&step_fraction));
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let initial=seed_cubic_shells_oriented(
        &star,center,V3::default(),n,seed).0;
    let original=if layout=="glass"{
        sph_glass::relax_antipodal_angles(&initial,&star,
            center,1.2,8).0
    }else{initial};
    let mut csv=String::from(
        "N,seed,layout,kernel,density,eta,analytic_s,finite_difference_s,analytic_vs_oracle_rms,conventional_vs_analytic_rms,rk_vs_analytic_rms,potential_relative_error,net_force_ratio,net_torque_ratio,maximum_mls_scalar_disagreement,mls_fallback_fraction,capped_fraction,omega_floor_fraction,min_omega,direct_displacement_rms_n,density_chain_rms_n,implicit_h_chain_rms_n,relative_hamiltonian_drift,linear_momentum_drift,angular_momentum_drift,reverse_position_error,reverse_velocity_error,steps,dt_s\n");
    for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4] {
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==DensityKernel::WendlandC4{cfg=cfg.with_wendland_c4();}
        cfg.eta=eta;
        for density in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
            let started=Instant::now();
            let analytic=sph_surface_variation::total_energy_gradient(
                &original,physics,cfg,star.radius_m,density);
            let analytic_s=started.elapsed().as_secs_f64();
            let (fd_s,fd_error,potential_error)=if oracle{
                let started=Instant::now();
                let fd=sph_rk_variational::total_gradient(
                    &original,physics,cfg,star.radius_m,
                    density,step_fraction);
                (started.elapsed().as_secs_f64(),
                    rms_diff(&analytic.forces_n,&fd.forces_n),
                    (analytic.potential_j-fd.potential_j).abs()/
                        fd.potential_j.abs().max(1.0))
            }else{(0.0,-1.0,0.0)};
            let empty=Massive3{pos_m:V3::default(),
                vel_m_s:V3::default(),mass_kg:0.0,radius_m:0.0};
            let mut state=AdaptiveState::new(
                &original,star.radius_m,cfg);
            let conventional=sph_adaptive::forces(
                &empty,&original,&mut state,physics,cfg);
            let mut rk_state=AdaptiveState::new(
                &original,star.radius_m,cfg);
            let rk=sph_rk_compatible::pressure_snapshot(
                &original,&mut rk_state,physics,cfg,density);
            let old_force=conventional.pressure_accelerations_m_s2.iter()
                .zip(&original).map(|(a,x)|*a*x.mass_kg).collect::<Vec<_>>();
            let rk_force=rk.accelerations_m_s2.iter().zip(&original)
                .map(|(a,x)|*a*x.mass_kg).collect::<Vec<_>>();
            let conventional_error=rms_diff(&old_force,&analytic.forces_n);
            let rk_error=rms_diff(&rk_force,&analytic.forces_n);
            let mut pilot=AnalyticPilot{particles:original.clone(),
                physics,cfg,reference_radius_m:star.radius_m,density};
            let energy0=pilot.energy_j();
            let p0=momentum(&pilot.particles);
            let l0=angular(&pilot.particles);
            for _ in 0..steps{pilot.advance(dt);}
            let energy1=pilot.energy_j();
            let mass=pilot.particles.iter().map(|x|x.mass_kg).sum::<f64>();
            let linear_drift=(momentum(&pilot.particles)-p0).len()/(mass*1e6);
            let angular_drift=(angular(&pilot.particles)-l0).len()/
                pilot.particles.iter().map(|x|
                    x.mass_kg*(x.pos_m-center).len()*1e6).sum::<f64>();
            let drift=(energy1-energy0).abs()/energy0.abs();
            for _ in 0..steps{pilot.advance(-dt);}
            let dx=pilot.particles.iter().zip(&original).map(|(a,b)|
                (a.pos_m-b.pos_m).len()).fold(0.0_f64,f64::max)/star.radius_m;
            let dv=pilot.particles.iter().zip(&original).map(|(a,b)|
                (a.vel_m_s-b.vel_m_s).len()).fold(0.0_f64,f64::max)/1e6;
            let kind=match density{
                DensityChoice::Original=>"original",
                DensityChoice::SurfaceBlended=>"surface_blended_mls",
            };
            let cols=vec![
                n.to_string(),seed.to_string(),layout.into(),
                kernel.label().into(),kind.into(),f(eta),
                f(analytic_s),f(fd_s),f(fd_error),
                f(conventional_error),f(rk_error),f(potential_error),
                f(analytic.net_force_ratio),f(analytic.net_torque_ratio),
                f(analytic.max_mls_relative_disagreement),
                f(analytic.mls_fallback_fraction),f(analytic.capped_fraction),
                f(analytic.omega_floor_fraction),f(analytic.minimum_omega),
                f(analytic.direct_force_rms_n),f(analytic.density_force_rms_n),
                f(analytic.smoothing_force_rms_n),
                f(drift),f(linear_drift),f(angular_drift),
                f(dx),f(dv),steps.to_string(),f(dt),
            ];
            let _=writeln!(csv,"{}",cols.join(","));
            println!("BG-PHY-22 N={n} {} {kind:18} FD={fd_error:.3e} old-vs-var={conventional_error:.3e} RK-vs-var={rk_error:.3e} netF={:.2e} netT={:.2e} drift={drift:.2e} P={linear_drift:.2e} L={angular_drift:.2e}, reverse={dx:.2e}/{dv:.2e}",
                kernel.label(),analytic.net_force_ratio,
                analytic.net_torque_ratio);
            if oracle {
                assert!(fd_error<5e-5,
                    "analytic total chain must match independent FD oracle: {fd_error}");
            }
            assert!(analytic.net_force_ratio<1e-9);
            assert!(analytic.net_torque_ratio<1e-9);
            assert!(analytic.max_mls_relative_disagreement<1e-9);
            if n==64 {assert!(analytic.capped_fraction<1e-12,
                "N64 benchmark should exercise uncapped implicit h derivatives");}
            assert!(drift<1e-5);
            assert!(linear_drift<1e-9&&angular_drift<1e-9);
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write analytic-variation CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn cross_product_rotation_identity(){
        assert_eq!(cross(V3::new(1.0,0.0,0.0),V3::new(0.0,1.0,0.0)),
            V3::new(0.0,0.0,1.0));
    }
}
