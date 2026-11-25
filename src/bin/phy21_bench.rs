//! BG-PHY-21: small-N pressure-only reproducing-kernel pilot.
//!
//! Reads original M4/C4 cold EOS, validates RK identities, compares frozen
//! pressure errors and runs reversible pressure-only compatible-work steps.
//! Momentum and K+U conservation are pass/fail; angular momentum and the
//! mismatch between evolved u and cold u_EOS(rho) are measured explicitly.
//!
//! No gravity, viscosity, star acceptance, or gameplay changes are implied.
//!
//! cargo run --release --bin phy21_bench -- --n 32 --steps 12 --dt 0.001
//!   --csv bg-phy21-small.csv
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{Fluid3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_kernel::DensityKernel;
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_glass;
use blighted_galaxy::sph_neighbor_calibration;
use blighted_galaxy::sph_force_audit;
use blighted_galaxy::sph_rk_compatible::{RkPilotState,DensityChoice,RkPressureSnapshot};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass,WhiteDwarfProfile};
use std::fmt::Write as _;
fn value<'a>(args:&'a [String],key:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==key).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],key:&str,default:T)->T{
    value(args,key).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn f(v:f64)->String{format!("{v:.10e}")}
fn cross(a:V3,b:V3)->V3{
    V3::new(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
fn pressure_error(
    p:&[Fluid3],star:&WhiteDwarfProfile,
    snap:&RkPressureSnapshot,
)->(f64,f64,f64){
    let mass=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let center=p.iter().fold(V3::default(),
        |s,x|s+x.pos_m*x.mass_kg)/mass;
    let mut index=(0..p.len()).collect::<Vec<_>>();
    index.sort_by(|&i,&j|
        (p[i].pos_m-center).len_sq()
            .total_cmp(&(p[j].pos_m-center).len_sq()).then(i.cmp(&j)));
    let mut err=0.0;
    let mut ref2=0.0;
    let mut err_outer=0.0;
    let mut ref2_outer=0.0;
    let mut outer_g=0.0;
    let mut outer_acc=0.0;
    let outer_start=p.len()-p.len()/4;
    for (k,i) in index.into_iter().enumerate(){
        let rvec=p[i].pos_m-center;
        let r=rvec.len();
        let g=G_SI*star.enclosed_mass_at(r)/(r*r);
        let u=rvec*(1.0/r);
        let a=snap.accelerations_m_s2[i];
        let m=p[i].mass_kg;
        err+=m*(a-u*g).len_sq();
        ref2+=m*g*g;
        if k>=outer_start{
            err_outer+=m*(a-u*g).len_sq();
            ref2_outer+=m*g*g;
            outer_g+=m*g;
            outer_acc+=m*a.dot(u);
        }
    }
    ((err/ref2).sqrt(),(err_outer/ref2_outer).sqrt(),
        outer_acc/outer_g)
}
fn torque_ratio(p:&[Fluid3],snap:&RkPressureSnapshot)->f64{
    let m=p.iter().map(|x|x.mass_kg).sum::<f64>();
    let c=p.iter().fold(V3::default(),|s,x|s+x.pos_m*x.mass_kg)/m;
    let mut torque=V3::default();
    let mut magnitude=0.0;
    for (i,x) in p.iter().enumerate(){
        let lever=x.pos_m-c;
        let force=snap.accelerations_m_s2[i]*x.mass_kg;
        torque+=cross(lever,force);
        magnitude+=lever.len()*force.len();
    }
    torque.len()/magnitude.max(1.0)
}
fn original_sp_err(
    p:&[Fluid3],physics:Sph3Config,
    cfg:AdaptiveConfig,star:&WhiteDwarfProfile
)->f64{
    let mut state=AdaptiveState::new(p,star.radius_m,cfg);
    sph_force_audit::frozen_force_audit(p,&mut state,physics,cfg,star,8)
        .global.pressure_ref_rms
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",32);
    let seed:u64=opt(&args,"--seed",17);
    let steps:usize=opt(&args,"--steps",16);
    let dt:f64=opt(&args,"--dt",0.001);
    let eta:f64=opt(&args,"--eta",1.7);
    let layout=value(&args,"--layout").unwrap_or("glass");
    assert!(n>=16&&n<=128&&n%2==0);
    assert!(steps>=1&&steps<=500);
    assert!(dt>0.0&&dt<=0.02&&dt.is_finite());
    assert!(eta>0.5&&eta<3.0);
    assert!(layout=="glass"||layout=="shell");
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,
        DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let center=V3::new(1e8,0.0,0.0);
    let source=seed_cubic_shells_oriented(&star,
        center,V3::default(),n,seed).0;
    let p=if layout=="glass" {
        sph_glass::relax_antipodal_angles(&source,&star,
            center,1.2,12).0
    }else{source};
    let mut csv=String::from(
        "N,seed,layout,kernel,density,eta,realized_neighbors,valid_rk_fraction,max_constant_error,max_linear_error,max_gradient_sum_error,legacy_pressure_ref_rms,rk_pressure_ref_rms,rk_outer_pressure_ref_rms,rk_outer_pressure_ratio,initial_pair_torque_ratio,energy_drift,linear_momentum_drift,angular_momentum_drift,eos_energy_mismatch,reverse_position_error,reverse_velocity_error,steps,dt_s\n");
    for kernel in [DensityKernel::CubicM4,DensityKernel::WendlandC4]{
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity().with_newton_density();
        if kernel==DensityKernel::WendlandC4 {
            cfg=cfg.with_wendland_c4();
        }
        cfg.eta=eta;
        let original_error=original_sp_err(&p,physics,cfg,&star);
        for density in [DensityChoice::Original,DensityChoice::SurfaceBlended]{
            let mut state=RkPilotState::new(
                p.clone(),physics,cfg,star.radius_m,density);
            let baseline=state.clone();
            let snap=state.forces();
            let reference=pressure_error(&state.particles,&star,&snap);
            let torque=torque_ratio(&state.particles,&snap);
            let support=sph_neighbor_calibration::measure(
                &state.particles,&state.sph_state,cfg).mean;
            let energy0=state.total_energy_j();
            let momentum0=state.momentum();
            let l0=state.angular_momentum();
            for _ in 0..steps{state.advance(dt);}
            let drift=(state.total_energy_j()-energy0).abs()/
                energy0.abs().max(1.0);
            let total_mass=state.particles.iter().map(|x|x.mass_kg).sum::<f64>();
            let dp=(state.momentum()-momentum0).len()/(total_mass*1e6);
            let l1=state.angular_momentum();
            let dl=(l1-l0).len()/
                state.particles.iter().map(|x|x.mass_kg*
                    (x.pos_m-center).len()*1e6).sum::<f64>().max(1.0);
            let eos=state.eos_consistency_rms();
            for _ in 0..steps{state.advance(-dt);}
            let dx=state.particles.iter().zip(&baseline.particles)
                .map(|(a,b)|(a.pos_m-b.pos_m).len())
                .fold(0.0_f64,f64::max)/star.radius_m;
            let dv=state.particles.iter().zip(&baseline.particles)
                .map(|(a,b)|(a.vel_m_s-b.vel_m_s).len())
                .fold(0.0_f64,f64::max)/1e6;
            assert!(drift<1e-9,"total compatible-work energy drift: {drift}");
            assert!(dp<1e-10,"momentum drift: {dp}");
            assert!(dx<1e-7&&dv<1e-7,
                "reverse Verlet mismatch: position={dx} velocity={dv}");
            let kind=match density{
                DensityChoice::Original=>"original",
                DensityChoice::SurfaceBlended=>"surface_blended_mls",
            };
            let cols=vec![
                n.to_string(),seed.to_string(),layout.into(),
                kernel.label().into(),kind.into(),f(eta),f(support),
                f(snap.valid_rk_fraction),
                f(snap.rk.maximum_constant_error),
                f(snap.rk.maximum_linear_error),
                f(snap.rk.maximum_query_gradient_error),
                f(original_error),f(reference.0),f(reference.1),
                f(reference.2),f(torque),f(drift),f(dp),f(dl),
                f(eos),f(dx),f(dv),steps.to_string(),f(dt),
            ];
            let _=writeln!(csv,"{}",cols.join(","));
            println!("BG-PHY-21 N={n} {} {kind:18} eta={eta:.3} ngb={support:.1} valid={:.0}% W0err={:.2e} W1err={:.2e} P-ref={:.4} outerP={:.4} torque={torque:.4} energy={drift:.2e} dp={dp:.2e} EOS={eos:.3e} reverse={dx:.2e}/{dv:.2e}",
                kernel.label(),100.0*snap.valid_rk_fraction,
                snap.rk.maximum_constant_error,
                snap.rk.maximum_linear_error,
                reference.0,reference.2);
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write BG-PHY-21 CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn signed_torque_is_finite(){
        let a=V3::new(1.0,0.0,0.0);
        let b=V3::new(0.0,1.0,0.0);
        assert_eq!(cross(a,b),V3::new(0.0,0.0,1.0));
    }
}
