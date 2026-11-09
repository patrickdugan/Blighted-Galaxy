//! Headless integration diagnostics for BG-PHY-04.
 //! cargo run --release --bin strongfield_bench
use blighted_galaxy::causal::{WorldlineHistory,WorldlineSample};
use blighted_galaxy::debris;
use blighted_galaxy::physics::{C_SI,SOLAR_MASS_KG};
use blighted_galaxy::schwarzschild::{
    gravitational_radius_m,SchwarzschildOrbit,OrbitStatus,periapsis_advance_weak_rad,
};
use blighted_galaxy::si_bridge::{Seconds,SiTidalEncounter};
use blighted_galaxy::V2;
use std::f64::consts::PI;

fn main() {
    let mass=10.0*SOLAR_MASS_KG;
    let rg=gravitational_radius_m(mass).unwrap();
    let mut orbit=SchwarzschildOrbit::bound_from_p_e(mass,30.0,0.2).unwrap();
    let mut apo=false;
    let mut precession=None;
    for _ in 0..90_000 {
        let old=orbit.radial_u;
        orbit.advance_proper(0.045*rg/C_SI).unwrap();
        if old>0.0&&orbit.radial_u<=0.0 {apo=true;}
        if apo&&old<0.0&&orbit.radial_u>=0.0 {
            precession=Some(orbit.phi_rad-2.0*PI);break;
        }
    }
    let precession=precession.expect("eccentric orbit must return to periapsis");
    println!("Schwarzschild 10 Msun p=30 rg e=0.2");
    println!("  observed periapsis advance {:.9} rad",precession);
    println!("  weak-field estimate         {:.9} rad",
        periapsis_advance_weak_rad(30.0).unwrap());
    println!("  radial energy constraint    {:.5e}",
        orbit.current_energy_constraint_error());
    assert!(precession>0.4&&precession<1.5);
    assert!(orbit.current_energy_constraint_error().abs()<1e-7);
    assert_eq!(orbit.status,OrbitStatus::Active);

    let mut history=WorldlineHistory::new(100);
    for i in 0..=30 {
        let t=i as f64/5.0;
        history.push(WorldlineSample {
            time_s:t,pos_m:V2::new(C_SI*(1.0+0.5*t),0.0),
            vel_m_s:V2::new(0.5*C_SI,0.0),
        }).unwrap();
    }
    let observation=history.observe(V2::default(),3.0).unwrap();
    println!("Retarded 0.5c emitter at t=3s");
    println!("  emission time {:.9}s, delay {:.9}s",
        observation.emitted_at_s,observation.delay_s);
    assert!((observation.emitted_at_s-4.0/3.0).abs()<1e-10);

    let mut tide=SiTidalEncounter::new();
    for _ in 0..2500 {
        tide.step(Seconds(0.1));
        if tide.disruption.is_some() {break;}
    }
    let event=tide.disruption.expect("tidal disruption required");
    let p0=debris::system_momentum_kg_m_s(&tide.black_hole,&tide.fragments);
    let mass0=debris::system_mass_kg(&tide.black_hole,&tide.fragments);
    for _ in 0..100 {tide.step(Seconds(0.02));}
    let p1=debris::system_momentum_kg_m_s(&tide.black_hole,&tide.fragments);
    let mass1=debris::system_mass_kg(&tide.black_hole,&tide.fragments);
    let scale=tide.fragments[0].mass_kg*tide.fragments[0].vel_m_s.len();
    let p_err=(p1-p0).len()/scale;
    let m_err=(mass1-mass0).abs()/mass0;
    println!("Self-gravitating 96-tracer breakup");
    println!("  trigger t={:.3}s q={:.5}",event.time_s,event.ratio);
    println!("  2-second post-breakup relative momentum drift {:.5e}",p_err);
    println!("  relative traced mass drift {:.5e}",m_err);
    assert!(p_err<1e-8 && m_err<1e-12);
}
