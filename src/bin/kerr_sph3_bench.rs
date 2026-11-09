//! BG-PHY-05 reproducible Kerr/3D fluid/particle-batch regression benchmark.
 //! cargo run --release --bin kerr_sph3_bench
use blighted_galaxy::kerr::{KerrGeometry,KerrPhoton,PhotonStatus};
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::particle_batch::{
    build_meshes,count_indexed_quads,ParticleBillboard,
};
use blighted_galaxy::sph3::{
    self,Massive3,Sph3Config,V3,
};
use blighted_galaxy::physics::SiBody;
use blighted_galaxy::V2;
use macroquad::prelude::{Color,Vec2};

fn trace(bh:KerrGeometry,b:f64)->(PhotonStatus,f64,f64) {
    let mut ray=KerrPhoton::inbound(bh,46.0,b).unwrap();
    let mut max_residual=0.0_f64;
    for _ in 0..3_000 {
        ray.advance(0.055).unwrap();
        // The Boyer-Lindquist time/radial ODE becomes numerically stiff
        // at the horizon chart cutoff; assess convergence outside the
        // excluded near-horizon annulus, not at chart termination.
        if ray.status==PhotonStatus::Active &&
            ray.radius_rg>bh.outer_horizon_rg()+0.5 {
            max_residual=max_residual.max(ray.radial_constraint_residual().abs());
        }
        if ray.status!=PhotonStatus::Active {break;}
    }
    (ray.status,ray.phi_rad,max_residual)
}
fn main() {
    let mass=10.0*SOLAR_MASS_KG;
    let kerr=KerrGeometry::new(mass,0.85).unwrap();
    println!("BG-PHY-05: rotating Kerr + photon null geodesics");
    println!("a*=0.85: horizon={:.6} rg prograde_ISCO={:.6} retrograde_ISCO={:.6}",
        kerr.outer_horizon_rg(),kerr.isco_rg(1).unwrap(),
        kerr.isco_rg(-1).unwrap());
    println!("equatorial photon critical b: prograde={:.6} retrograde={:.6} rg",
        kerr.photon_critical_b_rg(1).unwrap(),
        kerr.photon_critical_b_rg(-1).unwrap());
    println!("ZAMO frame drag at 5 rg: {:.6e} rad/s",
        kerr.zamo_omega_rad_s(5.0).unwrap());
    for b in [-8.0,-6.5,6.5,8.0] {
        let (status,phi,residual)=trace(kerr,b);
        println!("  b={b:+.1} rg : {:?} phi={phi:.5} null_res={residual:.3e}",status);
        assert!(residual<1e-4);
    }
    let radial=KerrGeometry::new(mass,0.0).unwrap();
    let (status,_,_) = trace(radial,0.0);
    assert_eq!(status,PhotonStatus::HorizonChartCutoff);

    let star=SiBody::new(0.65*SOLAR_MASS_KG,9e6,
        V2::new(1.0e8,0.0),V2::new(-3e5,2e4));
    let (mut particles,_)=sph3::seed_from_progenitor(star,0.0,0.0,96);
    let mut central=Massive3{pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:10.0*SOLAR_MASS_KG,radius_m:3e4};
    let cfg=Sph3Config::white_dwarf_proxy();
    let e0=sph3::energy_j(&central,&particles,cfg);
    let p0=sph3::total_momentum(&central,&particles);
    let m0=sph3::total_mass_kg(&central,&particles);
    for _ in 0..160 {
        sph3::advance(&mut central,&mut particles,cfg,0.0025);
    }
    let e=sph3::energy_j(&central,&particles,cfg);
    let p=sph3::total_momentum(&central,&particles);
    let m=sph3::total_mass_kg(&central,&particles);
    let energy_drift=(e-e0).abs()/e0.abs();
    let momentum_drift=(p-p0).len()/(m0*3e5);
    let mass_drift=(m-m0).abs()/m0;
    println!("3-D polytropic fluid: {} parcels",particles.len());
    println!("  0.4-s energy drift {:.6e}",energy_drift);
    println!("  linear momentum drift {:.6e}",momentum_drift);
    println!("  mass drift {:.6e}",mass_drift);
    assert!(energy_drift<0.005);
    assert!(momentum_drift<1e-10);
    assert!(mass_drift<1e-13);

    let billboard=ParticleBillboard {
        center_px:Vec2::new(100.0,100.0),radius_px:1.0,
        tint:Color::new(0.75,0.9,1.0,1.0),depth:0.0
    };
    let meshes=build_meshes(&vec![billboard;10_000]);
    println!("GPU mesh batching: {} quads into {} draw calls",
        count_indexed_quads(&meshes),meshes.len());
    assert_eq!(count_indexed_quads(&meshes),10_000);
    assert_eq!(meshes.len(),3);
}
