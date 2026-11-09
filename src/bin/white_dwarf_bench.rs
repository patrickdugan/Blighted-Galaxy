//! Calibrated cold white-dwarf structure and 3D SPH acceptance benchmark.
 //! cargo run --release --bin white_dwarf_bench -- --output wd-mass-radius.csv
use blighted_galaxy::white_dwarf::{
    DegenerateElectrons,profile_for_mass,solve_hydrostatic,
};
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{
    self,Massive3,Sph3Config,V3,
};

fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let output=args.windows(2).find(|x|x[0]=="--output").map(|x|x[1].as_str());
    let eos=DegenerateElectrons::new(2.0).unwrap();
    println!("Cold fully relativistic Fermi-electron EOS: mu_e=2");
    println!("  pressure at rho=1e9 kg/m^3: {:.7e} Pa",eos.pressure_pa(1e9));
    println!("  non-relativistic K: {:.7e}",eos.nonrel_k_si());

    let mut csv=String::from("rho_c_kg_m3,mass_solar,radius_km\n");
    for density in [1e8,3e8,1e9,3e9,1e10,3e10,1e11] {
        let star=solve_hydrostatic(density,eos,10_000.0).unwrap();
        println!("  rho_c={density:.2e} kg/m^3 -> M={:.6} Msun, R={:.3} km",
            star.mass_solar(),star.radius_m/1e3);
        csv.push_str(&format!("{density:.7e},{:.9},{:.6}\n",
            star.mass_solar(),star.radius_m/1e3));
    }
    if let Some(p)=output {
        std::fs::write(p,csv).expect("mass-radius output");
        println!("  wrote {p}");
    }
    let equilibrium=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    println!("Target 0.65 Msun: rho_c={:.5e} kg/m^3, radius={:.3} km",
        equilibrium.rho_c_kg_m3,equilibrium.radius_m/1e3);
    assert!((equilibrium.mass_solar()-0.65).abs()<1e-5);
    assert!((equilibrium.radius_m/1e3-8402.2).abs()<100.0);

    let mut particles=sph3::seed_from_hydrostatic(
        &equilibrium,V3::new(1.0e8,0.0,0.0),
        V3::new(-3.0e5,2.0e4,0.0),128);
    let mut central=Massive3 {
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:10.0*SOLAR_MASS_KG,radius_m:3.0e4,
    };
    let cfg=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let p0=sph3::total_momentum(&central,&particles);
    let e0=sph3::energy_j(&central,&particles,cfg);
    let m0=sph3::total_mass_kg(&central,&particles);
    for _ in 0..80 {sph3::advance(&mut central,&mut particles,cfg,0.001);}
    let p1=sph3::total_momentum(&central,&particles);
    let e1=sph3::energy_j(&central,&particles,cfg);
    let m1=sph3::total_mass_kg(&central,&particles);
    let momentum_drift=(p1-p0).len()/(m0*3.0e5);
    let energy_drift=(e1-e0).abs()/e0.abs();
    let mass_drift=(m1-m0).abs()/m0;
    println!("128-parcel 3D cold electron EOS (0.08 seconds)");
    println!("  relative momentum drift {momentum_drift:.6e}");
    println!("  relative energy drift   {energy_drift:.6e}");
    println!("  relative mass drift     {mass_drift:.6e}");
    assert!(momentum_drift<1e-9 && mass_drift<1e-13);
    assert!(energy_drift.is_finite());
}
