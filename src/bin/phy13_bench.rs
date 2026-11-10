//! BG-PHY-13: cold-white-dwarf smoothing length and radial gravity sweep.
//!
//! Run: cargo run --release --bin phy13_bench -- --n 256 \
//!   --etas 0.85,1.05,1.20,1.45 --relax-steps 450 --hold-steps 788 \
//!   --csv bg-phy13.csv --shell-csv bg-phy13-shells.csv
//!
//! All eta variants use the same near-cubic mass-CDF particle seed, 0.65-Msun
//! mu_e=2 cold-electron EOS, M4 Poisson self-gravity, dt, damping duration,
//! and subsequent undamped hold. Reference gravity is computed independently
//! from the 1D hydrostatic mass profile at each observed shell radius.
//! The physical 5% acceptance gate additionally demands >=3 t_dyn of hold.
//! A short smoke test MUST NOT be called a successful equilibrium.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Fluid3,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{self,AdaptiveConfig,AdaptiveState};
use blighted_galaxy::sph_sampling::{self,StarSampling};
use blighted_galaxy::sph_massprofile::{self,MassProfileDiagnostic};
use blighted_galaxy::sph_neighbors::CellNeighbors;
use blighted_galaxy::stellar_balance::{self,StellarBalance};
use blighted_galaxy::stellar_relax::{relax_isolated_adaptive,RelaxParameters};
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::time::Instant;
use std::fmt::Write as _;

fn option<'a>(a:&'a [String],k:&str)->Option<&'a str> {
    a.windows(2).find(|p|p[0]==k).map(|p|p[1].as_str())
}
fn parse<T:std::str::FromStr>(a:&[String],k:&str,default:T)->T {
    option(a,k).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn eta_list(args:&[String])->Vec<f64> {
    let raw=option(args,"--etas").unwrap_or("0.85,1.05,1.2,1.45");
    let values=raw.split(',').map(|s|s.parse::<f64>()
        .expect("--etas must be comma-separated finite numbers"))
        .collect::<Vec<_>>();
    assert!(!values.is_empty()&&values.len()<=12);
    assert!(values.iter().all(|x|x.is_finite()&&*x>0.5&&*x<3.0));
    for w in values.windows(2) {
        assert!(w[1]>w[0],"--etas must be strictly increasing");
    }
    values
}
fn gravity_reference_rms(balance:&StellarBalance)->f64 {
    let (numerator,denominator)=balance.shells.iter().fold(
        (0.0,0.0),|(a,b),s|{
            let weight=s.mass_fraction_hi-s.mass_fraction_lo;
            (a+weight*(s.gravity_in_m_s2-
                s.hydrostatic_gravity_in_m_s2).powi(2),
             b+weight*s.hydrostatic_gravity_in_m_s2.powi(2))
        });
    (numerator/denominator).sqrt()
}
fn mean_candidates(particles:&[Fluid3],state:&AdaptiveState)->f64 {
    let h=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    CellNeighbors::new(particles,h).mean_candidate_count()
}
struct Case {
    eta:f64,
    rmax_post_km:f64,
    radius_bias_pct:f64,
    r99_post_km:f64,
    mass_cdf_post:f64,
    gravity_reference_rms_post:f64,
    force_ratio_post:f64,
    virial_post:f64,
    inner_gravity_reference_ratio:f64,
    outer_gravity_reference_ratio:f64,
    omega_min:f64,
    cap_fraction:f64,
    mean_neighbors:f64,
    largest_outer_drift_pct:f64,
    largest_r99_drift_pct:f64,
    max_cdf_rms:f64,
    final_energy_drift:f64,
    final_momentum_drift:f64,
    elapsed_wall_s:f64,
    completed_three_dyn:bool,
    physical_gate_pass:bool,
}
fn record_shells(
    out:&mut String,n:usize,eta:f64,phase:&str,bal:&StellarBalance,
) {
    for (i,s) in bal.shells.iter().enumerate() {
        let g_ratio=s.gravity_in_m_s2/
            s.hydrostatic_gravity_in_m_s2.max(1e-20);
        let _=writeln!(out,
            "{n},{eta:.5},{phase},{i},{:.6},{:.6},{},{:.4},{:.8e},{:.8e},{:.8e},{g_ratio:.7},{:.8e},{:.8e},{:.8e}",
            s.mass_fraction_lo,s.mass_fraction_hi,s.particle_count,
            s.mean_radius_m/1000.0,s.pressure_out_m_s2,
            s.gravity_in_m_s2,s.hydrostatic_gravity_in_m_s2,
            s.net_radial_m_s2,s.net_radial_rms_m_s2,
            s.tangential_net_rms_m_s2,
        );
    }
}
fn main() {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=parse(&args,"--n",128);
    let dt:f64=parse(&args,"--dt",0.01);
    let relax_steps:usize=parse(&args,"--relax-steps",400);
    let hold_steps:usize=parse(&args,"--hold-steps",788);
    let tau_s:f64=parse(&args,"--tau",0.7);
    let etas=eta_list(&args);
    assert!(n>=32&&n<=4096&&n%2==0);
    assert!(relax_steps>0&&relax_steps<=30_000);
    assert!(hold_steps>0&&hold_steps<=30_000);
    assert!(dt.is_finite()&&dt>0.0&&dt<=0.1);
    assert!(tau_s.is_finite()&&tau_s>0.0);
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let seed=sph_sampling::seeded_star(&star,
        V3::new(1.0e8,0.0,0.0),V3::default(),n,
        StarSampling::ApproximatelyCubicShells);
    let three_dyn=hold_steps as f64*dt>=3.0*t_dyn;
    println!("BG-PHY-13 ETA / RESOLUTION / LOCAL HYDROSTATIC GRAVITY SWEEP");
    println!("  N={n} ref_R={:.3} km t_dyn={t_dyn:.6} s, relax={:.3} s hold={:.3} s ({:.3} t_dyn)",
        star.radius_m/1000.0,relax_steps as f64*dt,
        hold_steps as f64*dt,hold_steps as f64*dt/t_dyn);
    if !three_dyn {
        println!("  HOLD LENGTH INSUFFICIENT: numerical smoke only, equilibrium acceptance disabled");
    }
    let mut csv=String::from(
        "N,eta,relax_s,hold_s,hydro_radius_km,post_rmax_km,post_radius_bias_pct,post_r99_km,post_cdf_rms,post_gravity_reference_rms,post_net_gravity_rms,post_force_virial,inner_gravity_reference_ratio,outer_gravity_reference_ratio,min_omega,cap_fraction,mean_neighbors,worst_rmax_drift_pct,worst_r99_drift_pct,max_cdf_rms,final_energy_drift,final_momentum_drift,wall_seconds,complete_three_dyn,physical_gate_pass\n"
    );
    let mut shell_csv=String::from(
        "N,eta,phase,shell,mass_lo,mass_hi,n,mean_radius_km,pressure_out_m_s2,gravity_in_m_s2,hydro_gravity_in_m_s2,gravity_reference_ratio,net_radial_m_s2,net_radial_rms_m_s2,tangential_rms_m_s2\n"
    );
    let mut outcomes=Vec::<Case>::new();
    for eta in etas {
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity();
        cfg.eta=eta;
        assert!(cfg.valid());
        let mut particles=seed.clone();
        let mut state=AdaptiveState::new(&particles,star.radius_m,cfg);
        let start=Instant::now();
        let report=relax_isolated_adaptive(
            &mut particles,physics,&star,cfg,&mut state,
            RelaxParameters{iterations:relax_steps,dt_s:dt,
                damping_time_s:tau_s}
        ).expect("cold degenerate relaxation must remain finite");
        let profile=sph_massprofile::profile_diagnostics(&particles,&star);
        let balance=stellar_balance::radial_balance(
            &particles,&mut state,physics,cfg,&star,8);
        let gravity_rms=gravity_reference_rms(&balance);
        let neighbors=mean_candidates(&particles,&state);
        let inner=balance.shells[0].gravity_in_m_s2/
            balance.shells[0].hydrostatic_gravity_in_m_s2;
        let last=balance.shells.last().unwrap();
        let outer=last.gravity_in_m_s2/last.hydrostatic_gravity_in_m_s2;
        record_shells(&mut shell_csv,n,eta,"post_damp",&balance);
        let mut central=Massive3{
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:0.0,radius_m:0.0,
        };
        let m0=sph3::total_mass_kg(&central,&particles);
        let p0=sph3::total_momentum(&central,&particles);
        let e0=sph_adaptive::adaptive_energy_j(
            &central,&particles,&mut state,physics,cfg);
        let mut outer_drift=0.0_f64;
        let mut r99_drift=0.0_f64;
        let mut max_cdf=profile.radial_quantile_rms;
        let sample_interval=(hold_steps/32).max(1);
        let mut final_energy=0.0;
        let mut final_p=0.0;
        for step in 1..=hold_steps {
            sph_adaptive::advance(
                &mut central,&mut particles,&mut state,physics,cfg,dt);
            if step%sample_interval==0||step==hold_steps {
                let q=sph_massprofile::profile_diagnostics(&particles,&star);
                outer_drift=outer_drift.max(
                    (q.radius_max_m/profile.radius_max_m-1.0).abs());
                r99_drift=r99_drift.max(
                    (q.r99_m/profile.r99_m-1.0).abs());
                max_cdf=max_cdf.max(q.radial_quantile_rms);
                if step==hold_steps {
                    let f=stellar_balance::radial_balance(
                        &particles,&mut state,physics,cfg,&star,8);
                    record_shells(&mut shell_csv,n,eta,"final_hold",&f);
                }
                let energy=sph_adaptive::adaptive_energy_j(
                    &central,&particles,&mut state,physics,cfg);
                let momentum=sph3::total_momentum(&central,&particles);
                final_energy=(energy-e0).abs()/e0.abs().max(1.0);
                final_p=(momentum-p0).len()/(m0*3e5);
                assert!(final_energy.is_finite()&&final_energy<1e-3,
                    "energy drift {final_energy} at eta={eta}");
                assert!(final_p.is_finite()&&final_p<1e-8,
                    "momentum drift {final_p} at eta={eta}");
                assert!((sph3::total_mass_kg(&central,&particles)-m0).abs()/
                    m0<1e-12);
            }
        }
        let radius_bias=100.0*(profile.radius_max_m/star.radius_m-1.0);
        let pass=three_dyn&&radius_bias.abs()<5.0&&
            outer_drift<0.05&&r99_drift<0.05&&max_cdf<0.05;
        let wall=start.elapsed().as_secs_f64();
        println!("  eta={eta:.3}: R={:.2}km ({radius_bias:+.2}%), R99={:.2}km, CDF={:.4}, g_ref_RMS={gravity_rms:.4}, inner/outer_g={inner:.3}/{outer:.3}, virial={:.4}, force={:.4}",
            profile.radius_max_m/1e3,profile.r99_m/1e3,
            profile.radial_quantile_rms,
            balance.pressure_to_gravity_virial,
            balance.residual_to_gravity_rms);
        println!("    omega_min={:.3}, h_range=[{:.3},{:.3}]Mm, mean_neighbors={neighbors:.1}, hold_dRmax={:.2}%, dR99={:.2}%, maxCDF={max_cdf:.4}, dE={final_energy:.3e}, pass={pass}, wall={wall:.2}s",
            report.minimum_omega,report.minimum_h_m/1e6,
            report.maximum_h_m/1e6,100.0*outer_drift,100.0*r99_drift);
        outcomes.push(Case{
            eta,rmax_post_km:profile.radius_max_m/1e3,
            radius_bias_pct:radius_bias,r99_post_km:profile.r99_m/1e3,
            mass_cdf_post:profile.radial_quantile_rms,
            gravity_reference_rms_post:gravity_rms,
            force_ratio_post:balance.residual_to_gravity_rms,
            virial_post:balance.pressure_to_gravity_virial,
            inner_gravity_reference_ratio:inner,outer_gravity_reference_ratio:outer,
            omega_min:report.minimum_omega,cap_fraction:report.capped_fraction,
            mean_neighbors:neighbors,
            largest_outer_drift_pct:100.0*outer_drift,
            largest_r99_drift_pct:100.0*r99_drift,
            max_cdf_rms:max_cdf,final_energy_drift:final_energy,
            final_momentum_drift:final_p,elapsed_wall_s:wall,
            completed_three_dyn:three_dyn,physical_gate_pass:pass,
        });
    }
    for x in &outcomes {
        let _=writeln!(csv,
            "{n},{:.4},{:.5},{:.5},{:.5},{:.5},{:.6},{:.5},{:.8},{:.8},{:.8},{:.8},{:.6},{:.6},{:.8},{:.6},{:.3},{:.6},{:.6},{:.8},{:.9e},{:.9e},{:.4},{},{}",
            x.eta,relax_steps as f64*dt,hold_steps as f64*dt,
            star.radius_m/1e3,x.rmax_post_km,x.radius_bias_pct,
            x.r99_post_km,x.mass_cdf_post,x.gravity_reference_rms_post,
            x.force_ratio_post,x.virial_post,
            x.inner_gravity_reference_ratio,x.outer_gravity_reference_ratio,
            x.omega_min,x.cap_fraction,x.mean_neighbors,
            x.largest_outer_drift_pct,x.largest_r99_drift_pct,
            x.max_cdf_rms,x.final_energy_drift,x.final_momentum_drift,
            x.elapsed_wall_s,x.completed_three_dyn,x.physical_gate_pass);
    }
    let best=outcomes.iter().min_by(|a,b|{
        // "Best" uses the predeclared scientific score; this is exploratory
        // selection on the same runs, not an independent validation set.
        let score=|x:&&Case|x.radius_bias_pct.abs()+
            x.largest_outer_drift_pct+x.largest_r99_drift_pct+
            100.0*x.max_cdf_rms;
        score(a).total_cmp(&score(b))
    }).unwrap();
    println!("  exploratory score minimum: eta={:.3} (needs independent confirmation); physical pass={}",
        best.eta,best.physical_gate_pass);
    if let Some(path)=option(&args,"--csv") {
        std::fs::write(path,csv).expect("write eta sweep summary");
        println!("  wrote {path}");
    }
    if let Some(path)=option(&args,"--shell-csv") {
        std::fs::write(path,shell_csv).expect("write radial shell comparison");
        println!("  wrote {path}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn eta_sweep_parser_preserves_order_and_bounds() {
        let input=vec![
            "--etas".into(),"0.85,1.20,1.45".into()
        ];
        assert_eq!(eta_list(&input),vec![0.85,1.20,1.45]);
    }
    #[test]
    #[should_panic(expected="--etas must be strictly increasing")]
    fn duplicate_eta_is_rejected() {
        eta_list(&["--etas".into(),"1.2,1.2".into()]);
    }
    #[test]
    fn independent_hydrostatic_gravity_reference_is_finite() {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let p=sph_sampling::seeded_star(
            &star,V3::default(),V3::default(),64,
            StarSampling::ApproximatelyCubicShells
        );
        let mut cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_cubic_gravity();
        cfg.eta=1.2;
        let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
        let report=stellar_balance::radial_balance(
            &p,&mut state,Sph3Config::cold_white_dwarf(2.0).unwrap(),
            cfg,&star,8);
        let value=gravity_reference_rms(&report);
        assert!(value.is_finite()&&value>=0.0&&value<10.0);
    }
}
