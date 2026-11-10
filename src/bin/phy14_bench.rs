//! BG-PHY-14: growing-neighbor SPH scaling and independent shell layouts.
 //!
 //! cargo run --release --bin phy14_bench -- \
 //!  --n 512 --seeds 0,17 --relax-steps 300 --hold-steps 788 \
 //!  --csv bg-phy14-n512.csv
 //!
 //! Fixed eta=1.2 is compared to eta(N)=1.05*(N/256)^(1/6).
 //! This implements the N_neighbors~N^(1/2) compromise in 3D.
 //! A per-particle safeguarded Newton solve audits the implicit rho-h
 //! residual at the SAME final fluid positions. Root integration can be
 //! separately enabled with --root-dynamics for a small-N pilot.
 //!
 //! Physical acceptance: post-damping radius within 5% of independent
 //! hydrostatic R, worst Rmax and R99 drift <=5% over >=3 t_dyn, radial
 //! CDF RMS <=5%. A short CI smoke cannot pass the physical gate.
use blighted_galaxy::physics::{G_SI,SOLAR_MASS_KG};
use blighted_galaxy::sph3::{self,Massive3,Sph3Config,V3};
use blighted_galaxy::sph_adaptive::{
    self,AdaptiveConfig,AdaptiveState,DensityMethod,
};
use blighted_galaxy::sph_sampling::seed_cubic_shells_oriented;
use blighted_galaxy::sph_neighbors::CellNeighbors;
use blighted_galaxy::sph_massprofile;
use blighted_galaxy::sph_hroot::{self,RootLocation};
use blighted_galaxy::stellar_relax::{
    relax_isolated_adaptive,RelaxParameters,
};
use blighted_galaxy::stellar_balance;
use blighted_galaxy::white_dwarf::{DegenerateElectrons,profile_for_mass};
use std::fmt::Write as _;
use std::time::Instant;

fn value<'a>(args:&'a [String],flag:&str)->Option<&'a str>{
    args.windows(2).find(|w|w[0]==flag).map(|w|w[1].as_str())
}
fn opt<T:std::str::FromStr>(args:&[String],flag:&str,default:T)->T {
    value(args,flag).and_then(|x|x.parse().ok()).unwrap_or(default)
}
fn seeds(args:&[String])->Vec<u64>{
    let input=value(args,"--seeds").unwrap_or("0,17");
    let list=input.split(',').map(|x|x.parse::<u64>()
        .expect("invalid integer seed")).collect::<Vec<_>>();
    assert!(!list.is_empty()&&list.len()<=8);
    assert!(list.windows(2).all(|w|w[0]<w[1]));
    list
}
pub fn growth_eta(n:usize,base_n:usize,eta_at_base:f64)->f64 {
    assert!(n>0&&base_n>0&&eta_at_base>0.5&&eta_at_base.is_finite());
    eta_at_base*(n as f64/base_n as f64).powf(1.0/6.0)
}
fn true_neighbors(
    p:&[sph3::Fluid3],state:&AdaptiveState,
)->f64 {
    let hmax=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let candidates=CellNeighbors::new(p,hmax).all_candidates();
    candidates.iter().enumerate().map(|(i,list)|
        list.iter().filter(|&&j|(p[i].pos_m-p[j].pos_m).len_sq()<
            (2.0*state.smoothing_m[i]).powi(2)).count()
    ).sum::<usize>() as f64/p.len() as f64
}
fn hydro_gravity_error(b:&stellar_balance::StellarBalance)->f64{
    let (num,den)=b.shells.iter().fold((0.0,0.0),|(a,d),shell|{
        let weight=shell.mass_fraction_hi-shell.mass_fraction_lo;
        (a+weight*(shell.gravity_in_m_s2-
            shell.hydrostatic_gravity_in_m_s2).powi(2),
         d+weight*shell.hydrostatic_gravity_in_m_s2.powi(2))
    });
    (num/den).sqrt()
}
#[derive(Debug,Clone,Copy)]
struct ClosureAudit{
    max_free_before:f64,
    rms_free_before:f64,
    max_free_after:f64,
    cap_fraction:f64,
    max_delta_h_relative:f64,
    mean_iterations:f64,
    /// Relative RMS difference in SPH accelerations when using exact
    /// safeguarded root h values at the same frozen particle positions.
    force_discrepancy:f64,
}
fn closure_audit(
    p:&[sph3::Fluid3],state:&AdaptiveState,
    cfg:AdaptiveConfig,physics:Sph3Config,
)->ClosureAudit {
    let full=CellNeighbors::new(p,cfg.max_h_m).all_candidates();
    let mut initial_sum=0.0;
    let mut max_free=0.0_f64;
    let mut n_free=0_usize;
    let mut after_free=0.0_f64;
    let mut caps=0_usize;
    let mut max_delta=0.0_f64;
    let mut iterations=0_usize;
    let mut corrected=state.clone();
    for i in 0..p.len(){
        let h0=state.smoothing_m[i];
        let root=sph_hroot::solve_particle(
            i,p,&full[i],cfg.eta,cfg.min_h_m,cfg.max_h_m,h0,
            1e-9,64);
        let before=sph_hroot::closure_residual(
            i,p,&full[i],cfg.eta,h0);
        if root.location==RootLocation::Interior {
            n_free+=1;
            initial_sum+=before*before;
            max_free=max_free.max(before);
            after_free=after_free.max(root.residual);
        } else {caps+=1;}
        max_delta=max_delta.max((root.h_m-h0).abs()/h0);
        iterations+=root.iterations;
        corrected.smoothing_m[i]=root.h_m;
        corrected.densities_kg_m3[i]=root.density_kg_m3;
        corrected.omega[i]=root.omega;
    }
    let center=Massive3{
        pos_m:V3::default(),vel_m_s:V3::default(),
        mass_kg:0.0,radius_m:0.0,
    };
    // Evaluate the total SPH acceleration with both closures at the
    // same coordinates; any change is solely a closure numerical effect.
    let mut old_state=state.clone();
    let old_force=sph_adaptive::forces(
        &center,p,&mut old_state,physics,cfg);
    let mut root_state=corrected;
    let root_force=sph_adaptive::forces(
        &center,p,&mut root_state,physics,cfg.with_newton_density());
    let diff=old_force.accelerations_m_s2.iter()
        .zip(&root_force.accelerations_m_s2)
        .zip(p).map(|((a,b),x)|x.mass_kg*(*a-*b).len_sq())
        .sum::<f64>();
    let scale=root_force.accelerations_m_s2.iter().zip(p)
        .map(|(a,x)|x.mass_kg*a.len_sq()).sum::<f64>();
    ClosureAudit{
        max_free_before:max_free,
        rms_free_before:(initial_sum/n_free.max(1) as f64).sqrt(),
        max_free_after:after_free,
        cap_fraction:caps as f64/p.len() as f64,
        max_delta_h_relative:max_delta,
        mean_iterations:iterations as f64/p.len() as f64,
        force_discrepancy:(diff/scale.max(1e-30)).sqrt(),
    }
}
fn main(){
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    let n:usize=opt(&args,"--n",128);
    let relax_steps:usize=opt(&args,"--relax-steps",150);
    let hold_steps:usize=opt(&args,"--hold-steps",262);
    let dt:f64=opt(&args,"--dt",0.01);
    let fixed_eta:f64=opt(&args,"--fixed-eta",1.2);
    let base_eta:f64=opt(&args,"--growth-base-eta",1.05);
    let base_n:usize=opt(&args,"--growth-base-n",256);
    let root_dynamics=args.iter().any(|x|x=="--root-dynamics");
    assert!(n>=32&&n<=2048&&n%2==0&&base_n>0);
    assert!(relax_steps>0&&relax_steps<=30000);
    assert!(hold_steps>0&&hold_steps<=30000);
    assert!(dt.is_finite()&&dt>0.0&&dt<=0.1);
    let etas=[
        ("fixed",fixed_eta),
        ("growth",growth_eta(n,base_n,base_eta))
    ];
    assert!(etas.iter().all(|x|x.1>0.5&&x.1<3.0));
    let eos=DegenerateElectrons::new(2.0).unwrap();
    let star=profile_for_mass(
        0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
    let physics=Sph3Config::cold_white_dwarf(2.0).unwrap();
    let t_dyn=(star.radius_m.powi(3)/(G_SI*star.mass_kg)).sqrt();
    let valid_hold=hold_steps as f64*dt>=3.0*t_dyn;
    println!("BG-PHY-14: growing-neighbor scaling and implicit h/rho audit");
    println!("N={n}, Rref={:.3}km, t_dyn={t_dyn:.5}s, relax={:.3}s, hold={:.3}s, root_dynamics={root_dynamics}",
        star.radius_m/1e3,relax_steps as f64*dt,
        hold_steps as f64*dt);
    if !valid_hold {println!("SHORT HOLD: physical acceptance disabled");}
    let mut csv=String::from(
        "N,seed,law,eta,neighbor_count,relax_seconds,hold_seconds,post_Rmax_km,post_R99_km,post_radius_bias_pct,post_cdf_rms,post_gravity_ref_rms,post_net_gravity_ratio,max_Rmax_drift_pct,max_R99_drift_pct,max_CDF_rms,final_energy_drift,final_momentum_drift,max_free_h_residual_before,rms_free_h_residual_before,max_free_h_residual_after,h_capped_fraction,max_delta_h_fraction,mean_root_iterations,rms_root_force_delta,wall_seconds,full_hold,physical_gate_pass\n"
    );
    for seed in seeds(&args){
        let (base,_)=seed_cubic_shells_oriented(&star,
            V3::new(1e8,0.0,0.0),V3::default(),n,seed);
        for (law,eta) in etas {
            let mut p=base.clone();
            let mut cfg=AdaptiveConfig::from_star_radius(
                star.radius_m).with_cubic_gravity();
            cfg.eta=eta;
            if root_dynamics {cfg=cfg.with_newton_density();}
            let mut state=AdaptiveState::new(&p,star.radius_m,cfg);
            let started=Instant::now();
            let relaxed=relax_isolated_adaptive(
                &mut p,physics,&star,cfg,&mut state,
                RelaxParameters{
                    iterations:relax_steps,dt_s:dt,damping_time_s:0.7,
                }).expect("finite relaxation");
            let mass_profile=sph_massprofile::profile_diagnostics(&p,&star);
            let balance=stellar_balance::radial_balance(
                &p,&mut state,physics,cfg,&star,8);
            let g_error=hydro_gravity_error(&balance);
            let neighbors=true_neighbors(&p,&state);
            let closure=closure_audit(&p,&state,cfg,physics);
            let mut central=Massive3{
                pos_m:V3::default(),vel_m_s:V3::default(),
                mass_kg:0.0,radius_m:0.0,
            };
            let mass0=sph3::total_mass_kg(&central,&p);
            let momentum0=sph3::total_momentum(&central,&p);
            let energy0=sph_adaptive::adaptive_energy_j(
                &central,&p,&mut state,physics,cfg);
            let mut worst_rmax=0.0_f64;
            let mut worst_r99=0.0_f64;
            let mut worst_cdf=mass_profile.radial_quantile_rms;
            let mut final_de=0.0;
            let mut final_dp=0.0;
            let stride=(hold_steps/30).max(1);
            for step in 1..=hold_steps {
                sph_adaptive::advance(
                    &mut central,&mut p,&mut state,physics,cfg,dt);
                if step%stride==0 || step==hold_steps {
                    let q=sph_massprofile::profile_diagnostics(&p,&star);
                    worst_rmax=worst_rmax.max(
                        (q.radius_max_m/mass_profile.radius_max_m-1.0).abs());
                    worst_r99=worst_r99.max(
                        (q.r99_m/mass_profile.r99_m-1.0).abs());
                    worst_cdf=worst_cdf.max(q.radial_quantile_rms);
                    let e=sph_adaptive::adaptive_energy_j(
                        &central,&p,&mut state,physics,cfg);
                    final_de=(e-energy0).abs()/energy0.abs().max(1.0);
                    let pp=sph3::total_momentum(&central,&p);
                    final_dp=(pp-momentum0).len()/(mass0*3e5);
                    assert!(final_de.is_finite()&&final_de<1e-3);
                    assert!(final_dp.is_finite()&&final_dp<1e-8);
                    assert!((sph3::total_mass_kg(&central,&p)-mass0).abs()/
                        mass0<1e-12);
                }
            }
            let radius_bias=100.0*(mass_profile.radius_max_m/
                star.radius_m-1.0);
            let pass=valid_hold&&radius_bias.abs()<5.0&&
                worst_rmax<0.05&&worst_r99<0.05&&worst_cdf<0.05;
            let elapsed=started.elapsed().as_secs_f64();
            println!("seed={seed}, {law:6}, eta={eta:.4}, neighbors={neighbors:.1}, R={:.2}km ({radius_bias:+.2}%), g_ref_rms={g_error:.4}, net/grav={:.4}, hold Rmax/R99={:.2}%/{:.2}%, CDFmax={worst_cdf:.4}, gate={pass}, wall={elapsed:.1}s",
                mass_profile.radius_max_m/1e3,
                balance.residual_to_gravity_rms,
                100.0*worst_rmax,100.0*worst_r99);
            println!("  h closure max free residual {:.3e} -> {:.3e}, caps {:.1}%, max delta h {:.3e}, force discrepancy {:.3e}, root iterations {:.2}",
                closure.max_free_before,closure.max_free_after,
                100.0*closure.cap_fraction,
                closure.max_delta_h_relative,
                closure.force_discrepancy,closure.mean_iterations);
            let _=writeln!(csv,
                "{n},{seed},{law},{eta:.6},{neighbors:.4},{:.5},{:.5},{:.5},{:.5},{radius_bias:.6},{:.8},{g_error:.8},{:.8},{:.6},{:.6},{worst_cdf:.8},{final_de:.9e},{final_dp:.9e},{:.9e},{:.9e},{:.9e},{:.7},{:.9e},{:.5},{:.9e},{elapsed:.4},{valid_hold},{pass}",
                relax_steps as f64*dt,hold_steps as f64*dt,
                mass_profile.radius_max_m/1e3,mass_profile.r99_m/1e3,
                mass_profile.radial_quantile_rms,
                balance.residual_to_gravity_rms,
                100.0*worst_rmax,100.0*worst_r99,
                closure.max_free_before,closure.rms_free_before,
                closure.max_free_after,closure.cap_fraction,
                closure.max_delta_h_relative,closure.mean_iterations,
                closure.force_discrepancy);
            assert!(relaxed.report.mass_relative_error<1e-12);
        }
    }
    if let Some(path)=value(&args,"--csv"){
        std::fs::write(path,csv).expect("write PHY14 CSV");
        println!("wrote {path}");
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn neighbor_count_scaling_law_is_consistent() {
        let e0=growth_eta(256,256,1.05);
        let e1=growth_eta(512,256,1.05);
        let e2=growth_eta(1024,256,1.05);
        assert!((e0-1.05).abs()<1e-12);
        assert!((e1/e0-2.0_f64.powf(1.0/6.0)).abs()<1e-12);
        assert!((e2/e0-4.0_f64.powf(1.0/6.0)).abs()<1e-12);
        assert!(e0<e1&&e1<e2);
        // The expected homogeneous neighbor count is proportional to eta³.
        assert!(((e2/e0).powi(3)-2.0).abs()<1e-12);
    }
    #[test]
    #[should_panic]
    fn seed_order_must_be_unique() {
        let a=vec!["--seeds".to_string(),"17,17".to_string()];
        seeds(&a);
    }
}
