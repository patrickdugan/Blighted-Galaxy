//! BG-PHY-17: calibrate eta to REALIZED kernel-support neighbours.
//!
//! The same physical particles and kernel solve the implicit h/rho equation
//! at each trial eta. We use the safeguarded Newton root during calibration,
//! so a transient 24-pass closure residual cannot move a target neighbour
//! count. Strict support is r < 2*h_i and includes the self particle.
//! Targets are approximate because counts are discrete and particle h caps
//! can bind. A returned unreachable target is explicitly marked.
use crate::sph3::Fluid3;
use crate::sph_adaptive::{
    self,AdaptiveConfig,AdaptiveState,
};
use crate::sph_neighbors::CellNeighbors;

#[derive(Clone,Copy,Debug)]
pub struct NeighborStats {
    pub mean:f64,
    pub min:usize,
    pub max:usize,
    pub inner_mean:f64,
    pub outer_mean:f64,
    pub lower_h_caps:f64,
    pub upper_h_caps:f64,
    pub min_omega:f64,
    pub mean_h_m:f64,
}
#[derive(Clone,Copy,Debug)]
pub struct CalibratedEta {
    pub eta:f64,
    pub neighbors:NeighborStats,
    pub target:f64,
    pub abs_error:f64,
    pub bracketed:bool,
    pub evaluations:usize,
}
/// Counts support neighbours with the smoothing radii in the supplied state.
/// Shell stratification uses sorted particle radius about the mass centroid.
pub fn measure(p:&[Fluid3],state:&AdaptiveState,cfg:AdaptiveConfig)->NeighborStats {
    assert!(!p.is_empty()&&state.valid_for(p.len())&&cfg.valid());
    let hmax=state.smoothing_m.iter().copied().fold(0.0_f64,f64::max);
    let lists=CellNeighbors::new(p,hmax).all_candidates();
    let center=p.iter().fold(crate::sph3::V3::default(),|sum,x|
        sum+x.pos_m*x.mass_kg)/p.iter().map(|x|x.mass_kg).sum::<f64>();
    let mut rank=(0..p.len()).collect::<Vec<_>>();
    rank.sort_by(|&i,&j| {
        (p[i].pos_m-center).len_sq().total_cmp(&(p[j].pos_m-center).len_sq())
            .then(i.cmp(&j))
    });
    let mut count=vec![0_usize;p.len()];
    for (i,js) in lists.iter().enumerate() {
        let h=state.smoothing_m[i];
        count[i]=js.iter().filter(|&&j|
            (p[i].pos_m-p[j].pos_m).len_sq()<(2.0*h).powi(2)).count();
    }
    let n=p.len();
    let quarter=(n/4).max(1);
    let inner=rank.iter().take(quarter).map(|&i|count[i]).sum::<usize>();
    let outer=rank.iter().rev().take(quarter).map(|&i|count[i]).sum::<usize>();
    let lower=state.smoothing_m.iter().filter(|&&h|
        (h-cfg.min_h_m).abs()<=cfg.min_h_m*1e-8).count();
    let upper=state.smoothing_m.iter().filter(|&&h|
        (h-cfg.max_h_m).abs()<=cfg.max_h_m*1e-8).count();
    NeighborStats {
        mean:count.iter().sum::<usize>() as f64/n as f64,
        min:*count.iter().min().unwrap(),
        max:*count.iter().max().unwrap(),
        inner_mean:inner as f64/quarter as f64,
        outer_mean:outer as f64/quarter as f64,
        lower_h_caps:lower as f64/n as f64,
        upper_h_caps:upper as f64/n as f64,
        min_omega:state.omega.iter().copied().fold(f64::INFINITY,f64::min),
        mean_h_m:state.smoothing_m.iter().sum::<f64>()/n as f64,
    }
}
pub fn evaluate_eta(
    p:&[Fluid3],radius_m:f64,base_cfg:AdaptiveConfig,eta:f64
)->NeighborStats {
    let mut cfg=base_cfg.with_newton_density();
    cfg.eta=eta;
    assert!(cfg.valid());
    let mut state=AdaptiveState::new(p,radius_m,cfg);
    sph_adaptive::update_adaptive_densities(p,&mut state,cfg);
    measure(p,&state,cfg)
}
/// Bisection with a bounded number of eta evaluations. This also returns the
/// closest measured eta in case the target lies outside the attainable range.
/// Count is discrete and the mean can have small non-monotone steps.
pub fn calibrate(
    p:&[Fluid3],radius_m:f64,cfg:AdaptiveConfig,
    target:f64,eta_min:f64,eta_max:f64,rounds:usize,
)->CalibratedEta {
    assert!(target>=2.0&&target<=p.len() as f64);
    assert!(eta_min>0.5&&eta_max<3.0&&eta_max>eta_min);
    assert!(rounds>=1&&rounds<=20);
    let low_stats=evaluate_eta(p,radius_m,cfg,eta_min);
    let high_stats=evaluate_eta(p,radius_m,cfg,eta_max);
    let bracketed=low_stats.mean<=target&&high_stats.mean>=target;
    let mut best=CalibratedEta {
        eta:eta_min,neighbors:low_stats,target,
        abs_error:(low_stats.mean-target).abs(),
        bracketed,evaluations:2,
    };
    let high_error=(high_stats.mean-target).abs();
    if high_error<best.abs_error {
        best.eta=eta_max;
        best.neighbors=high_stats;
        best.abs_error=high_error;
    }
    if !bracketed {return best;}
    let mut lo=eta_min;
    let mut hi=eta_max;
    for _ in 0..rounds {
        let eta=0.5*(lo+hi);
        let stats=evaluate_eta(p,radius_m,cfg,eta);
        let error=(stats.mean-target).abs();
        best.evaluations+=1;
        if error<best.abs_error {
            best.eta=eta;
            best.abs_error=error;
            best.neighbors=stats;
        }
        if stats.mean<target {lo=eta;} else {hi=eta;}
    }
    best
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph_sampling::seed_cubic_shells_oriented;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    use crate::sph3::V3;
    #[test]
    fn calibrates_real_support_counts_with_exact_roots() {
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(
            &star,V3::new(1e8,0.0,0.0),V3::default(),96,17).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4();
        let a=evaluate_eta(&p,star.radius_m,cfg,1.15);
        let b=evaluate_eta(&p,star.radius_m,cfg,1.65);
        assert!(a.mean<b.mean,"support must grow in reference sample");
        let target=0.5*(a.mean+b.mean);
        let c=calibrate(&p,star.radius_m,cfg,target,1.15,1.65,10);
        assert!(c.bracketed);
        assert!(c.neighbors.mean>=a.mean&&c.neighbors.mean<=b.mean);
        assert!(c.abs_error<(b.mean-a.mean)*0.12+1.0);
        assert!(c.neighbors.inner_mean.is_finite());
        assert!(c.neighbors.outer_mean.is_finite());
        assert!(c.neighbors.min>=1);
        assert!(c.neighbors.max<=p.len());
    }
    #[test]
    fn unreachable_target_reports_bound_without_claiming_success() {
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let p=seed_cubic_shells_oriented(
            &star,V3::new(1e8,0.0,0.0),V3::default(),64,0).0;
        let cfg=AdaptiveConfig::from_star_radius(star.radius_m)
            .with_wendland_c4();
        let c=calibrate(&p,star.radius_m,cfg,64.0,0.8,0.9,6);
        assert!(!c.bracketed);
        assert!(c.eta==0.8||c.eta==0.9);
    }
}
