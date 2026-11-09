//! Shell-stratified equal-mass sampling of a cold 1-D hydrostatic star.
 //!
 //! Motivation: pairwise-antipodal Fibonacci directions correlated radius
 //! and angle in the BG-PHY-08 seed. Here shells are binned in physical
 //! radius, and each shell gets independently rotated equal-area Fibonacci
 //! directions. Each antipodal pair occupies the same radial mass quantile.
 //!
 //! This is a deterministic low-discrepancy initializer rather than a
 //! relaxed HEALPix/glass or an exact discrete equilibrium configuration.
use crate::sph3::{Fluid3,V3};
use crate::white_dwarf::WhiteDwarfProfile;
use std::f64::consts::PI;

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum StarSampling {
    LegacyAntipodalFibonacci,
    ShellBalancedAntipodal,
}
#[derive(Debug,Clone)]
pub struct ShellSeedReport {
    pub shell_count:usize,
    pub pair_counts:Vec<usize>,
    pub radius_min_m:f64,
    pub radius_max_m:f64,
    pub mass_relative_error:f64,
    pub center_error_m:f64,
}
fn shell_count_for_particles(count:usize)->usize {
    let pairs=count/2;
    (((pairs as f64).cbrt()*1.2).round() as usize)
        .clamp(2,pairs.max(2))
}
fn radical_inverse(mut n:usize,base:usize)->f64 {
    let mut f=1.0/base as f64;
    let mut result=0.0;
    while n>0 {
        result+=(n%base) as f64*f;
        n/=base;
        f/=base as f64;
    }
    result
}
pub fn seed_shell_balanced(
    profile:&WhiteDwarfProfile,
    center_m:V3,bulk_velocity_m_s:V3,count:usize,
)->(Vec<Fluid3>,ShellSeedReport) {
    assert!(count>=8 && count%2==0);
    let n_pairs=count/2;
    let shells=shell_count_for_particles(count);
    let mut radial_bins=vec![Vec::<(usize,f64)>::new();shells];
    for i in 0..n_pairs {
        let mass_fraction=(i as f64+0.5)/n_pairs as f64;
        let r=profile.radius_at_mass_fraction(mass_fraction);
        let idx=((r/profile.radius_m)*shells as f64) as usize;
        radial_bins[idx.min(shells-1)].push((i,r));
    }
    let each_mass=profile.mass_kg/count as f64;
    let golden=PI*(3.0-5.0_f64.sqrt());
    let mut pairs=Vec::<(usize,Fluid3,Fluid3)>::with_capacity(n_pairs);
    let mut pair_counts=Vec::with_capacity(shells);
    for (shell,group) in radial_bins.into_iter().enumerate() {
        let size=group.len();
        pair_counts.push(size);
        if size==0 {continue;}
        // Independent phases decorrelate the spherical placement from
        // the monotonically increasing radial mass CDF.
        let phase=2.0*PI*radical_inverse(shell+1,3);
        for (j,(i,r)) in group.into_iter().enumerate() {
            let z=1.0-2.0*(j as f64+0.5)/size as f64;
            let phi=golden*j as f64+phase;
            let xy=(1.0-z*z).sqrt();
            let direction=V3::new(xy*phi.cos(),xy*phi.sin(),z);
            let delta=direction*r;
            let plus=Fluid3{pos_m:center_m+delta,
                vel_m_s:bulk_velocity_m_s,mass_kg:each_mass};
            let minus=Fluid3{pos_m:center_m-delta,
                vel_m_s:bulk_velocity_m_s,mass_kg:each_mass};
            pairs.push((i,plus,minus));
        }
    }
    // Maintain radial quantile order for replay/benchmark comparisons.
    pairs.sort_by_key(|(i,_,_)|*i);
    let mut particles=Vec::with_capacity(count);
    for (_,a,b) in pairs {
        particles.push(a);particles.push(b);
    }
    let total=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    let weighted_pos=particles.iter().fold(V3::default(),|s,p|
        s+p.pos_m*p.mass_kg)/total;
    let rmax=particles.iter().map(|p|(p.pos_m-center_m).len())
        .fold(0.0_f64,f64::max);
    let rmin=particles.iter().map(|p|(p.pos_m-center_m).len())
        .fold(f64::INFINITY,f64::min);
    let report=ShellSeedReport {
        shell_count:shells,pair_counts,
        radius_min_m:rmin,radius_max_m:rmax,
        mass_relative_error:(total-profile.mass_kg).abs()/profile.mass_kg,
        center_error_m:(weighted_pos-center_m).len(),
    };
    (particles,report)
}
pub fn seeded_star(
    profile:&WhiteDwarfProfile,center:V3,velocity:V3,
    count:usize,method:StarSampling,
)->Vec<Fluid3> {
    match method {
        StarSampling::LegacyAntipodalFibonacci=>
            crate::sph3::seed_from_hydrostatic(profile,center,velocity,count),
        StarSampling::ShellBalancedAntipodal=>
            seed_shell_balanced(profile,center,velocity,count).0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn setup()->WhiteDwarfProfile {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap()
    }
    #[test]
    fn shell_seed_preserves_mass_center_and_bulk_momentum() {
        let profile=setup();
        let center=V3::new(1e8,4e7,-2e7);
        let bulk=V3::new(50_000.0,-10_000.0,200.0);
        for n in [64,128,512] {
            let (p,report)=seed_shell_balanced(&profile,center,bulk,n);
            assert_eq!(p.len(),n);
            assert_eq!(report.pair_counts.iter().sum::<usize>(),n/2);
            assert!(report.mass_relative_error<1e-13);
            assert!(report.center_error_m<1e-5);
            assert!(report.radius_min_m>0.0);
            assert!(report.radius_max_m<profile.radius_m);
            let mass=p.iter().map(|p|p.mass_kg).sum::<f64>();
            let momentum=p.iter().fold(V3::default(),
                |sum,x|sum+x.vel_m_s*x.mass_kg)/mass;
            assert!((momentum-bulk).len()<1e-6);
        }
    }
    #[test]
    fn shell_seed_keeps_interior_mass_quantiles() {
        let profile=setup();
        let n=128;
        let (p,_)=seed_shell_balanced(&profile,
            V3::default(),V3::default(),n);
        for i in 0..n/2 {
            let r=p[2*i].pos_m.len();
            let expected=profile.radius_at_mass_fraction(
                (i as f64+0.5)/(n/2) as f64);
            assert!((r-expected).abs()<1e-7);
            assert!((p[2*i].pos_m+p[2*i+1].pos_m).len()<1e-7);
        }
    }
    #[test]
    fn initializer_is_reproducible_at_all_reference_sizes() {
        let profile=setup();
        for n in [64,256,512] {
            let a=seeded_star(&profile,V3::default(),V3::default(),
                n,StarSampling::ShellBalancedAntipodal);
            let b=seeded_star(&profile,V3::default(),V3::default(),
                n,StarSampling::ShellBalancedAntipodal);
            assert_eq!(a[0].pos_m,b[0].pos_m);
            assert_eq!(a[n-1].pos_m,b[n-1].pos_m);
        }
    }
}
