//! Resolution-robust radial diagnostics for a discrete white-dwarf SPH star.
 //!
 //! The previous "radius" used the single outermost particle. This module
 //! retains that value and adds weighted radial mass quantiles, a mass-CDF
 //! mismatch and a robust r_99. The hydrostatic reference is the independent
 //! 1D cold-electron stellar structure solver.
use crate::sph3::{Fluid3,V3};
use crate::white_dwarf::WhiteDwarfProfile;
#[derive(Clone,Copy,Debug)]
pub struct MassProfileDiagnostic {
    pub radius_max_m:f64,
    pub r50_m:f64,
    pub r90_m:f64,
    pub r95_m:f64,
    pub r99_m:f64,
    /// RMS radius-quantile deviation normalized by reference R, evaluated
    /// at 19 mass-CDF fractions from 0.05 to 0.95.
    pub radial_quantile_rms:f64,
    pub center_m:V3,
    pub total_mass_kg:f64,
}
fn weighted_quantile(sorted:&[(f64,f64)],mass:f64,f:f64)->f64 {
    if f<=0.0 {return 0.0;}
    if f>=1.0 {return sorted.last().map(|x|x.0).unwrap_or(0.0);}
    let mut prev_mass_fraction=0.0;
    let mut prev_r=0.0;
    let mut cumulative=0.0;
    for (r,m) in sorted {
        let at=(cumulative+0.5*m)/mass;
        if f<=at {
            let w=if at>prev_mass_fraction {
                ((f-prev_mass_fraction)/(at-prev_mass_fraction)).clamp(0.0,1.0)
            } else {0.0};
            return prev_r+w*(*r-prev_r);
        }
        prev_mass_fraction=at;
        prev_r=*r;
        cumulative+=m;
    }
    sorted.last().unwrap().0
}
pub fn profile_diagnostics(
    particles:&[Fluid3],reference:&WhiteDwarfProfile,
)->MassProfileDiagnostic {
    assert!(!particles.is_empty());
    let mass=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    assert!(mass>0.0&&mass.is_finite());
    let center=particles.iter().fold(V3::default(),|sum,p|
        sum+p.pos_m*p.mass_kg)/mass;
    let mut radii=particles.iter().map(|p|
        ((p.pos_m-center).len(),p.mass_kg)).collect::<Vec<_>>();
    radii.sort_by(|a,b|a.0.total_cmp(&b.0));
    let q=|f|weighted_quantile(&radii,mass,f);
    let mut squared=0.0;
    for i in 1..=19 {
        let fraction=i as f64/20.0;
        let observed=q(fraction);
        let expected=reference.radius_at_mass_fraction(fraction);
        squared+=((observed-expected)/reference.radius_m).powi(2);
    }
    MassProfileDiagnostic {
        radius_max_m:radii.last().unwrap().0,
        r50_m:q(0.50),r90_m:q(0.90),
        r95_m:q(0.95),r99_m:q(0.99),
        radial_quantile_rms:(squared/19.0).sqrt(),
        center_m:center,total_mass_kg:mass,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sph_sampling::{self,StarSampling};
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn setup(n:usize)->(WhiteDwarfProfile,Vec<Fluid3>) {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        let p=sph_sampling::seeded_star(&star,V3::new(1e8,-5e7,0.0),
            V3::default(),n,StarSampling::ApproximatelyCubicShells);
        (star,p)
    }
    #[test]
    fn near_cubic_cdf_quantiles_match_initial_hydrostatic_mass_profile() {
        for n in [64,128,512,1024] {
            let (profile,p)=setup(n);
            let observed=profile_diagnostics(&p,&profile);
            assert!((observed.total_mass_kg-profile.mass_kg).abs()/
                profile.mass_kg<1e-12);
            assert!((observed.center_m-V3::new(1e8,-5e7,0.0)).len()<1e-5);
            assert!(observed.r50_m<observed.r90_m &&
                observed.r90_m<=observed.r95_m &&
                observed.r95_m<=observed.r99_m);
            assert!(observed.radial_quantile_rms<0.07);
        }
    }
    #[test]
    fn isotropic_uniform_ball_quantiles_monotonic() {
        let (profile,mut p)=setup(128);
        let before=profile_diagnostics(&p,&profile);
        for x in &mut p {
            x.pos_m=before.center_m+
                (x.pos_m-before.center_m)*1.1;
        }
        let after=profile_diagnostics(&p,&profile);
        assert!(after.r50_m>before.r50_m);
        assert!(after.r99_m>before.r99_m);
        assert!((after.radius_max_m/before.radius_max_m-1.1).abs()<1e-10);
        assert!(after.radial_quantile_rms>before.radial_quantile_rms);
    }
}
