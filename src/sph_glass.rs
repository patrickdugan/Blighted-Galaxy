//! BG-PHY-16: deterministic antipodal angular-glass relaxation.
//!
//! Equal-mass hydrostatic radii and center/bulk momentum remain unchanged.
//! Only directions of paired antipodal particles move on their spheres.
//! A smooth short-range Gaussian repulsion between physical ± pairs
//! lowers angular correlation at fixed radial mass CDF. This is a
//! constrained geometric seed optimizer, not hydrodynamic damping.
use crate::sph3::{Fluid3,V3};
use crate::white_dwarf::WhiteDwarfProfile;

#[derive(Clone,Copy,Debug)]
pub struct AngularGlassReport {
    pub starting_objective:f64,
    pub ending_objective:f64,
    pub accepted_sweeps:usize,
    pub maximum_rotation_degrees:f64,
}
fn pair_energy_and_gradient(
    directions:&[V3],radii:&[f64],hs:&[f64],
)->(f64,Vec<V3>) {
    let count=directions.len();
    assert_eq!(radii.len(),count);
    assert_eq!(hs.len(),count);
    let mut grads=vec![V3::default();count];
    let mut objective=0.0;
    for i in 0..count {
        for j in (i+1)..count {
            let ri=radii[i];
            let rj=radii[j];
            let dot=directions[i].dot(directions[j]).clamp(-1.0,1.0);
            let sigma=0.5*(hs[i]+hs[j]);
            let sigma_sq=sigma*sigma;
            let d0=ri*ri+rj*rj-2.0*ri*rj*dot;
            let d1=ri*ri+rj*rj+2.0*ri*rj*dot;
            let same=(-d0/(2.0*sigma_sq)).exp();
            let flipped=(-d1/(2.0*sigma_sq)).exp();
            objective+=same+flipped;
            let g=(ri*rj/sigma_sq)*(same-flipped);
            grads[i]+=directions[j]*g;
            grads[j]+=directions[i]*g;
        }
    }
    (objective,grads)
}
/// Minimize short-range two-sided pair correlations using projected-gradient
/// descent and deterministic global backtracking. Only the orientation of
/// each pair changes, so mass CDF and antipodal symmetry are exact.
pub fn relax_antipodal_angles(
    source:&[Fluid3],profile:&WhiteDwarfProfile,
    center:V3,eta:f64,requested_sweeps:usize,
)->(Vec<Fluid3>,AngularGlassReport) {
    assert!(source.len()>=8&&source.len()%2==0);
    assert!(eta.is_finite()&&eta>0.5&&eta<3.0);
    assert!(requested_sweeps<=200);
    let count=source.len()/2;
    let m=source[0].mass_kg;
    let mut initial=Vec::with_capacity(count);
    let mut radii=Vec::with_capacity(count);
    let mut hs=Vec::with_capacity(count);
    for pair in source.chunks_exact(2) {
        let v=pair[0].pos_m-center;
        let r=v.len();
        assert!(r>0.0&&r<profile.radius_m);
        assert!((pair[1].pos_m-center+v).len()<1e-7*profile.radius_m);
        assert!((pair[0].mass_kg-m).abs()/m<1e-12);
        assert!((pair[1].mass_kg-m).abs()/m<1e-12);
        initial.push(v*(1.0/r));
        radii.push(r);
        let rho=profile.density_at(r);
        assert!(rho>0.0);
        hs.push((eta*(m/rho).cbrt())
            .clamp(0.08*profile.radius_m,0.90*profile.radius_m));
    }
    let mut directions=initial.clone();
    let (mut objective,_)=pair_energy_and_gradient(&directions,&radii,&hs);
    let starting_objective=objective;
    let mut accepted=0;
    for _ in 0..requested_sweeps {
        let (_,gradients)=pair_energy_and_gradient(&directions,&radii,&hs);
        let tangents=directions.iter().zip(&gradients)
            .map(|(u,g)|*g-*u*g.dot(*u)).collect::<Vec<_>>();
        let mut step=0.24;
        let mut next=Vec::with_capacity(count);
        let mut accepted_candidate=false;
        for _ in 0..12 {
            next.clear();
            for (u,g) in directions.iter().zip(&tangents) {
                let scaled=step/(1.0+g.len());
                next.push((*u-*g*scaled).normalized());
            }
            let (trial,_)=pair_energy_and_gradient(&next,&radii,&hs);
            if trial+1e-13*objective.abs().max(1.0)<objective {
                objective=trial;
                accepted_candidate=true;
                break;
            }
            step*=0.5;
        }
        if !accepted_candidate {break;}
        std::mem::swap(&mut next,&mut directions);
        accepted+=1;
    }
    let mut maximum_rotation=0.0_f64;
    let mut p=source.to_vec();
    for i in 0..count {
        let angle=initial[i].dot(directions[i]).clamp(-1.0,1.0).acos();
        maximum_rotation=maximum_rotation.max(angle.to_degrees());
        let offset=directions[i]*radii[i];
        p[2*i].pos_m=center+offset;
        p[2*i+1].pos_m=center-offset;
    }
    (p,AngularGlassReport {
        starting_objective,ending_objective:objective,
        accepted_sweeps:accepted,maximum_rotation_degrees:maximum_rotation,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    use crate::sph_sampling::seed_cubic_shells_oriented;
    fn fixture(seed:u64,n:usize)->(WhiteDwarfProfile,Vec<Fluid3>,V3) {
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,
            DegenerateElectrons::new(2.0).unwrap(),10_000.0).unwrap();
        let center=V3::new(1e8,2e7,-7e7);
        let p=seed_cubic_shells_oriented(&star,center,
            V3::new(35.0,-12.0,11.0),n,seed).0;
        (star,p,center)
    }
    #[test]
    fn glass_preserves_every_mass_quantile_and_antipodal_invariant() {
        let (star,p,center)=fixture(17,128);
        let (after,report)=relax_antipodal_angles(&p,&star,center,1.2,32);
        assert!(report.ending_objective<=report.starting_objective);
        assert!(report.accepted_sweeps>0);
        for (before,now) in p.iter().zip(&after) {
            assert!(((before.pos_m-center).len()-
                (now.pos_m-center).len()).abs()/star.radius_m<1e-12);
            assert_eq!(before.mass_kg,now.mass_kg);
            assert_eq!(before.vel_m_s,now.vel_m_s);
        }
        for pair in after.chunks_exact(2) {
            assert!((pair[0].pos_m+pair[1].pos_m-center*2.0).len()<1e-7);
        }
    }
    #[test]
    fn glass_is_deterministic_and_zero_sweep_is_identity() {
        let (star,p,center)=fixture(0,64);
        let (a,ar)=relax_antipodal_angles(&p,&star,center,1.2,16);
        let (b,br)=relax_antipodal_angles(&p,&star,center,1.2,16);
        let (none,report)=relax_antipodal_angles(&p,&star,center,1.2,0);
        for i in 0..p.len() {
            assert_eq!(a[i].pos_m,b[i].pos_m);
            assert!((none[i].pos_m-p[i].pos_m).len()/star.radius_m<1e-12);
        }
        assert_eq!(ar.ending_objective,br.ending_objective);
        assert_eq!(report.accepted_sweeps,0);
    }
}
