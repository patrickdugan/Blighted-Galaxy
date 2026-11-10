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
    /// Mass-CDF shells sized so radial thickness is near the
    /// equal-area angular tile width.
    ApproximatelyCubicShells,
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
#[derive(Debug,Clone)]
pub struct VolumeShell {
    pub first_pair:usize,
    pub pair_count:usize,
    pub lower_radius_m:f64,
    pub upper_radius_m:f64,
    pub thickness_over_tile:f64,
}
#[derive(Debug,Clone)]
pub struct EqualVolumeSeedReport {
    pub shells:Vec<VolumeShell>,
    pub min_tile_aspect:f64,
    pub max_tile_aspect:f64,
    pub median_tile_aspect:f64,
    pub mass_relative_error:f64,
    pub center_error_m:f64,
}
/// Deterministic approximate equal-volume shells: each pair contributes
/// two equal-mass SPH parcels. Shell boundaries follow the continuum
/// enclosed-mass CDF. The number of pairs per shell is selected to make
/// shell thickness Delta r close to an angular tile edge
/// r_mid*sqrt(2*pi/k), with 2*k angular cells.
pub fn seed_approximately_cubic_shells(
    profile:&WhiteDwarfProfile,
    center:V3,velocity:V3,count:usize,
)->(Vec<Fluid3>,EqualVolumeSeedReport) {
    assert!(count>=8 && count%2==0);
    let pairs=count/2;
    let mass=profile.mass_kg/count as f64;
    let golden=PI*(3.0-5.0_f64.sqrt());
    let mut ordered=Vec::<(usize,Fluid3,Fluid3)>::with_capacity(pairs);
    let mut shells=Vec::<VolumeShell>::new();
    let mut first=0_usize;
    while first<pairs {
        let remaining=pairs-first;
        let lower=profile.radius_at_mass_fraction(
            first as f64/pairs as f64);
        let maximum=remaining.min(96);
        let mut best=(f64::INFINITY,1_usize,0.0,0.0,0.0);
        for k in 1..=maximum {
            // Prefer at least two directions per shell unless this is the
            // final residual bin.
            if k==1&&remaining>1 {continue;}
            if remaining-k==1 {continue;}
            let upper=profile.radius_at_mass_fraction(
                (first+k) as f64/pairs as f64);
            let middle=profile.radius_at_mass_fraction(
                (first as f64+0.5*k as f64)/pairs as f64);
            let thickness=upper-lower;
            let tile=middle*(2.0*PI/k as f64).sqrt();
            let ratio=thickness/tile.max(1.0);
            // The log-distance penalizes needles and flattened cells
            // symmetrically, while limiting very small outer shell sizes.
            let cost=ratio.ln().abs()+if k==2{0.05}else{0.0};
            if cost<best.0 {
                best=(cost,k,upper,middle,ratio);
            }
        }
        let (_,k,upper,_,ratio)=best;
        assert!(k>=1 && k<=remaining, "shell allocation failed");
        let index=shells.len();
        shells.push(VolumeShell {
            first_pair:first,pair_count:k,
            lower_radius_m:lower,
            upper_radius_m:upper,
            thickness_over_tile:ratio,
        });
        let phase=2.0*PI*radical_inverse(index+1,3);
        let tilt=PI*radical_inverse(index+1,5);
        let (sp,cp)=phase.sin_cos();
        let (st,ct)=tilt.sin_cos();
        for local in 0..k {
            let i=first+local;
            let f=(i as f64+0.5)/pairs as f64;
            let radius=profile.radius_at_mass_fraction(f);
            let z=1.0-2.0*(local as f64+0.5)/k as f64;
            let phi=golden*local as f64;
            let xy=(1.0-z*z).sqrt();
            let nx=xy*phi.cos();let ny=xy*phi.sin();
            let rx=ct*nx+st*z;
            let rz=-st*nx+ct*z;
            let dir=V3::new(cp*rx-sp*ny,sp*rx+cp*ny,rz);
            let d=dir*radius;
            ordered.push((i,
                Fluid3{pos_m:center+d,vel_m_s:velocity,mass_kg:mass},
                Fluid3{pos_m:center-d,vel_m_s:velocity,mass_kg:mass}));
        }
        first+=k;
    }
    ordered.sort_by_key(|(i,_,_)|*i);
    let mut particles=Vec::with_capacity(count);
    for (_,a,b) in ordered {
        particles.push(a);particles.push(b);
    }
    let total=particles.iter().map(|p|p.mass_kg).sum::<f64>();
    let mean=particles.iter().fold(V3::default(),|sum,p|
        sum+p.pos_m*p.mass_kg)/total;
    let mut aspects=shells.iter().map(|x|x.thickness_over_tile)
        .collect::<Vec<_>>();
    aspects.sort_by(|a,b|a.partial_cmp(b).unwrap());
    let report=EqualVolumeSeedReport {
        min_tile_aspect:*aspects.first().unwrap(),
        max_tile_aspect:*aspects.last().unwrap(),
        median_tile_aspect:aspects[aspects.len()/2],
        shells,
        mass_relative_error:(total-profile.mass_kg).abs()/profile.mass_kg,
        center_error_m:(mean-center).len(),
    };
    (particles,report)
}

/// Independent shell-wise angular orientations of the same hydrostatic
/// radial mass-CDF seed. Each radial shell receives a proper 3D rotation
/// from a deterministic splitmix64 stream. Antipodal pairs, center of
/// mass, velocity, mass, radial quantiles and within-shell angular geometry
/// are invariant. Inter-shell angular alignment changes meaningfully.
fn splitmix64(x:&mut u64)->u64 {
    *x=x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z=*x;
    z=(z^(z>>30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z=(z^(z>>27)).wrapping_mul(0x94d049bb133111eb);
    z^(z>>31)
}
fn random_unit(x:&mut u64)->f64 {
    ((splitmix64(x)>>11) as f64)/((1_u64<<53) as f64)
}
fn rotate_shell_vector(v:V3,cosa:f64,sina:f64,
    cosb:f64,sinb:f64,cosc:f64,sinc:f64)->V3 {
    let x1=cosa*v.x-sina*v.y;
    let y1=sina*v.x+cosa*v.y;
    let x2=cosb*x1+sinb*v.z;
    let z2=-sinb*x1+cosb*v.z;
    V3::new(cosc*x2-sinc*y1,sinc*x2+cosc*y1,z2)
}
pub fn seed_cubic_shells_oriented(
    profile:&WhiteDwarfProfile,
    center:V3,velocity:V3,count:usize,orientation_seed:u64,
)->(Vec<Fluid3>,EqualVolumeSeedReport) {
    let (mut particles,report)=seed_approximately_cubic_shells(
        profile,center,velocity,count);
    if orientation_seed==0 {return (particles,report);}
    let mut generator=orientation_seed;
    for shell in &report.shells {
        let a=2.0*PI*random_unit(&mut generator);
        let b=(1.0-2.0*random_unit(&mut generator)).acos();
        let c=2.0*PI*random_unit(&mut generator);
        let (sa,ca)=a.sin_cos();
        let (sb,cb)=b.sin_cos();
        let (sc,cc)=c.sin_cos();
        for pair in shell.first_pair..shell.first_pair+shell.pair_count {
            let offset=particles[2*pair].pos_m-center;
            let rotated=rotate_shell_vector(offset,ca,sa,cb,sb,cc,sc);
            particles[2*pair].pos_m=center+rotated;
            particles[2*pair+1].pos_m=center-rotated;
        }
    }
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
        StarSampling::ApproximatelyCubicShells=>
            seed_approximately_cubic_shells(profile,center,velocity,count).0,
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
    fn approximately_cubic_shells_cover_all_equal_mass_pairs() {
        let profile=setup();
        for n in [64,128,512,1024] {
            let (p,report)=seed_approximately_cubic_shells(
                &profile,V3::new(1e8,0.0,0.0),V3::default(),n);
            assert_eq!(p.len(),n);
            assert_eq!(report.shells.iter()
                .map(|s|s.pair_count).sum::<usize>(),n/2);
            assert!(report.shells.len()>1);
            assert!(report.min_tile_aspect>0.0&&
                report.max_tile_aspect.is_finite());
            assert!(report.mass_relative_error<1e-12);
            assert!(report.center_error_m<2e-5);
            assert!(p.iter().all(|x|x.pos_m.finite()));
            for pair in p.chunks_exact(2) {
                let pos=(pair[0].pos_m+pair[1].pos_m)*0.5;
                assert!((pos-V3::new(1e8,0.0,0.0)).len()<2e-7);
            }
        }
    }
    #[test]
    fn approximately_cubic_shells_respect_mass_cdf_quantiles() {
        let profile=setup();
        let count=128;
        let (samples,_)=seed_approximately_cubic_shells(
            &profile,V3::default(),V3::default(),count);
        for i in 0..count/2 {
            let r=samples[2*i].pos_m.len();
            let reference=profile.radius_at_mass_fraction(
                (i as f64+0.5)/(count/2) as f64);
            assert!((r-reference).abs()<1e-7);
        }
    }
    #[test]
    fn independent_shell_orientations_preserve_radius_cdf_and_center() {
        let profile=setup();
        let center=V3::new(1e8,-7e7,2e7);
        let a=seed_cubic_shells_oriented(
            &profile,center,V3::default(),256,1).0;
        let b=seed_cubic_shells_oriented(
            &profile,center,V3::default(),256,2).0;
        let c=seed_cubic_shells_oriented(
            &profile,center,V3::default(),256,1).0;
        let original=seed_cubic_shells_oriented(
            &profile,center,V3::default(),256,0).0;
        assert_eq!(a[0].pos_m,c[0].pos_m);
        assert!((a[0].pos_m-b[0].pos_m).len()>1e3);
        let mut mass=0.0;let mut cm=V3::default();
        for i in 0..256 {
            let da=(a[i].pos_m-center).len();
            let db=(b[i].pos_m-center).len();
            let dr=(original[i].pos_m-center).len();
            assert!((da-dr).abs()<1e-6);
            assert!((db-dr).abs()<1e-6);
            cm+=a[i].pos_m*a[i].mass_kg;
            mass+=a[i].mass_kg;
        }
        assert!((cm/mass-center).len()<1e-5);
        assert!((mass-profile.mass_kg).abs()/mass<1e-12);
    }
    #[test]
    fn orientation_varies_each_shell_without_breaking_antipodes() {
        let profile=setup();
        let (particles,shells)=seed_cubic_shells_oriented(
            &profile,V3::default(),V3::default(),128,987);
        assert!(shells.shells.len()>2);
        for pair in particles.chunks_exact(2) {
            assert!((pair[0].pos_m+pair[1].pos_m).len()<1e-7);
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
