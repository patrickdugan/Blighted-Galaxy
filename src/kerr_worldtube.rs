//! Local emitter-rest-frame target envelopes for Kerr optical acquisition.
 //!
 //! Given a source event recorded in equatorial Boyer-Lindquist coordinates,
 //! normalize its timelike coordinate velocity to u^mu, then project a
 //! same-coordinate-time radial/azimuthal displacement through
 //! h_{mu nu}=g_{mu nu}+u_mu*u_nu.
 //!
 //! The finite-radius spatial distance is a *local tangent-plane* estimate.
 //! It approaches Euclidean separation far from the hole and represents a
 //! physically motivated metric-aware improvement over a coordinate-circle.
 //! For targets comparable in size to their orbital radius, use a proper
 //! 4D worldtube and null geodesic boundary-value solver instead.
use crate::causal::{WorldlineHistory,WorldlineSample};
use crate::kerr::KerrGeometry;
use crate::physics::C_SI;
use crate::V2;
use std::f64::consts::{PI,TAU};

#[derive(Clone,Copy,Debug)]
pub struct KerrSourceFrame {
    pub radius_rg:f64,
    pub phi_rad:f64,
    pub radial_metric:f64,
    pub azimuth_metric:f64,
    pub cross_metric:f64,
    pub dproper_dcoordinate:f64,
}
impl KerrSourceFrame {
    pub fn from_sample(geometry:&KerrGeometry,s:WorldlineSample)->Option<Self> {
        let rg=geometry.r_g_m;
        let physical_r=s.pos_m.len();
        let r=physical_r/rg;
        if !r.is_finite()||
            r<=geometry.outer_horizon_rg()+0.04 {return None;}
        let phi=s.pos_m.y.atan2(s.pos_m.x);
        let radial=s.pos_m/physical_r;
        let azimuth=V2::new(-radial.y,radial.x);
        // Convert SI rates to dimensionless geometric coordinates:
        // dr/(c dt) and dphi/(c dt/rg).
        let urate=radial.dot(s.vel_m_s)/C_SI;
        let phi_rate=azimuth.dot(s.vel_m_s)/(r*C_SI);
        let a=geometry.spin;
        let delta=geometry.delta(r);
        let big_a=(r*r+a*a).powi(2)-a*a*delta;
        let gtt=-(1.0-2.0/r);
        let gtphi=-2.0*a/r;
        let grr=r*r/delta;
        let gphiphi=big_a/(r*r);
        // u^t = 1/sqrt(-g_{mu nu} (dx^mu/dt)(dx^nu/dt))
        let inverse_ut_sq=-(gtt+2.0*gtphi*phi_rate+
            gphiphi*phi_rate*phi_rate+grr*urate*urate);
        if !inverse_ut_sq.is_finite()||inverse_ut_sq<=0.0 {
            return None;
        }
        let dproper_dt=inverse_ut_sq.sqrt();
        let ut=1.0/dproper_dt;
        let ur_cov=grr*urate*ut;
        let uphi_cov=(gtphi+gphiphi*phi_rate)*ut;
        Some(Self {
            radius_rg:r,phi_rad:phi,
            radial_metric:grr+ur_cov*ur_cov,
            azimuth_metric:gphiphi+uphi_cov*uphi_cov,
            cross_metric:ur_cov*uphi_cov,
            dproper_dcoordinate:dproper_dt,
        })
    }
    pub fn rest_distance_squared_m2(
        &self,geometry:&KerrGeometry,point:V2,
    )->Option<f64> {
        if !point.x.is_finite()||!point.y.is_finite()||
            point.len()==0.0 {return None;}
        let r=point.len()/geometry.r_g_m;
        let d_r=r-self.radius_rg;
        let phi=point.y.atan2(point.x);
        let d_phi=(phi-self.phi_rad+PI).rem_euclid(TAU)-PI;
        let separation_geom=
            self.radial_metric*d_r*d_r+
            2.0*self.cross_metric*d_r*d_phi+
            self.azimuth_metric*d_phi*d_phi;
        if !separation_geom.is_finite() {return None;}
        Some(separation_geom.max(0.0)*geometry.r_g_m.powi(2))
    }
}

/// Approximate the earliest photon chord intersection with a finite-radius
/// emitter rest-space. Each chord carries past-directed coordinate times.
pub fn segment_hit_rest_worldtube(
    geometry:&KerrGeometry,history:&WorldlineHistory,
    photon0:V2,photon1:V2,t0_s:f64,t1_s:f64,
    envelope_m:f64,
)->Option<(f64,V2,V2)> {
    if !envelope_m.is_finite()||envelope_m<=0.0||
        t1_s>=t0_s {return None;}
    let sample_at=|f:f64| {
        let t=t0_s+(t1_s-t0_s)*f;
        let source=history.sample_at(t)?;
        let photon=photon0+(photon1-photon0)*f;
        let frame=KerrSourceFrame::from_sample(geometry,source)?;
        let squared=frame.rest_distance_squared_m2(geometry,photon)?;
        Some((squared,source.pos_m,photon))
    };
    // Euclidean closest approach is a good initial minimizer for nearby
    // timelike worldlines, with additional fractions to cover chart curvature.
    let source0=history.sample_at(t0_s)?;
    let source1=history.sample_at(t1_s)?;
    let d0=photon0-source0.pos_m;
    let dv=(photon1-photon0)-(source1.pos_m-source0.pos_m);
    let projection=if dv.len_sq()>0.0 {
        (-d0.dot(dv)/dv.len_sq()).clamp(0.0,1.0)
    } else {0.0};
    let mut fractions=[0.0,0.25,0.5,0.75,1.0,projection];
    fractions.sort_by(|a,b|a.partial_cmp(b).unwrap());
    let mut minimum=None;
    for f in fractions {
        if let Some((d2,source,ray))=sample_at(f) {
            if minimum.as_ref().is_none_or(|(dist,_,_,_)|d2<*dist) {
                minimum=Some((d2,f,source,ray));
            }
        }
    }
    let (_,f_min,src,ray)=minimum?;
    if minimum?.0>envelope_m*envelope_m {return None;}
    // Find the first entry into the target volume along the photon chord.
    // If already inside at the beginning, f=0 is correct.
    if sample_at(0.0)?.0<=envelope_m*envelope_m {
        return Some((0.0,source0.pos_m,photon0));
    }
    let mut lo=0.0;let mut hi=f_min;
    if hi>0.0 {
        for _ in 0..36 {
            let mid=0.5*(lo+hi);
            if sample_at(mid)?.0<=envelope_m*envelope_m {
                hi=mid;
            } else {lo=mid;}
        }
        let (_,s,p)=sample_at(hi)?;
        Some((hi,s,p))
    } else {Some((f_min,src,ray))}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    fn geometry()->KerrGeometry {
        KerrGeometry::new(10.0*SOLAR_MASS_KG,0.85).unwrap()
    }
    #[test]
    fn stationary_flat_far_limit_has_euclidean_distance() {
        let g=geometry();
        let radius=2000.0*g.r_g_m;
        let s=WorldlineSample{
            time_s:0.0,pos_m:V2::new(radius,0.0),
            vel_m_s:V2::default(),
        };
        let f=KerrSourceFrame::from_sample(&g,s).unwrap();
        let r=5000.0;
        let d=f.rest_distance_squared_m2(
            &g,V2::new(radius+r,0.0)).unwrap().sqrt();
        assert!((d-r).abs()/r<0.002);
    }
    #[test]
    fn local_frame_circular_orbit_is_timelike_and_distance_positive() {
        let g=geometry();
        let r=20.0;
        let orbit=g.circular_timelike(r,1).unwrap();
        let x=V2::new(r*g.r_g_m,0.0);
        let vel=V2::new(0.0,r*g.r_g_m*orbit.omega_rad_s);
        let frame=KerrSourceFrame::from_sample(&g,WorldlineSample{
            time_s:0.0,pos_m:x,vel_m_s:vel
        }).unwrap();
        assert!(frame.dproper_dcoordinate>0.0&&
            frame.dproper_dcoordinate<1.0);
        assert!(frame.azimuth_metric>r*r);
        let dr=V2::new(x.x+0.7*g.r_g_m,0.0);
        let d=frame.rest_distance_squared_m2(&g,dr).unwrap().sqrt();
        assert!(d>0.7*g.r_g_m);
    }
    #[test]
    fn metric_rejects_non_timelike_coordinate_velocity() {
        let g=geometry();
        let p=V2::new(20.0*g.r_g_m,0.0);
        let frame=KerrSourceFrame::from_sample(&g,WorldlineSample{
            time_s:0.0,pos_m:p,vel_m_s:V2::new(2.0*C_SI,0.0),
        });
        assert!(frame.is_none());
    }
    #[test]
    fn local_worldtube_detects_stationary_chord() {
        let g=geometry();
        let x=V2::new(30.0*g.r_g_m,0.0);
        let mut h=WorldlineHistory::new(3);
        h.push(WorldlineSample{
            time_s:0.0,pos_m:x,vel_m_s:V2::default()
        }).unwrap();
        h.push(WorldlineSample{
            time_s:1.0,pos_m:x,vel_m_s:V2::default()
        }).unwrap();
        let p0=x+V2::new(2.0*g.r_g_m,0.0);
        let p1=x-V2::new(2.0*g.r_g_m,0.0);
        let hit=segment_hit_rest_worldtube(&g,&h,p0,p1,
            0.9,0.8,0.9*g.r_g_m).unwrap();
        assert!(hit.0>0.0&&hit.0<0.5);
    }
}
