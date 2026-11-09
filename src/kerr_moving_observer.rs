//! Causal Kerr optical acquisition with a moving equatorial observer.
 //!
 //! Observer motion is a circular equatorial timelike Kerr test geodesic
 //! (or a local ZAMO reference). Incoming camera directions are aberrated
 //! from the moving observer's local tetrad to the ZAMO launch frame.
 //! Axial Killing symmetry rotates time-labelled photon paths to the
 //! observer's *current* Boyer-Lindquist azimuth. The emitter can be any
 //! recorded, subluminal piecewise-linear equatorial source history.
 //!
 //! The physical limitations remain explicit: the finite target radius
 //! is a Euclidean proxy in BL chart x/y, and only the equatorial subset
 //! of Kerr's null geodesics is used to intersect source worldlines.
 //! Orbital sources can be checked exactly at the camera's observation
 //! event; realistic 3D target volumes require a 4D spacetime worldtube.
use crate::causal::WorldlineHistory;
use crate::kerr_camera::{
    EquatorialNullPath,KerrCamera,NullPathEnd,
};
use crate::kerr_optics::OrbitingCombatant;
use crate::physics::C_SI;
use crate::V2;
use crate::kerr_worldtube::segment_hit_rest_worldtube;
use std::f64::consts::{FRAC_PI_2,TAU};

#[derive(Clone,Copy,Debug)]
pub enum ObserverMotion {
    Zamo,
    CircularGeodesic {direction:i8},
}
#[derive(Clone,Debug)]
pub struct KerrMovingObserver {
    pub camera:KerrCamera,
    pub motion:ObserverMotion,
    pub omega_rad_s:f64,
    pub phi_at_zero_rad:f64,
    /// Velocity relative to the local ZAMO along phi; |beta|<1.
    pub beta_zamo:f64,
    pub lapse:f64,
    pub gamma_zamo:f64,
}
impl KerrMovingObserver {
    pub fn new(
        camera:KerrCamera,motion:ObserverMotion,phi_at_zero_rad:f64,
    )->Option<Self> {
        if (camera.inclination_rad-FRAC_PI_2).abs()>1e-8||
            !phi_at_zero_rad.is_finite() {return None;}
        let r=camera.observer_rg;
        let a=camera.geometry.spin;
        let delta=camera.geometry.delta(r);
        let big_a=(r*r+a*a).powi(2)-a*a*delta;
        let lapse=(r*r*delta/big_a).sqrt();
        let gphiphi=big_a/(r*r);
        let omega_zamo=camera.geometry.zamo_omega_rad_s(r)?;
        let omega=match motion {
            ObserverMotion::Zamo=>omega_zamo,
            ObserverMotion::CircularGeodesic{direction}=>
                camera.geometry.circular_timelike(r,direction)?.omega_rad_s,
        };
        let beta=(omega-omega_zamo)*camera.geometry.r_g_m/C_SI*
            gphiphi.sqrt()/lapse;
        if !beta.is_finite()||beta.abs()>=1.0 {return None;}
        let gamma=1.0/(1.0-beta*beta).sqrt();
        Some(Self {camera,motion,
            omega_rad_s:omega,phi_at_zero_rad,beta_zamo:beta,
            lapse,gamma_zamo:gamma})
    }
    pub fn phi_at(&self,coordinate_time_s:f64)->f64 {
        self.phi_at_zero_rad+self.omega_rad_s*coordinate_time_s
    }
    pub fn position_m(&self,coordinate_time_s:f64)->V2 {
        V2::from_angle(self.phi_at(coordinate_time_s))*
            self.camera.observer_rg*self.camera.geometry.r_g_m
    }
    /// Proper time per unit Boyer-Lindquist time for this observer.
    pub fn dproper_dcoordinate(&self)->f64 {
        self.lapse/self.gamma_zamo
    }
    /// Special-relativistic aberration from a moving camera's local
    /// past-directed photon direction to the ZAMO's local direction.
    pub fn observer_alpha_to_zamo_alpha(&self,alpha_rg:f64)->Option<f64> {
        if !alpha_rg.is_finite() {return None;}
        let r=self.camera.observer_rg;
        let n=alpha_rg/(r*r+alpha_rg*alpha_rg).sqrt();
        let denom=1.0-self.beta_zamo*n;
        if denom<=0.0 {return None;}
        let n_zamo=(n-self.beta_zamo)/denom;
        if !n_zamo.is_finite()||n_zamo.abs()>=1.0 {return None;}
        Some(r*n_zamo/(1.0-n_zamo*n_zamo).sqrt())
    }
    /// Kerr axial symmetry transports a ray launched at phi=0 to
    /// the moving observer event's instantaneous azimuth.
    pub fn ray_at(
        &self,observer_time_s:f64,alpha_camera_rg:f64,
        inner_stop_rg:f64,max_steps:usize,
    )->Option<EquatorialNullPath> {
        if !observer_time_s.is_finite() {return None;}
        let alpha_zamo=self.observer_alpha_to_zamo_alpha(alpha_camera_rg)?;
        let mut ray=self.camera.equatorial_null_path(
            alpha_zamo,inner_stop_rg,max_steps)?;
        let phi=self.phi_at(observer_time_s);
        let (sin,cos)=phi.sin_cos();
        for sample in &mut ray.samples {
            let p=sample.chart_xy_m;
            sample.chart_xy_m=V2::new(
                cos*p.x-sin*p.y,
                sin*p.x+cos*p.y
            );
            sample.phi_rad+=phi;
        }
        // Preserve the observer-camera alpha in the source ray metadata;
        // its corresponding ZAMO direction is obtained through the method.
        ray.impact_alpha_rg=alpha_camera_rg;
        Some(ray)
    }
}
#[derive(Debug)]
pub struct MovingOpticalObservation {
    pub observer_time_s:f64,
    pub observer_pos_m:V2,
    pub emitter_time_s:f64,
    pub delay_s:f64,
    pub apparent_pos_m:V2,
    pub current_pos_m:V2,
    pub intersection_chart_m:V2,
    pub alpha_camera_rg:f64,
    pub ray:EquatorialNullPath,
}
/// For every sampled observer-local camera direction, find intersections
/// between the time-labelled Kerr null path and an *arbitrary recorded*
/// source history, using BL chart x/y as the finite target envelope proxy.
/// Return the shortest-delay detected path from the sampled ray bundle.
pub fn observe_recorded_source(
    observer:&KerrMovingObserver,
    source_history:&WorldlineHistory,
    observer_time_s:f64,
    source_envelope_rg:f64,
    alpha_halfspan_rg:f64,
    ray_count:usize,
    inner_stop_rg:f64,
)->Option<MovingOpticalObservation> {
    if !observer_time_s.is_finite()||source_envelope_rg<=0.0||
        !source_envelope_rg.is_finite()||
        !alpha_halfspan_rg.is_finite()||alpha_halfspan_rg<0.0||
        ray_count==0||ray_count>1025 {return None;}
    let observer_pos=observer.position_m(observer_time_s);
    let now=source_history.sample_at(observer_time_s)?.pos_m;
    let envelope_m=source_envelope_rg*observer.camera.geometry.r_g_m;
    let mut best:Option<MovingOpticalObservation>=None;
    for i in 0..ray_count {
        let alpha=if ray_count==1 {0.0} else {
            -alpha_halfspan_rg+
                2.0*alpha_halfspan_rg*i as f64/(ray_count-1) as f64
        };
        let ray=observer.ray_at(observer_time_s,alpha,
            inner_stop_rg,3600)?;
        if ray.end==NullPathEnd::HorizonChartCutoff {continue;}
        let mut matched=None;
        for seg in ray.samples.windows(2) {
            let (p0,p1)=(seg[0],seg[1]);
            let (t0,t1)=(
                observer_time_s+p0.t_offset_s,
                observer_time_s+p1.t_offset_s
            );
            let e0=match source_history.sample_at(t0) {
                Some(v)=>v,None=>continue,
            };
            let e1=match source_history.sample_at(t1) {
                Some(v)=>v,None=>continue,
            };
            // The source's rest-space metric sets its finite collision
            // envelope. The photon/source chord remains a local linear
            // approximation between consecutive Kerr ray samples.
            if let Some((f,observed,hit))=segment_hit_rest_worldtube(
                &observer.camera.geometry,source_history,
                p0.chart_xy_m,p1.chart_xy_m,t0,t1,envelope_m,
            ) {
                let time=t0+(t1-t0)*f;
                matched=Some((time,observed,hit));
                break;
            }
        }
        if let Some((emitter_time,observed,hit))=matched {
            let delay=observer_time_s-emitter_time;
            if best.as_ref().is_none_or(|v|delay<v.delay_s) {
                best=Some(MovingOpticalObservation{
                    observer_time_s,observer_pos_m:observer_pos,
                    emitter_time_s:emitter_time,delay_s:delay,
                    apparent_pos_m:observed,current_pos_m:now,
                    intersection_chart_m:hit,alpha_camera_rg:alpha,ray,
                });
            }
        }
    }
    best
}
/// Place a circular target at the phase where a selected moving-observer
/// past-directed ray crosses its orbital radius. This is a deterministic
/// encounter generator, with the target's subsequent states recorded as
/// an independent timelike history.
pub fn source_aligned_to_observer_ray(
    observer:&KerrMovingObserver,
    local_alpha_rg:f64,source_radius_rg:f64,at_time_s:f64,
    source_direction:i8,envelope_rg:f64,
)->Option<OrbitingCombatant> {
    let path=observer.ray_at(at_time_s,local_alpha_rg,
        source_radius_rg-2.0,3600)?;
    let segment=path.samples.windows(2).find(|p|
        p[0].radius_rg>=source_radius_rg&&
        p[1].radius_rg<=source_radius_rg)?;
    let f=(segment[0].radius_rg-source_radius_rg)/
        (segment[0].radius_rg-segment[1].radius_rg);
    let ph=segment[0].phi_rad+
        (segment[1].phi_rad-segment[0].phi_rad)*f;
    let ts=at_time_s+segment[0].t_offset_s+
        (segment[1].t_offset_s-segment[0].t_offset_s)*f;
    let orbit=observer.camera.geometry.circular_timelike(
        source_radius_rg,source_direction)?;
    let phase=(ph-orbit.omega_rad_s*ts).rem_euclid(TAU);
    OrbitingCombatant::new(&observer.camera,
        source_radius_rg,phase,source_direction,envelope_rg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kerr::KerrGeometry;
    use crate::physics::SOLAR_MASS_KG;
    fn cam(spin:f64)->KerrCamera {
        KerrCamera::new(
            KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap(),
            85.0,FRAC_PI_2,36.0).unwrap()
    }
    #[test]
    fn zamo_reference_has_no_local_aberration() {
        let observer=KerrMovingObserver::new(
            cam(0.8),ObserverMotion::Zamo,0.0).unwrap();
        assert!(observer.beta_zamo.abs()<1e-12);
        for alpha in [-8.0,0.0,5.0] {
            assert!((observer.observer_alpha_to_zamo_alpha(alpha).unwrap()-
                alpha).abs()<1e-10);
        }
        let a=observer.ray_at(0.0,3.0,17.0,3600).unwrap();
        let b=observer.camera.equatorial_null_path(
            3.0,17.0,3600).unwrap();
        assert!((a.samples.last().unwrap().phi_rad-
            b.samples.last().unwrap().phi_rad).abs()<1e-12);
    }
    #[test]
    fn orbiting_observer_is_timelike_and_aberrates_rays() {
        let observer=KerrMovingObserver::new(
            cam(0.85),ObserverMotion::CircularGeodesic{direction:1},
            0.4).unwrap();
        assert!(observer.beta_zamo>0.0&&observer.beta_zamo<1.0);
        assert!(observer.dproper_dcoordinate()>0.0 &&
            observer.dproper_dcoordinate()<1.0);
        assert!(observer.observer_alpha_to_zamo_alpha(0.0).unwrap()<0.0);
        let t=0.08;
        let p=observer.position_m(t);
        let ray=observer.ray_at(t,0.0,17.0,3600).unwrap();
        assert!((ray.samples[0].chart_xy_m-p).len()<
            1e-5*observer.camera.geometry.r_g_m);
        assert!(ray.samples.windows(2).all(|w|
            w[1].t_offset_s<w[0].t_offset_s));
    }
    #[test]
    fn moving_observer_gets_delayed_optical_detection() {
        for spin in [0.0,0.75,-0.75] {
            let ob=KerrMovingObserver::new(
                cam(spin),ObserverMotion::CircularGeodesic{direction:1},
                0.27).unwrap();
            let now=0.06;
            let source=source_aligned_to_observer_ray(
                &ob,0.0,20.0,now,1,1.0).unwrap();
            let history=source.record(&ob.camera,
                now-0.015,now+0.008,0.00008).unwrap();
            let observed=observe_recorded_source(&ob,&history,
                now,1.0,4.0,65,17.0).unwrap();
            assert!(observed.delay_s>0.0&&observed.delay_s<0.015);
            assert!((observed.apparent_pos_m-observed.current_pos_m)
                .len()>ob.camera.geometry.r_g_m);
            assert!((observed.intersection_chart_m-
                observed.apparent_pos_m).len()<
                1.01*ob.camera.geometry.r_g_m);
        }
    }
    #[test]
    fn insufficient_past_history_is_unobservable() {
        let ob=KerrMovingObserver::new(
            cam(0.0),ObserverMotion::CircularGeodesic{direction:1},
            0.0).unwrap();
        let now=0.06;
        let source=source_aligned_to_observer_ray(
            &ob,0.0,20.0,now,1,1.0).unwrap();
        let h=source.record(&ob.camera,
            now-0.0003,now+0.001,0.0001).unwrap();
        assert!(observe_recorded_source(&ob,&h,now,
            1.0,0.0,1,17.0).is_none());
    }
    #[test]
    fn generic_accelerating_source_history_is_accepted() {
        let ob=KerrMovingObserver::new(
            cam(0.65),ObserverMotion::CircularGeodesic{direction:1},
            0.25).unwrap();
        let now=0.06;
        let star=source_aligned_to_observer_ray(
            &ob,0.0,20.0,now,1,1.0).unwrap();
        let mut history=WorldlineHistory::new(300);
        for i in 0..=250 {
            let t=now-0.012+i as f64*0.00008;
            let mut p=star.position_m(&ob.camera,t);
            // Smooth, low-amplitude nongeodesic perturbation. The recorded
            // source need not be a strictly circular orbit.
            let d=0.015*ob.camera.geometry.r_g_m;
            p+=V2::new(
                d*(200.0*t).sin(),d*(150.0*t).cos());
            let v=star.velocity_m_s(&ob.camera,t)+V2::new(
                d*200.0*(200.0*t).cos(),
                -d*150.0*(150.0*t).sin());
            history.push(crate::causal::WorldlineSample{
                time_s:t,pos_m:p,vel_m_s:v
            }).unwrap();
        }
        assert!(history.validate_segments());
        let found=observe_recorded_source(
            &ob,&history,now,1.0,4.0,65,17.0).unwrap();
        assert!(found.emitter_time_s<now);
    }
}
