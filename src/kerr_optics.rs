//! Curved-spacetime optical acquisition of a *moving* posthuman body.
 //!
 //! A ZAMO equatorial Kerr camera traces photon null geodesics backwards in
 //! Boyer-Lindquist coordinate time, then tests the resulting worldline
 //! against recorded circular-orbit combatant trajectories.
 //!
 //! The null geodesic is physical for a test photon in stationary Kerr.
 //! The craft is a timelike test body on a circular equatorial Kerr orbit.
 //! Intersection uses a finite coordinate-envelope radius in the shared
 //! BL chart, a gameplay approximation to a locally measured worldtube.
 //! This is sensor/retarded-vision integration, not outgoing laser fire.
use crate::causal::{WorldlineHistory,WorldlineSample};
use crate::kerr_camera::{EquatorialNullPath,KerrCamera,NullPathEnd};
use crate::physics::C_SI;
use crate::{swept_circle_hit,V2};
use std::f64::consts::{FRAC_PI_2,TAU};

#[derive(Debug,Clone,Copy)]
pub struct OrbitingCombatant {
    pub orbital_radius_rg:f64,
    pub orbit_rate_rad_s:f64,
    pub phase_at_zero_rad:f64,
    /// Collision envelope, *coordinate radius* in gravitational radii.
    pub envelope_rg:f64,
}
impl OrbitingCombatant {
    pub fn new(
        camera:&KerrCamera,orbital_radius_rg:f64,phase_at_zero_rad:f64,
        direction:i8,envelope_rg:f64,
    )->Option<Self> {
        if (camera.inclination_rad-FRAC_PI_2).abs()>1e-8 ||
            envelope_rg<=0.0||!envelope_rg.is_finite()||
            orbital_radius_rg <= camera.geometry.isco_rg(direction)? + 1.0 ||
            orbital_radius_rg >= camera.observer_rg-5.0 {return None;}
        let orbit=camera.geometry.circular_timelike(orbital_radius_rg,direction)?;
        Some(Self{
            orbital_radius_rg,orbit_rate_rad_s:orbit.omega_rad_s,
            phase_at_zero_rad,envelope_rg,
        })
    }
    pub fn position_m(&self,camera:&KerrCamera,time_s:f64)->V2 {
        let theta=self.phase_at_zero_rad+self.orbit_rate_rad_s*time_s;
        V2::from_angle(theta)*self.orbital_radius_rg*camera.geometry.r_g_m
    }
    pub fn velocity_m_s(&self,camera:&KerrCamera,time_s:f64)->V2 {
        let theta=self.phase_at_zero_rad+self.orbit_rate_rad_s*time_s;
        let tangential=V2::new(-theta.sin(),theta.cos());
        tangential*(self.orbital_radius_rg*camera.geometry.r_g_m*self.orbit_rate_rad_s)
    }
    pub fn record(
        &self,camera:&KerrCamera,start_s:f64,end_s:f64,step_s:f64,
    )->Option<WorldlineHistory> {
        if !start_s.is_finite()||!end_s.is_finite()||end_s<=start_s||
            !step_s.is_finite()||step_s<=0.0 {return None;}
        let n=((end_s-start_s)/step_s).ceil() as usize;
        if n<1||n>8190 {return None;}
        let mut history=WorldlineHistory::new(n+2);
        for i in 0..=n {
            let t=start_s+(end_s-start_s)*i as f64/n as f64;
            history.push(WorldlineSample{
                time_s:t,
                pos_m:self.position_m(camera,t),
                vel_m_s:self.velocity_m_s(camera,t),
            }).ok()?;
        }
        Some(history)
    }
}

/// A ray reconstructed from the observer's past light cone that intersects
/// a stored source worldline. The source's present position is shown
/// separately, using that *same historical record* at observation time.
#[derive(Debug)]
pub struct OpticalAcquisition {
    pub alpha_rg:f64,
    pub observer_time_s:f64,
    pub emission_time_s:f64,
    pub light_delay_s:f64,
    pub observed_pos_m:V2,
    pub present_pos_m:V2,
    pub photon_chart_pos_m:V2,
    pub radial_error:f64,
    pub polar_error:f64,
    pub ray:EquatorialNullPath,
}

pub fn observe_combatant(
    camera:&KerrCamera,
    source_history:&WorldlineHistory,
    observed_time_s:f64,
    source_envelope_rg:f64,
    alpha_span_rg:f64,
    ray_count:usize,
    inner_stop_rg:f64,
)->Option<OpticalAcquisition> {
    if !observed_time_s.is_finite()||source_envelope_rg<=0.0||
        !source_envelope_rg.is_finite()||alpha_span_rg<0.0||
        !alpha_span_rg.is_finite()||ray_count<1||ray_count>1025 {
        return None;
    }
    let present_pos=source_history.sample_at(observed_time_s)?.pos_m;
    let radius=source_envelope_rg*camera.geometry.r_g_m;
    let mut result:Option<OpticalAcquisition>=None;
    for i in 0..ray_count {
        let alpha=if ray_count==1 {0.0} else {
            -alpha_span_rg+2.0*alpha_span_rg*i as f64/
                (ray_count-1) as f64
        };
        let path=camera.equatorial_null_path(alpha,inner_stop_rg,3600)?;
        if path.end==NullPathEnd::HorizonChartCutoff||
            path.samples.len()<2 {continue;}
        let mut observation=None;
        for segment in path.samples.windows(2) {
            let photon_prev=segment[0];
            let photon_next=segment[1];
            let t0=observed_time_s+photon_prev.t_offset_s;
            let t1=observed_time_s+photon_next.t_offset_s;
            let src0=match source_history.sample_at(t0) {
                Some(p)=>p,None=>continue,
            };
            let src1=match source_history.sample_at(t1) {
                Some(p)=>p,None=>continue,
            };
            if let Some(frac)=swept_circle_hit(
                photon_prev.chart_xy_m,photon_next.chart_xy_m,0.0,
                src0.pos_m,src1.pos_m,radius,
            ) {
                let t_emit=t0+(t1-t0)*frac;
                let source=source_history.sample_at(t_emit)?;
                let ray_pos=photon_prev.chart_xy_m+
                    (photon_next.chart_xy_m-photon_prev.chart_xy_m)*frac;
                observation=Some((t_emit,source.pos_m,ray_pos));
                break;
            }
        }
        if let Some((t_emit,observed,ray_pos))=observation {
            let delay=observed_time_s-t_emit;
            if result.as_ref().is_none_or(|old|delay<old.light_delay_s) {
                result=Some(OpticalAcquisition{
                    alpha_rg:alpha,
                    observer_time_s:observed_time_s,
                    emission_time_s:t_emit,
                    light_delay_s:delay,
                    observed_pos_m:observed,
                    present_pos_m:present_pos,
                    photon_chart_pos_m:ray_pos,
                    radial_error:path.max_radial_residual,
                    polar_error:path.max_polar_residual,
                    ray:path,
                });
            }
        }
    }
    result
}

/// Deterministically choose a circular source phase that aligns its earlier
/// position with the past-directed null ray at a selected orbital radius.
/// Useful for reproducible game/eval encounters. No future observation is
/// needed to construct the source's subsequent worldline.
pub fn source_aligned_to_past_ray(
    camera:&KerrCamera,
    alpha_rg:f64,
    source_radius_rg:f64,
    observer_time_s:f64,
    direction:i8,
    source_envelope_rg:f64,
)->Option<OrbitingCombatant> {
    let path=camera.equatorial_null_path(
        alpha_rg,source_radius_rg-2.0,3600)?;
    let pair=path.samples.windows(2).find(|p|
        p[0].radius_rg>=source_radius_rg &&
        p[1].radius_rg<=source_radius_rg)?;
    let frac=((pair[0].radius_rg-source_radius_rg)/
        (pair[0].radius_rg-pair[1].radius_rg)).clamp(0.0,1.0);
    let phi=pair[0].phi_rad+
        (pair[1].phi_rad-pair[0].phi_rad)*frac;
    let event_offset=pair[0].t_offset_s+
        (pair[1].t_offset_s-pair[0].t_offset_s)*frac;
    let orbit=camera.geometry.circular_timelike(source_radius_rg,direction)?;
    let emitter_time=observer_time_s+event_offset;
    let phase=(phi-orbit.omega_rad_s*emitter_time).rem_euclid(TAU);
    OrbitingCombatant::new(camera,source_radius_rg,phase,
        direction,source_envelope_rg)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::kerr::KerrGeometry;
    use crate::physics::SOLAR_MASS_KG;
    fn cam(spin:f64)->KerrCamera {
        KerrCamera::new(
            KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap(),
            85.0,FRAC_PI_2,35.0).unwrap()
    }
    #[test]
    fn equatorial_kerr_path_carries_strictly_past_directed_time() {
        let camera=cam(0.75);
        let ray=camera.equatorial_null_path(1.5,16.0,3600).unwrap();
        assert!(ray.samples.len()>10);
        assert!(ray.samples.windows(2).all(|w|
            w[1].t_offset_s<w[0].t_offset_s));
        assert!(ray.max_radial_residual<1e-4);
        assert!(ray.max_polar_residual<1e-4);
    }
    #[test]
    fn schwarzschild_delay_exceeds_flat_light_distance() {
        let camera=cam(0.0);
        let path=camera.equatorial_null_path(0.0,18.0,3600).unwrap();
        let pair=path.samples.windows(2).find(|w|
            w[0].radius_rg>=20.0&&w[1].radius_rg<=20.0).unwrap();
        let f=(pair[0].radius_rg-20.0)/
            (pair[0].radius_rg-pair[1].radius_rg);
        let emit_offset=pair[0].t_offset_s+
            f*(pair[1].t_offset_s-pair[0].t_offset_s);
        let geometric_delay=(85.0-20.0)*camera.geometry.r_g_m/C_SI;
        assert!(-emit_offset>geometric_delay);
        let analytical=(65.0+2.0*(83.0/18.0f64).ln())*
            camera.geometry.r_g_m/C_SI;
        let rel=(-emit_offset-analytical).abs()/analytical;
        assert!(rel<0.002,"relative Kerr null time error {rel}");
    }
    #[test]
    fn curved_photon_finds_moving_combatant_in_recorded_past() {
        let camera=cam(0.75);
        let now=0.06;
        let source=source_aligned_to_past_ray(
            &camera,0.0,20.0,now,1,1.0).unwrap();
        let past=source.record(&camera,0.0,0.08,0.0001).unwrap();
        assert!(past.validate_segments());
        let hit=observe_combatant(
            &camera,&past,now,source.envelope_rg,
            5.0,81,17.0).expect("optical image of moving circular combatant");
        assert!(hit.light_delay_s>0.0 && hit.light_delay_s<0.04);
        assert!((hit.observed_pos_m-hit.present_pos_m).len()>
            camera.geometry.r_g_m*1.0);
        assert!(hit.radial_error<1e-4&&hit.polar_error<1e-4);
        assert!((hit.observed_pos_m-hit.photon_chart_pos_m).len()<=
            source.envelope_rg*camera.geometry.r_g_m*1.001);
    }
    #[test]
    fn finite_history_prevents_future_information_leak() {
        let camera=cam(0.0);
        let source=source_aligned_to_past_ray(
            &camera,0.0,20.0,0.06,1,1.0).unwrap();
        let incomplete=source.record(&camera,0.059,0.08,0.0001).unwrap();
        assert!(observe_combatant(
            &camera,&incomplete,0.06,1.0,0.0,1,17.0).is_none());
    }
    #[test]
    fn spin_changes_curved_path_and_time_delay() {
        let zero=cam(0.0);
        let spin=cam(0.85);
        let s0=zero.equatorial_null_path(4.0,20.0,3600).unwrap();
        let s1=spin.equatorial_null_path(4.0,20.0,3600).unwrap();
        let p0=s0.samples.last().unwrap();
        let p1=s1.samples.last().unwrap();
        assert!((p0.phi_rad-p1.phi_rad).abs()>1e-5);
        assert!((p0.t_offset_s-p1.t_offset_s).abs()>1e-8);
    }
}
