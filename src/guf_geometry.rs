//! BG-STR-05: finite-angle Guf relay mounts and beam footprints.
//! Angles and steering delays are scenario engineering parameters.
//! Light travels at one light-year per galactic year.
use crate::V2;
use serde::{Serialize,Deserialize};
use std::f64::consts::{PI,TAU};

pub const SLEW_RATE_RAD_PER_YEAR:f64=0.00025;
pub const BEAM_HALF_ANGLE_RAD:f64=0.0008;
pub const TARGET_RADIUS_LY:f64=8.0;

pub fn signed_angle_delta(to:f64,from:f64)->f64 {
    (to-from+PI).rem_euclid(TAU)-PI
}
pub fn bearing(source:V2,target:V2)->f64 {
    let delta=target-source;
    delta.y.atan2(delta.x)
}
pub fn beam_intersects(source:V2,target:V2,aim:f64)->bool{
    let range=(target-source).len();
    if !range.is_finite()||range<1.0{return false;}
    let target_half=(TARGET_RADIUS_LY/range).atan();
    signed_angle_delta(bearing(source,target),aim).abs()
        <=BEAM_HALF_ANGLE_RAD+target_half
}
#[derive(Clone,Copy,Debug,Serialize,Deserialize)]
pub struct SlewMount{
    pub initial_azimuth_rad:f64,
    pub target_azimuth_rad:f64,
    pub slew_start_year:f64,
    pub slew_end_year:f64,
    pub assigned_target:Option<usize>,
}
impl Default for SlewMount {
    fn default()->Self{
        Self{initial_azimuth_rad:0.0,target_azimuth_rad:0.0,
            slew_start_year:0.0,slew_end_year:0.0,assigned_target:None}
    }
}
impl SlewMount {
    pub fn bearing_at(&self,year:f64)->f64{
        if year<=self.slew_start_year{return self.initial_azimuth_rad;}
        if year>=self.slew_end_year{return self.target_azimuth_rad;}
        let elapsed=(year-self.slew_start_year)/
            (self.slew_end_year-self.slew_start_year).max(1e-10);
        self.initial_azimuth_rad+signed_angle_delta(
            self.target_azimuth_rad,self.initial_azimuth_rad)*elapsed
    }
    /// Reorientation starts once a light-speed command arrives from a relay
    /// controller. Its initial bearing matches the already-settled mount.
    pub fn schedule(&self,year:f64,command_received:f64,
        requested_azimuth:f64,target:usize)->Option<Self>{
        if !year.is_finite()||!command_received.is_finite()||
            !requested_azimuth.is_finite()||command_received<year||
            self.slew_end_year>year{return None;}
        let from=self.bearing_at(year);
        let angle=signed_angle_delta(requested_azimuth,from);
        let finish=command_received+angle.abs()/SLEW_RATE_RAD_PER_YEAR;
        Some(Self{initial_azimuth_rad:from,
            target_azimuth_rad:from+angle,
            slew_start_year:command_received,
            slew_end_year:finish,
            assigned_target:Some(target)})
    }
    pub fn aligned(&self,year:f64,target:usize,desired_azimuth:f64)->bool{
        self.assigned_target==Some(target)&&year>=self.slew_end_year
            &&signed_angle_delta(self.bearing_at(year),desired_azimuth).abs()
                <=BEAM_HALF_ANGLE_RAD
    }
    pub fn valid(&self,year:f64,site_count:usize)->bool{
        self.initial_azimuth_rad.is_finite()&&self.target_azimuth_rad.is_finite()
            &&self.initial_azimuth_rad.abs()<=20.0*TAU
            &&self.target_azimuth_rad.abs()<=20.0*TAU
            &&self.slew_start_year.is_finite()&&self.slew_end_year.is_finite()
            &&self.slew_start_year>=0.0&&self.slew_end_year>=self.slew_start_year
            &&self.slew_end_year<=200_000.0
            &&year>=0.0&&self.assigned_target.is_none_or(|t|t<site_count)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn shortest_way_round_azimuth_wrap(){
        let from=179.0_f64.to_radians();
        let to=-179.0_f64.to_radians();
        assert!((signed_angle_delta(to,from)-2.0_f64.to_radians()).abs()<1e-12);
    }
    #[test] fn angular_steering_takes_finite_time_and_can_hold_prior_bearing(){
        let m=SlewMount::default().schedule(
            100_000.0,105_000.0,PI,2).unwrap();
        assert!(m.slew_end_year>117_000.0);
        assert_eq!(m.bearing_at(100_001.0),0.0);
        assert!(!m.aligned(106_000.0,2,PI));
        assert!(m.aligned(m.slew_end_year,2,PI));
        assert!(m.schedule(106_000.0,106_001.0,0.0,1).is_none());
    }
    #[test] fn displaced_star_leaves_the_locked_beam_cone(){
        let from=V2::new(0.0,0.0);
        let target=V2::new(4000.0,0.0);
        assert!(beam_intersects(from,target,0.0));
        assert!(!beam_intersects(from,target+V2::new(0.0,100.0),0.0));
    }
}
