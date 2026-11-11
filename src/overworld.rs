//! Procedural, deterministic galactic overworld model.
//!
//! The 100,000-year game scenario is bounded by a 0.10c effective
//! Jenkins Swarm wavefront: every blighted/contested site is within
//! 10,000 light-years of the Orion seed at T+100,000 years.
//!
//! Resistance creates lagged *inward* notches and narrow outward
//! corridors inside that causal envelope. Time, path length, and
//! proper time use light-years / years (c=1 ly/year).
//!
//! Canon: the 22-episode macro-history centers the Blight in the
//! Orion Arm and follows Amicus, Perseus/Olympus, Outer Arm refugees,
//! Bel-Air, the Nebular Wastes, and diverged Patricia Trie branches.
//! Exact system coordinates and dates are not in the manuscripts:
//! these waypoints are explicitly schematic gameplay placements.
use crate::V2;
use std::f64::consts::{PI,TAU};

pub const DISK_RADIUS_LY:f64=52_000.0;
pub const MAP_EXTENT_LY:f64=59_000.0;
pub const SOL:V2=V2::new(0.0,-26_000.0);
pub const ORION_OUTBREAK:V2=V2::new(-1_100.0,-25_900.0);
pub const FRONT_BETA:f64=0.10;
pub const TIMELINE_START_YR:f64=0.0;
pub const TIMELINE_PRESENT_YR:f64=100_000.0;
pub const TIMELINE_MAX_YR:f64=200_000.0;
pub const FRONT_AT_100_KYR_LY:f64=10_000.0;
pub const GALACTIC_SECTOR_DEG:f64=30.0;
pub const GALACTIC_RING_LY:f64=15_000.0;
pub const LOCAL_SECTOR_DEG:f64=15.0;
pub const LOCAL_RING_LY:f64=2_000.0;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum SectorStatus {
    Unaffected,
    Held,
    Contested,
    Defended,
    Reclaimed,
}
impl SectorStatus {
    pub fn name(self)->&'static str {
        match self{
            Self::Unaffected=>"UNREACHED",
            Self::Held=>"BLIGHT-HELD",
            Self::Contested=>"CONTESTED",
            Self::Defended=>"DEFENDED / DENIED",
            Self::Reclaimed=>"RECLAIMED",
        }
    }
}
#[derive(Clone,Copy,Debug)]
pub struct Front {
    pub since_outbreak_years:f64,
    pub beta:f64,
}
impl Front {
    pub fn new(years:f64)->Self{
        Self{since_outbreak_years:years.clamp(0.0,TIMELINE_MAX_YR),
            beta:FRONT_BETA}
    }
    /// Causal geometric reach, ignoring defense. At 100 kyr <=10,000 ly.
    pub fn reach_ly(self)->f64 {
        self.since_outbreak_years.max(0.0)*self.beta
    }
    /// Angular difference in (-pi, pi].
    fn angle_delta(a:f64,b:f64)->f64{
        ((a-b+PI).rem_euclid(TAU))-PI
    }
    /// Delay in years from distributed defensive resistance.
    /// The slow sectors are precisely why the leading tendrils only
    /// *approach* the causal circle; they cannot overrun it.
    pub fn resistance_delay_yr(self,theta:f64)->f64{
        let age=(self.since_outbreak_years/100_000.0).clamp(0.0,1.0);
        let delayed=(age*age*(3.0-2.0*age));
        let seed=3_000.0+1_000.0*(3.0*theta+0.4).sin().abs();
        let pockets=[
            (0.15,14_000.0,0.29),
            (1.68,22_000.0,0.38),
            (3.00,16_000.0,0.27),
            (4.42,19_000.0,0.34),
            (5.55,12_000.0,0.24),
        ];
        let mut counter=0.0;
        for (center,power,width) in pockets {
            let delta=Self::angle_delta(theta,center)/width;
            counter+=power*(-0.5*delta*delta).exp();
        }
        let thin_lanes=(theta*9.0+0.47).sin().powi(12);
        let resistance=(seed+counter)*(1.0-0.84*thin_lanes);
        delayed*resistance
    }
    /// Most advanced *settled* sector front after millennia of
    /// defensive delays. Always inside the free 0.10c envelope.
    pub fn sector_reach_ly(self,theta:f64)->f64 {
        (self.reach_ly()-self.beta*self.resistance_delay_yr(theta))
            .max(0.0)
            .min(self.reach_ly())
    }
    /// Blight occupation is intermittent behind the leading front.
    /// Returned classification is a strategic indicator, not a
    /// probability of each individual star's capture.
    pub fn status(self,p:V2)->SectorStatus {
        let delta=p-ORION_OUTBREAK;
        let r=delta.len();
        let causal=self.reach_ly();
        if r>causal||causal<1.0{return SectorStatus::Unaffected;}
        let theta=delta.y.atan2(delta.x);
        let achieved=self.sector_reach_ly(theta);
        if r>achieved {
            return if self.since_outbreak_years>=50_000.0 &&
                r<causal && r>achieved*0.83 {
                SectorStatus::Defended
            }else{SectorStatus::Unaffected};
        }
        // A few surviving post-human zones break continuity in the
        // occupied interior. These appear only in the late era.
        if self.since_outbreak_years>70_000.0 {
            let protected=[
                (V2::new(-6_900.0,-29_400.0),1_650.0),
                (V2::new(3_500.0,-20_200.0),1_350.0),
            ];
            for (center,radius) in protected {
                if (p-center).len()<radius {
                    return SectorStatus::Reclaimed;
                }
            }
        }
        let ragged=280.0+250.0*(11.0*theta).sin().abs();
        if r>(achieved-ragged).max(0.0) {
            SectorStatus::Contested
        }else{
            SectorStatus::Held
        }
    }
}
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Sector {
    pub galactic_angle_index:usize,
    pub galactic_ring_index:usize,
    pub local_angle_index:usize,
    pub local_ring_index:usize,
    pub from_core_ly:f64,
    pub from_outbreak_ly:f64,
}
pub fn sector_of(p:V2)->Sector {
    let galactic_r=p.len();
    let local=p-ORION_OUTBREAK;
    let to_clock_deg=|v:V2| {
        // Screen compass: 0 deg at Galactic north/up, clockwise
        // toward the right-hand side of the overhead atlas.
        v.x.atan2(v.y).to_degrees().rem_euclid(360.0)
    };
    Sector{
        galactic_angle_index:(to_clock_deg(p)/GALACTIC_SECTOR_DEG)
            .floor() as usize%12,
        galactic_ring_index:(galactic_r/GALACTIC_RING_LY)
            .floor() as usize,
        local_angle_index:(to_clock_deg(local)/LOCAL_SECTOR_DEG)
            .floor() as usize%24,
        local_ring_index:(local.len()/LOCAL_RING_LY)
            .floor() as usize,
        from_core_ly:galactic_r,
        from_outbreak_ly:local.len(),
    }
}
#[derive(Clone,Copy,Debug)]
pub struct Waypoint {
    pub name:&'static str,
    pub description:&'static str,
    pub pos:V2,
    pub source:&'static str,
}
/// Region placement is for navigable game staging only.
pub const WAYPOINTS:[Waypoint;9]=[
    Waypoint{name:"SOL / ORION SPUR",
        description:"Human continuum and the old outbreak corridor",
        pos:SOL,source:"Story XV–XVII"},
    Waypoint{name:"AMICUS",
        description:"Group-mind protector; last baseline-human evacuation",
        pos:V2::new(-6_000.0,-29_000.0),source:"Episodes VIII–IX"},
    Waypoint{name:"OLYMPUS",
        description:"Encircled Titan system in the Perseus Arm",
        pos:V2::new(-15_000.0,-34_000.0),source:"Episodes XI–XIII"},
    Waypoint{name:"THE GREAT LIBRARY",
        description:"Ancestor archives, Dyson industry and transcendence",
        pos:V2::new(-20_000.0,-35_000.0),source:"Episode XIV"},
    Waypoint{name:"OUTER ARM EXODUS",
        description:"Refugees and wild-black-hole Halo Drives",
        pos:V2::new(-29_000.0,-40_000.0),source:"Episodes XVI–XVII"},
    Waypoint{name:"BEL-AIR",
        description:"Ringworld of Alfonso and Tatyana",
        pos:V2::new(-38_500.0,-35_000.0),source:"Episode XVIII"},
    Waypoint{name:"NEBULAR WASTES",
        description:"Weena's Chamber of Guf / Formoi and Patricia Tries",
        pos:V2::new(-38_000.0,-30_000.0),source:"Episodes XIX–XX"},
    Waypoint{name:"PERSEUS FRONT",
        description:"Dim-sky warning, Sauronoids and changing fronts",
        pos:V2::new(-11_000.0,-34_500.0),source:"Episodes XI, XXI"},
    Waypoint{name:"CANIS MAJOR EXODUS",
        description:"Extragalactic journey: schematic projected exit",
        pos:V2::new(18_000.0,-49_000.0),source:"Episode X"},
];
#[derive(Clone,Copy,Debug)]
pub struct TravelLeg{
    pub start:V2,
    pub destination:V2,
    pub cruise_beta:f64,
    pub total_distance_ly:f64,
    pub distance_travelled_ly:f64,
    pub proper_elapsed_yr:f64,
}
impl TravelLeg {
    pub fn new(start:V2,destination:V2,beta:f64)->Option<Self>{
        if !(0.01..1.0).contains(&beta){return None;}
        let distance=(destination-start).len();
        if !distance.is_finite()||distance<1.0{return None;}
        Some(Self{start,destination,cruise_beta:beta,
            total_distance_ly:distance,distance_travelled_ly:0.0,
            proper_elapsed_yr:0.0})
    }
    pub fn coordinate_duration_yr(self)->f64{
        self.total_distance_ly/self.cruise_beta
    }
    pub fn proper_duration_yr(self)->f64{
        self.coordinate_duration_yr()*
            (1.0-self.cruise_beta*self.cruise_beta).sqrt()
    }
    pub fn remaining_yr(self)->f64{
        (self.total_distance_ly-self.distance_travelled_ly)
            /self.cruise_beta
    }
    pub fn position(self)->V2{
        let u=(self.distance_travelled_ly/self.total_distance_ly)
            .clamp(0.0,1.0);
        self.start+(self.destination-self.start)*u
    }
    /// Returns excess coordinate time after arrival.
    pub fn advance(&mut self,years:f64)->f64{
        assert!(years>=0.0&&years.is_finite());
        let consumed=years.min(self.remaining_yr());
        self.distance_travelled_ly=(self.distance_travelled_ly+
            consumed*self.cruise_beta).min(self.total_distance_ly);
        self.proper_elapsed_yr+=consumed*
            (1.0-self.cruise_beta*self.cruise_beta).sqrt();
        years-consumed
    }
    pub fn arrived(self)->bool {
        self.distance_travelled_ly>=self.total_distance_ly-1e-8
    }
}
#[derive(Clone,Debug)]
pub struct GalacticNavModel {
    pub year:f64,
    pub rate_yr_per_second:f64,
    pub paused:bool,
    pub position:V2,
    pub destination:Option<V2>,
    pub travel:Option<TravelLeg>,
    pub proper_time_total_yr:f64,
    pub last_arrival:Option<V2>,
}
impl Default for GalacticNavModel{
    fn default()->Self{
        Self{year:100_000.0,rate_yr_per_second:700.0,
            paused:true,position:SOL,destination:None,
            travel:None,proper_time_total_yr:0.0,last_arrival:None}
    }
}
impl GalacticNavModel{
    pub fn front(&self)->Front{Front::new(self.year)}
    pub fn choose(&mut self,p:V2){
        self.destination=Some(p);
    }
    pub fn depart(&mut self,beta:f64)->bool{
        let Some(goal)=self.destination else{return false;};
        let Some(leg)=TravelLeg::new(self.position,goal,beta)
            else{return false;};
        self.travel=Some(leg);
        self.last_arrival=None;
        self.paused=false;
        true
    }
    pub fn tick(&mut self,frame_seconds:f64){
        if self.paused{return;}
        let years=(frame_seconds.clamp(0.0,0.1)*
            self.rate_yr_per_second)
            .min(TIMELINE_MAX_YR-self.year);
        if years<=0.0{self.paused=true;return;}
        self.year+=years;
        if let Some(mut leg)=self.travel {
            let proper_before=leg.proper_elapsed_yr;
            leg.advance(years);
            self.position=leg.position();
            self.proper_time_total_yr+=leg.proper_elapsed_yr-proper_before;
            if leg.arrived(){
                self.last_arrival=Some(leg.destination);
                self.travel=None;
                self.paused=true;
            }else{self.travel=Some(leg);}
        }
    }
    pub fn seek(&mut self,years:f64){
        self.year=years.clamp(0.0,TIMELINE_MAX_YR);
        self.paused=true;
        // Explicit epoch jumps reset any active travel whose
        // coordinate-time history would otherwise be ambiguous.
        self.travel=None;
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn no_blight_beyond_10000_ly_at_100000_years(){
        let front=Front::new(100_000.0);
        assert!((front.reach_ly()-10_000.0).abs()<1e-9);
        for a in 0..720 {
            let theta=a as f64*TAU/720.0;
            let local=front.sector_reach_ly(theta);
            assert!(local>=0.0&&local<=10_000.0+1e-9);
            let out=ORION_OUTBREAK+
                V2::from_angle(theta)*(10_000.0+0.01);
            assert_eq!(front.status(out),SectorStatus::Unaffected);
        }
    }
    #[test]
    fn directional_resistance_produces_inward_notches_and_tendrils(){
        let f=Front::new(100_000.0);
        let r=(0..1440).map(|i|
            f.sector_reach_ly(i as f64*TAU/1440.0))
            .collect::<Vec<_>>();
        let min=r.iter().copied().fold(f64::INFINITY,f64::min);
        let max=r.iter().copied().fold(0.0_f64,f64::max);
        assert!(max-min>500.0);
        assert!(max<=10_000.0);
        assert!(min>6_000.0);
    }
    #[test]
    fn chronology_and_relativistic_travel_have_separate_clocks(){
        let mut leg=TravelLeg::new(V2::default(),
            V2::new(5_000.0,0.0),0.76).unwrap();
        assert!((leg.coordinate_duration_yr()-6578.94737).abs()<0.001);
        let gamma=1.0/(1.0-0.76_f64.powi(2)).sqrt();
        assert!((leg.proper_duration_yr()*
            gamma-leg.coordinate_duration_yr()).abs()<1e-8);
        let extra=leg.advance(7_000.0);
        assert!(leg.arrived()&&extra>420.0);
        assert!((leg.position().x-5000.0).abs()<1e-8);
    }
    #[test]
    fn sector_bins_are_stable_in_both_reference_frames(){
        let s=sector_of(ORION_OUTBREAK+V2::new(4_000.0,0.0));
        assert_eq!(s.local_ring_index,2);
        assert_eq!(s.local_angle_index,6);
        assert!(s.galactic_angle_index<12);
    }
    #[test]
    fn chronological_seek_clears_in_flight_route(){
        let mut m=GalacticNavModel::default();
        m.choose(V2::new(5000.0,-27000.0));
        assert!(m.depart(0.76));
        m.tick(1.0);
        assert!(m.year>100_000.0);
        m.seek(55_000.0);
        assert!(m.travel.is_none()&&m.paused);
    }
}
