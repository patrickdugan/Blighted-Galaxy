//! BG-STR-05: versioned JSON campaign persistence.
//! Save the authoritative model state, not GPU textures or derived map pixels.
//! All untrusted save content is checked before it is admitted into gameplay.
use crate::overworld::{GalacticNavModel,TravelLeg,TIMELINE_MAX_YR};
use crate::strategy::{CargoConvoy,Faction,Signal,SignalKind,StrategicState,System};
use crate::V2;
use crate::guf_geometry::SlewMount;
use serde::{Deserialize,Serialize};
use std::{error::Error,fmt,fs,path::Path};

pub const SAVE_FORMAT:&str="BLIGHTED-GALAXY-CAMPAIGN";
pub const SAVE_VERSION:u32=2;
const MAX_PACKETS:usize=4096;
const MAX_EVENTS:usize=256;

#[derive(Debug)]
pub enum SaveError{
    Io(std::io::Error),
    Json(serde_json::Error),
    Invalid(String),
}
impl fmt::Display for SaveError {
    fn fmt(&self,f:&mut fmt::Formatter<'_>)->fmt::Result{
        match self {
            Self::Io(e)=>write!(f,"save I/O: {e}"),
            Self::Json(e)=>write!(f,"save JSON: {e}"),
            Self::Invalid(e)=>write!(f,"save invalid: {e}"),
        }
    }
}
impl Error for SaveError {}
impl From<std::io::Error> for SaveError{
    fn from(e:std::io::Error)->Self{Self::Io(e)}
}
impl From<serde_json::Error> for SaveError{
    fn from(e:serde_json::Error)->Self{Self::Json(e)}
}
fn bad(msg:&str)->SaveError{SaveError::Invalid(msg.into())}
fn finite(v:f64)->bool{v.is_finite()}
fn world(p:V2)->bool{finite(p.x)&&finite(p.y)&&p.len()<1.0e8}
fn interval(v:f64,lo:f64,hi:f64)->bool{finite(v)&&v>=lo&&v<=hi}

#[derive(Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
struct SystemSave{
    owner:Faction,industry:f64,energy:f64,material:f64,
    kugelblitz:u32,array_node:bool,defense:f64,guf_ward:f64,
}
impl From<&System> for SystemSave{
    fn from(s:&System)->Self{
        Self{owner:s.owner,industry:s.industry,energy:s.energy,
            material:s.material,kugelblitz:s.kugelblitz,array_node:s.array_node,
            defense:s.defense,guf_ward:s.guf_ward}
    }
}
impl SystemSave {
    fn apply(self,site:&mut System){
        site.owner=self.owner;site.industry=self.industry;
        site.energy=self.energy;site.material=self.material;
        site.kugelblitz=self.kugelblitz;site.array_node=self.array_node;
        site.defense=self.defense;site.guf_ward=self.guf_ward;
    }
    fn valid(&self)->bool{
        interval(self.industry,0.0,1000.0)&&interval(self.energy,0.0,1000.0)
        &&interval(self.material,0.0,600.0)&&self.kugelblitz<=1_000_000
        &&interval(self.defense,-100_000.0,100_000.0)
        &&interval(self.guf_ward,0.0,1_000_000.0)
    }
}
#[derive(Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
struct StrategySave{
    year:f64,systems:Vec<SystemSave>,signals:Vec<Signal>,
    convoys:Vec<CargoConvoy>,next_convoy_id:u64,
    events:Vec<String>,array_charged:f64,pulses:u32,
    historical:bool,blight_autonomy:bool,guf_mounts:Vec<SlewMount>,
    territory_revision:u64,last_tactical_year:Vec<Option<f64>>,
}
#[derive(Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
struct CampaignFile{
    format:String,version:u32,site_fingerprint:u64,
    nav:GalacticNavModel,strategy:StrategySave,
}
/// A stable fingerprint of canonical waypoint identifiers and geometry.
/// Versioned save slots deliberately reject incompatible maps.
fn fingerprint(sites:&[System])->u64{
    let mut h:u64=0xcbf29ce484222325;
    for site in sites {
        for byte in site.name.as_bytes().iter().copied()
            .chain(site.pos.x.to_bits().to_le_bytes())
            .chain(site.pos.y.to_bits().to_le_bytes()) {
            h^=u64::from(byte);
            h=h.wrapping_mul(0x100000001b3);
        }
    }
    h
}
fn capture(nav:&GalacticNavModel,strategy:&StrategicState)->CampaignFile{
    CampaignFile{
        format:SAVE_FORMAT.into(),version:SAVE_VERSION,
        site_fingerprint:fingerprint(&strategy.systems),nav:nav.clone(),
        strategy:StrategySave{
            year:strategy.year,
            systems:strategy.systems.iter().map(SystemSave::from).collect(),
            signals:strategy.signals.clone(),convoys:strategy.convoys.clone(),
            next_convoy_id:strategy.next_convoy_id,
            events:strategy.events.clone(),array_charged:strategy.array_charged,
            pulses:strategy.pulses,historical:strategy.historical,
            blight_autonomy:strategy.blight_autonomy,
            guf_mounts:strategy.guf_mounts.clone(),
            territory_revision:strategy.territory_revision,
            last_tactical_year:strategy.last_tactical_year.iter().map(|&t|
                if t==f64::NEG_INFINITY{None}else{Some(t)}).collect(),
        }
    }
}
fn valid_leg(t:TravelLeg)->bool{
    world(t.start)&&world(t.destination)
        &&interval(t.cruise_beta,0.001,0.999999)
        &&interval(t.total_distance_ly,1.0,200_000.0)
        &&interval(t.distance_travelled_ly,0.0,t.total_distance_ly)
        &&interval(t.proper_elapsed_yr,0.0,TIMELINE_MAX_YR)
        &&(t.total_distance_ly-(t.destination-t.start).len()).abs()<1e-6
}
fn validate_nav(nav:&GalacticNavModel)->Result<(),SaveError>{
    if !interval(nav.year,0.0,TIMELINE_MAX_YR)
        ||!interval(nav.rate_yr_per_second,0.01,100_000.0)
        ||!interval(nav.proper_time_total_yr,0.0,1e9)
        ||!world(nav.position)
        ||nav.destination.is_some_and(|p|!world(p))
        ||nav.last_arrival.is_some_and(|p|!world(p)){
        return Err(bad("invalid navigation clock, position or rate"));
    }
    if let Some(t)=nav.travel {
        if !valid_leg(t)|| (nav.position-t.position()).len()>1e-4{
            return Err(bad("invalid in-flight relativistic leg"));
        }
    }
    Ok(())
}
fn validate_packet(p:&Signal,systems:&[System],year:f64)->bool{
    if p.source>=systems.len()||p.target>=systems.len()
        ||!interval(p.launch_year,0.0,TIMELINE_MAX_YR)
        ||!interval(p.arrival_year,p.launch_year,TIMELINE_MAX_YR)
        ||p.arrival_year<year||!interval(p.power,0.0,1e7){return false;}
    let length=(systems[p.source].pos-systems[p.target].pos).len();
    let duration=p.arrival_year-p.launch_year;
    let limit=match p.kind{
        SignalKind::ArrayPulse|SignalKind::WeenaClaim=>length,
        SignalKind::BlightSeed=>length/0.10,
        SignalKind::HaloStrike=>length/0.76,
    };
    let aim_valid=match p.kind{
        SignalKind::ArrayPulse=>p.beam_azimuth_rad.is_some_and(|az|
            finite(az)&&az.abs()<100.0),
        _=>p.beam_azimuth_rad.is_none(),
    };
    aim_valid&&duration+1e-6>=limit
}
fn validate_convoy(c:&CargoConvoy,systems:&[System],year:f64)->bool{
    if c.source>=systems.len()||c.target>=systems.len()
        ||c.source==c.target||c.id==0
        ||!interval(c.launch_year,0.0,year)
        ||!interval(c.arrival_year,year,TIMELINE_MAX_YR)
        ||!interval(c.beta,0.01,0.999999)
        ||!interval(c.material,0.0,600.0)||!interval(c.energy,0.0,1000.0){
        return false;
    }
    let duration=(systems[c.source].pos-systems[c.target].pos).len()/c.beta;
    (c.arrival_year-c.launch_year-duration).abs()<=1e-6
}
pub fn encode(nav:&GalacticNavModel,
    strategy:&StrategicState)->Result<String,SaveError>{
    validate_nav(nav)?;
    if (nav.year-strategy.year).abs()>1e-7 {
        return Err(bad("navigation and strategic clocks differ"));
    }
    let raw=serde_json::to_string_pretty(&capture(nav,strategy))?;
    // Validate the emitted save too, including its finite fields.
    let _=decode(&raw)?;
    Ok(raw)
}
pub fn decode(raw:&str)->Result<(GalacticNavModel,StrategicState),SaveError>{
    if raw.len()>8*1024*1024{return Err(bad("oversized campaign file"));}
    let file:CampaignFile=serde_json::from_str(raw)?;
    if file.format!=SAVE_FORMAT||file.version!=SAVE_VERSION{
        return Err(bad("unsupported save format/version"));
    }
    validate_nav(&file.nav)?;
    let data=file.strategy;
    if !interval(data.year,0.0,TIMELINE_MAX_YR)
        ||(file.nav.year-data.year).abs()>1e-7{
        return Err(bad("mismatched galactic clocks"));
    }
    let mut state=StrategicState::new(data.year);
    if file.site_fingerprint!=fingerprint(&state.systems){
        return Err(bad("incompatible canonical star map"));
    }
    if data.systems.len()!=state.systems.len()
        ||data.last_tactical_year.len()!=state.systems.len()
        ||data.guf_mounts.len()!=state.systems.len(){
        return Err(bad("incorrect number of industrial systems"));
    }
    if data.signals.len()>MAX_PACKETS||data.convoys.len()>MAX_PACKETS
        ||data.events.len()>MAX_EVENTS||data.events.iter().any(|e|e.len()>512){
        return Err(bad("oversized packets, convoy queue or event log"));
    }
    if !interval(data.array_charged,0.0,100.0) ||data.next_convoy_id==0 {
        return Err(bad("invalid Guf charge or convoy sequence"));
    }
    for (saved,site) in data.systems.into_iter().zip(state.systems.iter_mut()){
        if !saved.valid(){return Err(bad("invalid industrial resource balance"));}
        saved.apply(site);
    }
    if !data.signals.iter().all(|p|validate_packet(p,&state.systems,state.year)){
        return Err(bad("invalid attack route or light-cone time"));
    }
    if !data.convoys.iter().all(|c|validate_convoy(c,&state.systems,state.year)){
        return Err(bad("invalid convoy route or sublight time"));
    }
    let newest=data.convoys.iter().map(|c|c.id).max().unwrap_or(0);
    if data.next_convoy_id<=newest{
        return Err(bad("convoy sequence may duplicate existing cargo"));
    }
    if !data.guf_mounts.iter().all(|m|
        m.valid(state.year,state.systems.len())){
        return Err(bad("Guf mount targeting or slew clock invalid"));
    }
    if !data.last_tactical_year.iter().all(|t|
        t.is_none_or(|v|interval(v,0.0,state.year))){
        return Err(bad("tactical cooldown timestamp is invalid"));
    }
    state.signals=data.signals;
    state.convoys=data.convoys;
    state.next_convoy_id=data.next_convoy_id;
    state.events=data.events;
    state.array_charged=data.array_charged;
    state.pulses=data.pulses;
    state.historical=data.historical;
    state.blight_autonomy=data.blight_autonomy;
    state.territory_revision=data.territory_revision;
    state.guf_mounts=data.guf_mounts;
    state.last_tactical_year=data.last_tactical_year.into_iter()
        .map(|x|x.unwrap_or(f64::NEG_INFINITY)).collect();
    Ok((file.nav,state))
}
pub fn save_file(path:&Path,nav:&GalacticNavModel,
    strategy:&StrategicState)->Result<(),SaveError>{
    let text=encode(nav,strategy)?;
    let tmp=path.with_extension("json.tmp");
    let backup=path.with_extension("json.bak");
    fs::write(&tmp,text.as_bytes())?;
    // Keep the previous valid slot as a recovery point on Windows too.
    if path.exists(){
        if backup.exists(){fs::remove_file(&backup)?;}
        fs::rename(path,&backup)?;
    }
    if let Err(err)=fs::rename(&tmp,path){
        if backup.exists(){let _=fs::rename(&backup,path);}
        return Err(SaveError::Io(err));
    }
    Ok(())
}
pub fn load_file(path:&Path)->Result<(GalacticNavModel,StrategicState),SaveError>{
    let data=fs::read_to_string(path)?;
    decode(&data)
}

#[cfg(test)]
mod tests{
    use super::*;
    use crate::strategy::{TacticalKind,TacticalOutcome};
    #[test]
    fn round_trip_preserves_resources_in_flight_convoys_guf_and_cooldowns(){
        let mut nav=GalacticNavModel::default();
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.systems[3].material=350.0;s.systems[3].energy=500.0;
        s.send_convoy(3,5).unwrap();
        for i in [2,3,6]{
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        s.blight_autonomy=false;
        let ready=s.aim_guf(5).unwrap();
        s.advance_to(ready);
        s.array_charged=100.0;s.fire_guf(5).unwrap();
        s.systems[5].energy=500.0;
        s.resolve_tactical(5,TacticalKind::GufDefense,TacticalOutcome::Victory).unwrap();
        nav.seek(s.year);
        nav.choose(crate::overworld::ORION_OUTBREAK);
        assert!(nav.depart(0.76));
        nav.tick(0.05);
        s.advance_to(nav.year);
        let original=encode(&nav,&s).unwrap();
        let (n,t)=decode(&original).unwrap();
        let restored=encode(&n,&t).unwrap();
        assert_eq!(original,restored);
        assert_eq!(t.convoys.len(),1);
        assert_eq!(t.incoming_guf(5).unwrap().1,3);
        assert_eq!(t.last_tactical_year[5],s.last_tactical_year[5]);
        assert!(n.travel.is_some());
    }
    #[test]
    fn save_rejects_forged_resource_and_bad_version(){
        let nav=GalacticNavModel::default();
        let s=StrategicState::default();
        let data=encode(&nav,&s).unwrap();
        assert!(decode(&data.replace("\"version\": 1","\"version\": 9")).is_err());
        assert!(decode(&data.replace("\"energy\": 100.0","\"energy\": -100.0")).is_err());
    }
    #[test]
    fn saving_before_any_tactical_engagement_handles_negative_infinity(){
        let nav=GalacticNavModel::default();
        let s=StrategicState::default();
        let raw=encode(&nav,&s).unwrap();
        assert!(raw.contains("\"last_tactical_year\""));
        let (_,r)=decode(&raw).unwrap();
        assert_eq!(r.last_tactical_year[0],f64::NEG_INFINITY);
    }
    #[test]
    fn convoy_id_cannot_be_replayed_by_file_edit(){
        let nav=GalacticNavModel::default();let mut s=StrategicState::default();
        s.systems[3].material=400.0;s.systems[3].energy=400.0;
        s.send_convoy(3,5).unwrap();
        s.next_convoy_id=1;
        assert!(encode(&nav,&s).is_err());
    }
}
