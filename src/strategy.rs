//! BG-STR-01: deterministic thousand-year strategic scenario layer.
//! Abstract resource units are game balance parameters, not astrophysical yields.
//! All messages and attacks propagate no faster than light in galactic coordinates.
use crate::overworld::{Front,ORION_OUTBREAK,WAYPOINTS,SectorStatus};
use crate::V2;
use crate::guf_geometry::{self,SlewMount};
use crate::episodes::EpisodeTracker;
use crate::patricia_trie::{PatriciaTrieWorld,TrieEventKind};
use crate::morph_boxing::BoxingDebrief;
use serde::{Serialize,Deserialize};

#[derive(Clone,Copy,Debug,PartialEq,Eq,Serialize,Deserialize)]
pub enum Faction { Blight, Drummers, Weena, Independent }
impl Faction {
    pub fn label(self)->&'static str {
        match self { Self::Blight=>"JENKINS BLIGHT",Self::Drummers=>"DRUMMERS / REMNANT",
            Self::Weena=>"WEENA / FORMOI",Self::Independent=>"INDEPENDENT" }
    }
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Project { Industry, Kugelblitz, GufNode }
#[derive(Clone,Copy,Debug,PartialEq,Eq,Serialize,Deserialize)]
pub enum TacticalKind { Intercept, Tidal, GufDefense, ConvoyEscort, HaloStrike }
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum TacticalOutcome { Victory, Defeat, Withdrawn }
#[derive(Clone,Debug)]
pub struct System {
    pub name:&'static str, pub pos:V2, pub owner:Faction,
    pub industry:f64, pub energy:f64, pub material:f64,
    pub kugelblitz:u32, pub array_node:bool, pub defense:f64,
    /// Site-side shield capacity against incoming Guf beams (game units).
    pub guf_ward:f64,
}
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Signal {
    pub kind:SignalKind, pub source:usize, pub target:usize,
    pub launch_year:f64, pub arrival_year:f64, pub power:f64,
    pub beam_azimuth_rad:Option<f64>,
}
#[derive(Clone,Copy,Debug,PartialEq,Eq,Serialize,Deserialize)]
pub enum SignalKind { BlightSeed, ArrayPulse, WeenaClaim, HaloStrike }
/// Material and stored power in a finite-speed sublight transport.
#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct CargoConvoy {
    pub id:u64,pub source:usize,pub target:usize,pub owner:Faction,
    pub material:f64,pub energy:f64,pub launch_year:f64,
    pub arrival_year:f64,pub beta:f64,
    #[serde(default)]
    pub escorted:bool,
}
impl CargoConvoy {
    pub fn position(&self,year:f64,systems:&[System])->V2{
        let p=systems[self.source].pos;
        let q=systems[self.target].pos;
        let fraction=((year-self.launch_year)/
            (self.arrival_year-self.launch_year).max(1e-9)).clamp(0.0,1.0);
        p+(q-p)*fraction
    }
}
#[derive(Clone,Debug)]
pub struct StrategicState {
    pub year:f64, pub systems:Vec<System>,
    pub signals:Vec<Signal>,pub events:Vec<String>,
    pub convoys:Vec<CargoConvoy>,pub next_convoy_id:u64,
    pub array_charged:f64,pub pulses:u32,
    pub guf_mounts:Vec<SlewMount>,
    pub historical:bool,
    pub blight_autonomy:bool,
    pub weena_autonomy:bool,
    pub episodes:EpisodeTracker,
    pub patricia_trie:PatriciaTrieWorld,
    /// Incremented whenever a world changes its territorial footprint.
    pub territory_revision:u64,
    /// One settled result per industrial system per 250 galactic years.
    pub last_tactical_year:Vec<f64>,
}
impl Default for StrategicState { fn default()->Self {Self::new(100_000.0)} }
impl StrategicState {
    pub fn new(year:f64)->Self {
        // Waypoint positions and project ownership are schematic campaign setup.
        let owners=[Faction::Independent,Faction::Drummers,Faction::Independent,
            Faction::Drummers,Faction::Drummers,Faction::Drummers,Faction::Weena,
            Faction::Independent,Faction::Independent];
        let mut systems:Vec<System>=WAYPOINTS.iter().enumerate().map(|(i,w)|System {
            name:w.name,pos:w.pos,owner:owners[i],
            industry:if i==6 {3.0} else if matches!(i,3|4|5){2.0}else{0.5},
            energy:100.0,material:150.0,kugelblitz:0,guf_ward:0.0,
            array_node:false, defense:if i==5{4.0}else{2.0},
        }).collect();
        systems.push(System{name:"ORION SEED",pos:ORION_OUTBREAK,
            owner:Faction::Blight,industry:4.0,energy:250.0,material:250.0,
            kugelblitz:0,array_node:false,defense:25.0,guf_ward:0.0});
        let count=systems.len();
        let patricia_trie=PatriciaTrieWorld::new_for_systems(&systems,year);
        Self{year,systems,signals:vec![],convoys:vec![],next_convoy_id:1,
            events:vec![
            "STRATEGIC SCENARIO: industrial networks online".into()],
            array_charged:0.0,pulses:0,
            guf_mounts:vec![SlewMount::default();count],historical:false,
            blight_autonomy:true,weena_autonomy:false,
            episodes:EpisodeTracker::default(),patricia_trie,
            territory_revision:0,
            last_tactical_year:vec![f64::NEG_INFINITY;count]}
    }
    pub fn system_positions(&self)->Vec<V2>{
        self.systems.iter().map(|s|s.pos).collect()
    }
    pub fn evaluate_episodes(&mut self){
        let dispatches=self.episodes.evaluate(
            &self.systems,self.pulses,self.year,&self.patricia_trie);
        self.events.extend(dispatches);
    }
    pub fn broadcast_memetic_record(&mut self,source:usize,
        label:impl Into<String>)->Result<f64,&'static str>{
        let positions=self.system_positions();
        let eta=self.patricia_trie.broadcast_memetic_record(
            source,&positions,self.year,label)?;
        self.events.push(format!(
            "MEMETIC BROADCAST from {} / wave horizon T+{eta:.0}",
            self.systems[source].name));
        self.evaluate_episodes();
        Ok(eta)
    }
    pub fn record_boxing_debrief(&mut self,debrief:BoxingDebrief)->String{
        let site=debrief.site.filter(|&i|i<self.systems.len());
        let positions=self.system_positions();
        let msg=if let Some(idx)=site {
            let label=format!("CONTACT DUEL @ {} (pulses:{})",
                self.systems[idx].name,debrief.memetic_pulses_landed);
            if debrief.memetic_pulses_landed>0 || debrief.outcome==TacticalOutcome::Victory {
                let _=self.patricia_trie.record_event(
                    idx,&positions,self.year,TrieEventKind::ContactDuel,label);
            }
            match debrief.outcome {
                TacticalOutcome::Victory=>{
                    self.episodes.boxing_duels_won+=1;
                    self.episodes.tactical_victories+=1;
                    if self.tactical_ready(idx){
                        let s=&mut self.systems[idx];
                        if matches!(s.owner,Faction::Blight|Faction::Weena){
                            s.defense=(s.defense-2.5).max(0.0);
                            if s.defense<=0.0 {
                                s.owner=Faction::Drummers;s.array_node=false;
                                s.defense=1.0;
                                self.guf_mounts[idx]=SlewMount::default();
                            }
                        }else{
                            s.defense=(s.defense+1.0).min(18.0);
                            s.energy=(s.energy+25.0).min(1000.0);
                        }
                        self.last_tactical_year[idx]=self.year;
                        self.territory_revision+=1;
                    }
                    format!("{}: CONTACT BOXING victory / Patricia Trie record logged",
                        self.systems[idx].name)
                }
                TacticalOutcome::Defeat=>{
                    if self.tactical_ready(idx){
                        self.systems[idx].defense=(self.systems[idx].defense-1.0).max(0.0);
                        self.last_tactical_year[idx]=self.year;
                        self.territory_revision+=1;
                    }
                    format!("{}: CONTACT BOXING defeat / chassis withdrawn",
                        self.systems[idx].name)
                }
                TacticalOutcome::Withdrawn=>{
                    format!("{}: CONTACT BOXING disengaged{}",
                        self.systems[idx].name,
                        if debrief.memetic_pulses_landed>0{" / memetic contact logged"}else{""})
                }
            }
        }else{
            match debrief.outcome {
                TacticalOutcome::Victory=>"SANDBOX CONTACT BOXING victory".into(),
                TacticalOutcome::Defeat=>"SANDBOX CONTACT BOXING defeat".into(),
                TacticalOutcome::Withdrawn=>"SANDBOX CONTACT BOXING disengaged".into(),
            }
        };
        self.events.push(msg.clone());
        self.evaluate_episodes();
        msg
    }
    pub fn plan_weena_expansion(&mut self)->Option<String>{
        let action=crate::episodes::plan_weena_expansion(self);
        if let Some(ref msg)=action {
            self.events.push(msg.clone());
        }
        action
    }
    pub fn convoy_threatened(&self,convoy:&CargoConvoy)->bool{
        if convoy.owner!=Faction::Drummers {return false;}
        let a=self.systems[convoy.source].pos;
        let b=self.systems[convoy.target].pos;
        let cur=convoy.position(self.year,&self.systems);
        let ab=b-a;
        let ab_len_sq=ab.len_sq().max(1e-9);
        self.systems.iter().any(|star|{
            if !matches!(star.owner,Faction::Blight|Faction::Weena){return false;}
            if (star.pos-cur).len()<=10_000.0 {return true;}
            let t=((star.pos-a).dot(ab)/ab_len_sq).clamp(0.0,1.0);
            let closest=a+ab*t;
            (star.pos-closest).len()<=10_000.0
        })
    }
    /// Finds an active Drummer convoy associated with `site` (origin, destination,
    /// or within 2,500 ly of `site`), preferring unescorted threatened convoys.
    pub fn active_convoy_at_site(&self,site:usize)->Option<usize>{
        let site_pos=self.systems.get(site)?.pos;
        let mut matches:Vec<(u8,f64,usize)>=self.convoys.iter().enumerate()
            .filter(|(_,c)|c.owner==Faction::Drummers && (
                c.source==site || c.target==site ||
                (c.position(self.year,&self.systems)-site_pos).len()<=2_500.0
            ))
            .map(|(idx,c)|{
                let prio=match (!c.escorted,self.convoy_threatened(c)){
                    (true,true)=>0,
                    (true,false)=>1,
                    (false,true)=>2,
                    (false,false)=>3,
                };
                let dist=(c.position(self.year,&self.systems)-site_pos).len();
                (prio,dist,idx)
            }).collect();
        matches.sort_by(|a,b|a.0.cmp(&b.0).then_with(||a.1.total_cmp(&b.1)));
        matches.first().map(|x|x.2)
    }
    /// Finds an escortable Drummer convoy within intercept range of `ship_pos`
    /// or tied to `selected_site` when the ship is stationed there.
    pub fn escortable_convoy_for(&self,ship_pos:V2,selected_site:Option<usize>)->Option<(usize,usize)>{
        if let Some(site)=selected_site {
            if let Some(s)=self.systems.get(site) {
                if (ship_pos-s.pos).len()<=500.0 {
                    if let Some(cidx)=self.active_convoy_at_site(site){
                        return Some((site,cidx));
                    }
                }
            }
        }
        // Also allow direct mid-route convoy escort if the ship is within 2,500 ly of the convoy.
        let mut best:Option<(f64,usize,usize)>=None;
        for (cidx,c) in self.convoys.iter().enumerate(){
            if c.owner!=Faction::Drummers{continue;}
            let pos=c.position(self.year,&self.systems);
            let d=(pos-ship_pos).len();
            if d<=2_500.0 && best.is_none_or(|b|d<b.0){
                best=Some((d,c.target,cidx));
            }
        }
        best.map(|(_,site,cidx)|(site,cidx))
    }
    /// Site-scale occupation carved into the procedural historical atlas.
    /// Every resulting Blight pixel remains inside Orion's causal circle.
    /// This represents local control bubbles, not a new galactic wave speed.
    pub fn territory_status(&self,p:V2,front:Front)->SectorStatus{
        if (p-ORION_OUTBREAK).len()>front.reach_ly(){
            return SectorStatus::Unaffected;
        }
        if let Some(site)=self.systems.iter()
            .filter(|site| {
                let radius=(500.0+site.industry.max(0.0)*70.0
                    +site.defense.max(0.0)*50.0).clamp(300.0,1550.0);
                (p-site.pos).len()<=radius
            })
            .min_by(|a,b|(a.pos-p).len().total_cmp(&(b.pos-p).len())){
            return match site.owner {
                Faction::Blight=>SectorStatus::Held,
                Faction::Drummers=>SectorStatus::Reclaimed,
                Faction::Weena=>SectorStatus::Defended,
                Faction::Independent=>front.status(p),
            };
        }
        front.status(p)
    }
    pub fn tactical_ready(&self,site:usize)->bool{
        site<self.systems.len() && self.year-self.last_tactical_year[site]>=250.0
    }
    pub fn tactical_rearm_remaining(&self,site:usize)->Option<f64>{
        self.systems.get(site)?;
        Some((250.0-(self.year-self.last_tactical_year[site])).max(0.0))
    }
    pub fn resolve_tactical(&mut self,site:usize,kind:TacticalKind,
        outcome:TacticalOutcome)->Result<String,&'static str>{
        if site>=self.systems.len(){return Err("UNKNOWN TACTICAL SYSTEM");}
        if outcome!=TacticalOutcome::Withdrawn && !self.tactical_ready(site){
            return Err("SYSTEM REARMING AFTER PREVIOUS BATTLE");
        }
        let guf_threat=self.signals.iter().any(|packet|
            packet.kind==SignalKind::ArrayPulse&&packet.target==site&&
            packet.arrival_year>self.year);
        if kind==TacticalKind::GufDefense && !guf_threat {
            return Err("NO ACTIVE GUF THREAT");
        }
        let convoy_index=if kind==TacticalKind::ConvoyEscort {
            Some(self.active_convoy_at_site(site).ok_or("NO ACTIVE DRUMMER CONVOY TO ESCORT")?)
        }else{None};
        if kind==TacticalKind::HaloStrike &&
            !matches!(self.systems[site].owner,Faction::Blight|Faction::Weena) {
            return Err("HALO STRIKE REQUIRES HOSTILE TARGET");
        }
        let halo_supplier=if kind==TacticalKind::HaloStrike && outcome==TacticalOutcome::Victory {
            let target_pos=self.systems[site].pos;
            Some(self.systems.iter().enumerate()
                .filter(|(_,sys)|sys.owner==Faction::Drummers&&sys.kugelblitz>0)
                .min_by(|a,b|(a.1.pos-target_pos).len().total_cmp(&(b.1.pos-target_pos).len()))
                .map(|(i,_)|i)
                .ok_or("NO STORED DRUMMER KUGELBLITZ")?)
        }else{None};

        let name=self.systems[site].name;
        let descriptor=match kind {
            TacticalKind::Intercept=>"RKM INTERCEPT",
            TacticalKind::Tidal=>"WHITE-DWARF INTERDICTION",
            TacticalKind::GufDefense=>"GUF RELAY INTERCEPTION",
            TacticalKind::ConvoyEscort=>"CONVOY ESCORT",
            TacticalKind::HaloStrike=>"HALO STRIKE",
        };
        let message=match outcome {
            TacticalOutcome::Withdrawn=>{
                if kind==TacticalKind::HaloStrike {
                    format!("{name}: {descriptor} disengaged / kugelblitz retained")
                }else{
                    format!("{name}: {descriptor} disengaged / no strategic change")
                }
            },
            TacticalOutcome::Victory=>{
                if kind==TacticalKind::ConvoyEscort{
                    let cidx=convoy_index.unwrap();
                    self.convoys[cidx].escorted=true;
                    let cid=self.convoys[cidx].id;
                    let s=&mut self.systems[site];
                    s.energy=(s.energy+20.0).min(1000.0);
                    s.material=(s.material+15.0).min(600.0);
                    s.defense=(s.defense+0.5).min(18.0);
                    self.episodes.convoy_escorts_won+=1;
                    self.episodes.tactical_victories+=1;
                    format!("{name}: CONVOY #{cid} escorted / bandit raiders broken (+escort ward)")
                }else if kind==TacticalKind::HaloStrike{
                    let supplier_idx=halo_supplier.unwrap();
                    self.systems[supplier_idx].kugelblitz-=1;
                    let s=&mut self.systems[site];
                    s.defense-=4.5;
                    s.industry=(s.industry-0.75).max(0.25);
                    s.material=(s.material-45.0).max(0.0);
                    self.episodes.halo_strikes_won+=1;
                    self.episodes.halo_strikes_launched+=1;
                    self.episodes.tactical_victories+=1;
                    if s.defense<=0.0{
                        s.owner=Faction::Drummers;s.array_node=false;
                        self.guf_mounts[site]=SlewMount::default();
                        s.defense=1.0;
                        format!("{name}: HALO STRIKE victory / hostile stronghold shattered")
                    }else{
                        format!("{name}: HALO STRIKE victory / enemy defense {:.1}",s.defense)
                    }
                }else if kind==TacticalKind::GufDefense{
                    let s=&mut self.systems[site];
                    if s.owner!=Faction::Drummers{
                        return Err("DEFENSIVE INTERCEPTION REQUIRES DRUMMER CONTROL");
                    }
                    if s.material<25.0||s.energy<40.0{
                        return Err("INSUFFICIENT MATERIAL / ENERGY FOR GUF SCREEN");
                    }
                    s.material-=25.0;s.energy-=40.0;
                    let screen=(4.0+s.industry*1.5+s.defense.max(0.0)*0.4)
                        .clamp(4.0,18.0);
                    s.guf_ward+=screen;
                    self.episodes.guf_defenses_won+=1;
                    self.episodes.tactical_victories+=1;
                    format!("{name}: relay neutralized / GUF interception screen +{screen:.1}")
                }else if matches!(self.systems[site].owner,Faction::Blight|Faction::Weena){
                    let s=&mut self.systems[site];
                    let damage=if kind==TacticalKind::Intercept{3.0}else{2.0};
                    s.defense-=damage;
                    s.material=(s.material-25.0).max(0.0);
                    if kind==TacticalKind::Tidal {self.episodes.tidal_victories+=1;}
                    self.episodes.tactical_victories+=1;
                    if s.defense<=0.0{
                        s.owner=Faction::Drummers;s.array_node=false;
                        self.guf_mounts[site]=SlewMount::default();
                        s.industry=(s.industry*0.6).max(0.25);
                        s.defense=1.0;
                        format!("{name}: {descriptor} victory / hostile control broken")
                    }else{
                        format!("{name}: {descriptor} victory / enemy defense {:.1}",s.defense)
                    }
                }else{
                    let s=&mut self.systems[site];
                    s.defense=(s.defense+1.25).min(18.0);
                    s.energy=(s.energy+20.0).min(1000.0);
                    if kind==TacticalKind::Tidal {self.episodes.tidal_victories+=1;}
                    self.episodes.tactical_victories+=1;
                    format!("{name}: {descriptor} victory / defenses reinforced")
                }
            },
            TacticalOutcome::Defeat=>{
                if kind==TacticalKind::ConvoyEscort{
                    let cidx=convoy_index.unwrap();
                    let lost=self.convoys.remove(cidx);
                    let s=&mut self.systems[site];
                    s.defense=(s.defense-1.0).max(0.0);
                    format!("{name}: CONVOY #{} destroyed by bandit raiders",lost.id)
                }else{
                    let s=&mut self.systems[site];
                    s.defense=(s.defense-1.5).max(0.0);
                    s.material=(s.material-18.0).max(0.0);
                    if s.defense<=0.0 && s.owner==Faction::Drummers &&
                        (s.pos-ORION_OUTBREAK).len()<=Front::new(self.year).reach_ly(){
                        s.owner=Faction::Blight;s.array_node=false;
                        self.guf_mounts[site]=SlewMount::default();
                        s.defense=1.0;
                        format!("{name}: {descriptor} defeat / site lost to Blight")
                    }else{
                        format!("{name}: {descriptor} defeat / defenses damaged")
                    }
                }
            },
        };
        if outcome!=TacticalOutcome::Withdrawn{
            self.territory_revision+=1;
            self.last_tactical_year[site]=self.year;
        }
        self.events.push(message.clone());
        self.evaluate_episodes();
        Ok(message)
    }
    /// Orders a sublight transport between two Drummer-held industrial
    /// systems. Cargo is escrowed at launch, then delivered exactly once.
    pub fn send_convoy(&mut self,source:usize,target:usize)->Result<u64,&'static str>{
        if source==target{return Err("CONVOY NEEDS DISTINCT SYSTEMS");}
        let a=self.systems.get(source).ok_or("INVALID CONVOY SOURCE")?;
        let b=self.systems.get(target).ok_or("INVALID CONVOY DESTINATION")?;
        if a.owner!=Faction::Drummers||b.owner!=Faction::Drummers{
            return Err("BOTH CONVOY SYSTEMS MUST BE DRUMMER HELD");
        }
        if a.material<60.0||a.energy<65.0{
            return Err("CONVOY NEEDS 60 MATERIAL AND 65 ENERGY");
        }
        let beta=(0.35+0.04*a.industry).clamp(0.35,0.75);
        let distance=(a.pos-b.pos).len();
        if distance<1.0{return Err("INVALID CONVOY DISTANCE");}
        let arrival=self.year+distance/beta;
        if arrival>200_000.0{return Err("CONVOY ARRIVES AFTER CAMPAIGN HORIZON");}
        self.systems[source].material-=60.0;
        self.systems[source].energy-=65.0; // 40 cargo + 25 propulsion
        let id=self.next_convoy_id;
        self.next_convoy_id+=1;
        self.convoys.push(CargoConvoy{
            id,source,target,owner:Faction::Drummers,
            material:60.0,energy:40.0,launch_year:self.year,arrival_year:arrival,beta,
            escorted:false,
        });
        self.events.push(format!("DRUMMER CONVOY #{id}: {} to {} / ETA T+{arrival:.0}",
            self.systems[source].name,self.systems[target].name));
        Ok(id)
    }
    pub fn incoming_guf(&self,target:usize)->Option<(f64,usize,f64)>{
        self.systems.get(target)?;
        let pulses=self.signals.iter().filter(|packet|
            packet.kind==SignalKind::ArrayPulse &&
            packet.target==target&&packet.arrival_year>self.year);
        let mut earliest=f64::INFINITY;let mut count=0;let mut power=0.0;
        for packet in pulses {
            earliest=earliest.min(packet.arrival_year);count+=1;power+=packet.power;
        }
        if count>0{Some((earliest,count,power))}else{None}
    }
    pub fn nearest(&self,p:V2,max_distance_ly:f64)->Option<usize>{
        self.systems.iter().enumerate()
            .filter(|(_,s)|(s.pos-p).len()<=max_distance_ly)
            .min_by(|a,b|(a.1.pos-p).len().total_cmp(&(b.1.pos-p).len()))
            .map(|(i,_)|i)
    }
    pub fn project(&mut self,site:usize,project:Project)->Result<String,&'static str>{
        let s=self.systems.get_mut(site).ok_or("UNKNOWN SYSTEM")?;
        let name=s.name;
        match project {
            Project::Industry=>{
                if s.owner!=Faction::Drummers && s.owner!=Faction::Weena {
                    return Err("INDUSTRIAL CONTROL REQUIRED");
                }
                if s.material<45.0||s.energy<20.0{return Err("INSUFFICIENT MATERIAL / ENERGY");}
                s.material-=45.0;s.energy-=20.0;s.industry+=0.5;
                Ok(format!("{name}: industrial capacity expanded"))
            },
            Project::Kugelblitz=>{
                if s.owner!=Faction::Drummers{return Err("DRUMMER SYSTEM REQUIRED");}
                if s.industry<2.0{return Err("BUILD INDUSTRY FIRST");}
                if s.material<75.0||s.energy<115.0{return Err("KUGELBLITZ ENERGY / MATERIAL SHORTAGE");}
                s.material-=75.0;s.energy-=115.0;s.kugelblitz+=1;
                Ok(format!("{name}: kugelblitz contained / Halo Drive stock +1"))
            },
            Project::GufNode=>{
                if s.owner!=Faction::Weena{return Err("WEENA CONTROL REQUIRED");}
                if s.industry<2.0||s.material<65.0||s.energy<75.0 {
                    return Err("ARRAY NODE NEEDS INDUSTRY / STORES");
                }
                s.material-=65.0;s.energy-=75.0;s.array_node=true;
                Ok(format!("{name}: synchronized Guf relay prepared"))
            },
        }.map(|msg|{
            if project==Project::GufNode{
                self.guf_mounts[site]=SlewMount::default();
            }
            self.territory_revision+=1;self.events.push(msg.clone());
            self.evaluate_episodes();
            msg
        })
    }
    pub fn set_owner(&mut self,site:usize,owner:Faction){
        if let Some(s)=self.systems.get_mut(site){
            s.owner=owner;
            if owner!=Faction::Weena{
                s.array_node=false;
                self.guf_mounts[site]=SlewMount::default();
            }
            self.territory_revision+=1;
        }
    }
    pub fn schedule_blight_seed(&mut self,source:usize,target:usize)->Result<(),&'static str>{
        if self.signals.iter().any(|s|s.kind==SignalKind::BlightSeed&&s.target==target){
            return Err("BLIGHT PACKET ALREADY IN FLIGHT");
        }
        let (a,b)=(self.systems.get(source).ok_or("INVALID SOURCE")?,
                   self.systems.get(target).ok_or("INVALID TARGET")?);
        if a.owner!=Faction::Blight{return Err("SOURCE IS NOT BLIGHT");}
        if source==target||b.owner==Faction::Blight{return Err("TARGET ALREADY BLIGHT-HELD");}
        // Blight tendrils have effective expansion 0.10c; extra delay
        // encodes resistance. This also enforces the causal front geometry.
        let delay=(b.pos-a.pos).len()/0.10+1_500.0;
        let arrival=self.year+delay;
        if arrival>200_000.0{return Err("BLIGHT ARRIVAL AFTER CAMPAIGN HORIZON");}
        if (b.pos-ORION_OUTBREAK).len()>Front::new(arrival).reach_ly(){
            return Err("OUTSIDE ORION CAUSAL ENVELOPE");
        }
        self.signals.push(Signal{kind:SignalKind::BlightSeed,source,target,
            launch_year:self.year,arrival_year:arrival,power:1.0,
            beam_azimuth_rad:None});
        Ok(())
    }
    /// Use a Drummer-made micro black hole to propel a strike to a
    /// hostile star system. The signal is a 0.76c sublight mission.
    pub fn launch_halo_strike(&mut self,target:usize)->Result<f64,&'static str>{
        let dst=self.systems.get(target).ok_or("INVALID TARGET")?;
        if !matches!(dst.owner,Faction::Blight|Faction::Weena){
            return Err("HALO STRIKE REQUIRES HOSTILE TARGET");
        }
        let source=self.systems.iter().enumerate()
            .filter(|(_,s)|s.owner==Faction::Drummers&&s.kugelblitz>0&&s.energy>=20.0)
            .min_by(|a,b|(a.1.pos-dst.pos).len().total_cmp(&(b.1.pos-dst.pos).len()))
            .map(|(i,_)|i).ok_or("NO CHARGED DRUMMER KUGELBLITZ")?;
        let arrival=self.year+
            (self.systems[source].pos-dst.pos).len()/0.76;
        if arrival>200_000.0{return Err("STRIKE ARRIVES AFTER CAMPAIGN HORIZON");}
        self.systems[source].kugelblitz-=1;
        self.systems[source].energy-=20.0;
        self.signals.push(Signal{kind:SignalKind::HaloStrike,
            source,target,launch_year:self.year,
            arrival_year:arrival,power:5.0,beam_azimuth_rad:None});
        self.episodes.halo_strikes_launched+=1;
        self.events.push(format!("HALO STRIKE dispatched / target {} / T+{arrival:.0}",
            self.systems[target].name));
        self.evaluate_episodes();
        Ok(arrival)
    }
    /// Schedule political/industrial takeover, with a light-speed order
    /// from the nearest Weena-held system. Resistance adds 2 kyr.
    pub fn claim_for_weena(&mut self,target:usize)->Result<f64,&'static str>{
        let dst=self.systems.get(target).ok_or("INVALID TARGET")?;
        if dst.owner!=Faction::Independent && dst.owner!=Faction::Drummers {
            return Err("TARGET IS ALREADY OCCUPIED");
        }
        let source=self.systems.iter().enumerate()
            .filter(|(_,s)|s.owner==Faction::Weena&&s.energy>=80.0)
            .min_by(|a,b|(a.1.pos-dst.pos).len().total_cmp(&(b.1.pos-dst.pos).len()))
            .map(|(i,_)|i).ok_or("NO WEENA RELAY POWER")?;
        let arrival=self.year+(self.systems[source].pos-dst.pos).len()+2_000.0;
        if arrival>200_000.0{return Err("RECRUITMENT ARRIVES AFTER CAMPAIGN HORIZON");}
        self.systems[source].energy-=80.0;
        self.signals.push(Signal{kind:SignalKind::WeenaClaim,source,target,
            launch_year:self.year,arrival_year:arrival,power:1.0,
            beam_azimuth_rad:None});
        self.events.push(format!("WEENA influence relay sent to {}",self.systems[target].name));
        Ok(arrival)
    }
    /// (assigned, ready, latest finish) for the selected target.
    pub fn guf_alignment_status(&self,target:usize)->Option<(usize,usize,f64)>{
        let destination=self.systems.get(target)?.pos;
        let mut assigned=0;let mut aligned=0;let mut latest=self.year;
        for (idx,site) in self.systems.iter().enumerate(){
            if site.owner!=Faction::Weena||!site.array_node{continue;}
            let mount=&self.guf_mounts[idx];
            if mount.assigned_target==Some(target){
                assigned+=1;
                latest=latest.max(mount.slew_end_year);
                if mount.aligned(self.year,target,
                    guf_geometry::bearing(site.pos,destination)){aligned+=1;}
            }
        }
        Some((assigned,aligned,latest))
    }
    /// Point independent Guf relays toward a star. The order first
    /// propagates from Weena's controller at c, then each mount slews
    /// through the shortest angular path at a finite rotation rate.
    pub fn aim_guf(&mut self,target:usize)->Result<f64,&'static str>{
        let dest=self.systems.get(target).ok_or("INVALID GUF TARGET")?.pos;
        let emitters:Vec<usize>=self.systems.iter().enumerate()
            .filter(|(_,s)|s.owner==Faction::Weena&&s.array_node&&s.energy>=50.0)
            .map(|(i,_)|i).collect();
        if emitters.len()<3{return Err("THREE POWERED GUF NODES REQUIRED");}
        let controller=if self.systems.get(6).is_some_and(|s|s.owner==Faction::Weena) {
            6
        }else{emitters[0]};
        let origin=self.systems[controller].pos;
        let mut updates=Vec::new();
        let mut ready=self.year;
        for source in emitters {
            if source==target{return Err("TARGET OVERLAPS GUF EMITTER");}
            let desired=guf_geometry::bearing(self.systems[source].pos,dest);
            if self.guf_mounts[source].aligned(self.year,target,desired) {continue;}
            let receive=self.year+(self.systems[source].pos-origin).len();
            let next=self.guf_mounts[source]
                .schedule(self.year,receive,desired,target)
                .ok_or("GUF MOUNT BUSY WITH EXISTING SLEW")?;
            if next.slew_end_year>200_000.0{return Err("GUF ALIGNMENT BEYOND CAMPAIGN");}
            ready=ready.max(next.slew_end_year);
            updates.push((source,next));
        }
        for (source,mount) in updates{self.guf_mounts[source]=mount;}
        self.events.push(format!("GUF ALIGNMENT {} / READY T+{ready:.0}",
            self.systems[target].name));
        Ok(ready)
    }
    /// Requires three independently aligned industrial nodes, each with
    /// a narrow footprint and finite light-travel-time beam packet.
    /// The remaining beam is game-scale and may be intercepted on arrival.
    pub fn fire_guf(&mut self,target:usize)->Result<f64,&'static str>{
        let destination=self.systems.get(target).ok_or("INVALID TARGET")?;
        let target_pos=destination.pos;
        let emitters:Vec<usize>=self.systems.iter().enumerate()
            .filter(|(i,s)|s.owner==Faction::Weena&&s.array_node&&s.energy>=50.0&&
                self.guf_mounts[*i].aligned(self.year,target,
                    guf_geometry::bearing(s.pos,target_pos)))
            .map(|(i,_)|i).collect();
        if emitters.len()<3{return Err("THREE ALIGNED GUF NODES REQUIRED / SHIFT+P");}
        if self.array_charged<100.0{return Err("ARRAY CHARGING");}
        // Each emitter's light cone is independent. A distributed strike
        // reaches the destination at the latest packet arrival.
        let arrival=emitters.iter().map(|&i|
            self.year+(self.systems[i].pos-target_pos).len())
            .fold(self.year,f64::max);
        if arrival>200_000.0{return Err("GUF BEAMS ARRIVE AFTER CAMPAIGN HORIZON");}
        for &i in &emitters {self.systems[i].energy-=50.0;}
        self.array_charged=0.0;self.pulses+=1;
        // Delayed launches align all beam arrivals at the same target time.
        // Each separate flight is visible in the overworld and resolves once.
        for &source in &emitters {
            let light_time=(self.systems[source].pos-target_pos).len();
            self.signals.push(Signal{kind:SignalKind::ArrayPulse,
                source,target,launch_year:arrival-light_time,
                arrival_year:arrival,power:4.0,
                beam_azimuth_rad:Some(self.guf_mounts[source].bearing_at(self.year))});
        }
        let positions=self.system_positions();
        let _=self.patricia_trie.record_event(
            emitters[0],&positions,self.year,TrieEventKind::GufShock,
            format!("CHAMBER OF GUF DISCHARGE -> {}",self.systems[target].name));
        self.events.push(format!("CHAMBER OF GUF fired; target arrival T+{arrival:.0} yr"));
        self.evaluate_episodes();
        Ok(arrival)
    }
    /// Slow autonomous swarm campaign. Commands travel on existing 0.10c
    /// propagation rules and never exceed the Orion light-cone envelope.
    fn plan_blight_expansion(&mut self){
        if self.signals.iter().filter(|s|s.kind==SignalKind::BlightSeed).count()>=2{return;}
        let mut candidates=Vec::new();
        for (target,site) in self.systems.iter().enumerate(){
            if site.owner==Faction::Blight ||
                self.signals.iter().any(|p|p.kind==SignalKind::BlightSeed&&p.target==target){
                continue;
            }
            for (source,src) in self.systems.iter().enumerate(){
                if src.owner!=Faction::Blight||source==target{continue;}
                let eta=self.year+(src.pos-site.pos).len()/0.10+1_500.0;
                if eta>200_000.0 || (site.pos-ORION_OUTBREAK).len()>Front::new(eta).reach_ly(){
                    continue;
                }
                candidates.push(((src.pos-site.pos).len(),source,target));
            }
        }
        candidates.sort_by(|a,b|a.0.total_cmp(&b.0));
        if let Some((_,source,target))=candidates.first().copied(){
            if self.schedule_blight_seed(source,target).is_ok(){
                self.events.push(format!("BLIGHT autonomous growth order: {}",
                    self.systems[target].name));
            }
        }
    }
    pub fn advance_to(&mut self,year:f64){
        if !year.is_finite()||year<=self.year{return;}
        let end=year.min(200_000.0);
        while self.year<end {
            let previous_bin=(self.year/2_500.0).floor() as u64;
            let mut next=(self.year+250.0).min(end);
            // Resolve arrivals at their exact galactic coordinate time.
            // Economic integration uses shorter substeps at arrival events.
            for p in &self.signals{
                if p.arrival_year>self.year && p.arrival_year<next{
                    next=p.arrival_year;
                }
            }
            for c in &self.convoys{
                if c.arrival_year>self.year && c.arrival_year<next{
                    next=c.arrival_year;
                }
            }
            if let Some(trie_eta)=self.patricia_trie.next_arrival_after(self.year){
                if trie_eta<next {next=trie_eta;}
            }
            let next_plan=((previous_bin+1) as f64)*2500.0;
            if next_plan>self.year && next_plan<next{next=next_plan;}
            let dt=next-self.year;
            for s in &mut self.systems {
                let occupation=if s.owner==Faction::Blight{1.30}else{1.0};
                s.energy=(s.energy+s.industry*dt*0.035*occupation).min(1000.0);
                s.material=(s.material+s.industry*dt*0.013*occupation).min(600.0);
                if s.owner==Faction::Blight{s.industry=(s.industry+dt*0.00025).min(10.0);}
            }
            self.array_charged=(self.array_charged+dt*0.05*
                self.systems.iter().filter(|s|s.owner==Faction::Weena&&s.array_node)
                    .count() as f64).min(100.0);
            self.year=next;
            let mut pending=vec![];
            for packet in self.signals.drain(..){
                if packet.arrival_year>self.year {pending.push(packet);continue;}
                let intersects=packet.beam_azimuth_rad.is_some_and(|az|
                    guf_geometry::beam_intersects(
                        self.systems[packet.source].pos,
                        self.systems[packet.target].pos,az));
                if let Some(s)=self.systems.get_mut(packet.target){
                    match packet.kind {
                        SignalKind::WeenaClaim=>{
                            if s.owner!=Faction::Blight {
                                s.owner=Faction::Weena;s.array_node=false;
                                self.guf_mounts[packet.target]=SlewMount::default();
                                self.territory_revision+=1;
                                self.events.push(format!("WEENA recruited {}",s.name));
                            }
                        }
                        SignalKind::BlightSeed=>{
                            if (s.pos-ORION_OUTBREAK).len()<=Front::new(self.year).reach_ly(){
                                s.defense-=packet.power*2.0;
                                if s.defense<=0.0{
                                    s.owner=Faction::Blight;s.array_node=false;
                                    self.guf_mounts[packet.target]=SlewMount::default();
                                }
                                self.territory_revision+=1;
                                self.events.push(format!("BLIGHT wave reached {}",s.name));
                            }
                        }
                        SignalKind::HaloStrike=>{
                            if matches!(s.owner,Faction::Blight|Faction::Weena){
                                s.defense-=packet.power;
                                s.material=(s.material-35.0).max(0.0);
                                if s.defense<=0.0{
                                    s.owner=Faction::Drummers;
                                    s.array_node=false;
                                    self.guf_mounts[packet.target]=SlewMount::default();
                                    s.industry=(s.industry*0.6).max(0.25);
                                    s.defense=1.0;
                                }
                                self.territory_revision+=1;
                                self.events.push(format!("HALO STRIKE hit {}",s.name));
                            }
                        }
                        SignalKind::ArrayPulse=>{
                            if !intersects{
                                self.events.push(format!(
                                    "GUF beam missed {} / off-axis trajectory",s.name));
                                continue;
                            }
                            let absorbed=s.guf_ward.min(packet.power);
                            s.guf_ward=(s.guf_ward-absorbed).max(0.0);
                            let damage=packet.power-absorbed;
                            s.defense-=damage;
                            if damage>0.0{
                                self.territory_revision+=1;
                                if s.defense<=0.0 {
                                    s.industry=(s.industry*0.25).max(0.0);
                                    s.array_node=false;
                                    self.guf_mounts[packet.target]=SlewMount::default();
                                }
                            }
                            self.events.push(format!(
                                "GUF beam reached {} / absorbed {:.1} / damage {:.1}",
                                s.name,absorbed,damage));
                        }
                    }
                }
            }
            self.signals=pending;
            // Sample traffic at each 250-year substep; a star held by the
            // Blight may destroy a nearby convoy before it can deliver.
            let mut surviving=Vec::new();
            for mut convoy in std::mem::take(&mut self.convoys) {
                let p=convoy.position(self.year,&self.systems);
                let enemy_near=self.systems.iter().any(|star|
                    star.owner==Faction::Blight&&(star.pos-p).len()<850.0);
                if enemy_near && !convoy.escorted {
                    self.events.push(format!("DRUMMER CONVOY #{} interdicted",
                        convoy.id));
                }else{
                    if enemy_near && convoy.escorted {
                        convoy.escorted=false;
                        self.events.push(format!("DRUMMER CONVOY #{} escort repelled interdiction",
                            convoy.id));
                    }
                    if convoy.arrival_year<=self.year {
                        let dest=&mut self.systems[convoy.target];
                        if dest.owner==convoy.owner {
                            dest.material=(dest.material+convoy.material).min(600.0);
                            dest.energy=(dest.energy+convoy.energy).min(1000.0);
                            self.episodes.convoys_delivered+=1;
                            self.events.push(format!("CONVOY #{} delivered to {}",
                                convoy.id,dest.name));
                        }else{
                            self.events.push(format!("CONVOY #{} lost at hostile destination",
                                convoy.id));
                        }
                    }else{surviving.push(convoy);}
                }
            }
            self.convoys=surviving;
            let positions=self.system_positions();
            let trie_events=self.patricia_trie.advance_to(self.year,&positions);
            self.events.extend(trie_events);
            let crossed_plan_bin=(self.year/2_500.0).floor() as u64>previous_bin;
            if self.blight_autonomy && crossed_plan_bin {
                self.plan_blight_expansion();
            }
            if self.weena_autonomy && crossed_plan_bin {
                self.plan_weena_expansion();
            }
            self.evaluate_episodes();
        }
        if self.events.len()>14 {
            let extra=self.events.len()-14;
            self.events.drain(..extra);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn ready_array(s:&mut StrategicState,target:usize){
        s.blight_autonomy=false;
        let ready=s.aim_guf(target).unwrap();
        s.advance_to(ready);
        s.array_charged=100.0;
    }
    #[test] fn production_is_resource_limited(){
        let mut s=StrategicState::new(100_000.0);
        assert!(s.project(5,Project::Kugelblitz).is_err());
        s.advance_to(103_000.0);
        assert!(s.project(5,Project::Kugelblitz).is_ok());
        assert_eq!(s.systems[5].kugelblitz,1);
    }
    #[test] fn superluminal_array_fire_is_impossible(){
        let mut s=StrategicState::new(100_000.0);
        for i in [2,3,6] {s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=400.0;}
        ready_array(&mut s,5);
        let launch=s.year;
        let latest=[2,3,6].iter().map(|&i|(s.systems[i].pos-s.systems[5].pos).len())
            .fold(0.0_f64,f64::max);
        let arrival=s.fire_guf(5).unwrap();
        assert!((arrival-launch-latest).abs()<1e-6);
        assert_eq!(s.systems[5].owner,Faction::Drummers);
        s.advance_to(arrival);
        assert!(s.systems[5].defense<4.0);
    }
    #[test] fn independent_light_cone_limits_territory(){
        let mut s=StrategicState::new(100_000.0);
        s.set_owner(5,Faction::Blight); // outside the Orion causal radius
        assert_eq!(s.territory_status(s.systems[5].pos,Front::new(100_000.0)),
            SectorStatus::Unaffected);
        s.set_owner(1,Faction::Blight);
        assert_eq!(s.territory_status(s.systems[1].pos,Front::new(100_000.0)),
            SectorStatus::Held);
        s.resolve_tactical(1,TacticalKind::Intercept,TacticalOutcome::Victory).unwrap();
        assert_eq!(s.systems[1].owner,Faction::Drummers);
        assert_eq!(s.territory_status(s.systems[1].pos,Front::new(100_000.0)),
            SectorStatus::Reclaimed);
    }
    #[test] fn battles_cannot_be_farmed_by_withdrawing(){
        let mut s=StrategicState::default();
        let original=(s.systems[5].defense,s.systems[5].energy,s.territory_revision);
        s.resolve_tactical(5,TacticalKind::Tidal,TacticalOutcome::Withdrawn).unwrap();
        assert_eq!(original,(s.systems[5].defense,s.systems[5].energy,s.territory_revision));
    }
    #[test] fn distributed_pulse_respects_every_emitters_light_cone(){
        let mut s=StrategicState::new(100_000.0);
        for i in [2,3,6] {
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        ready_array(&mut s,5);
        let arrival=s.fire_guf(5).unwrap();
        assert_eq!(s.signals.len(),3);
        for packet in &s.signals {
            assert_eq!(packet.kind,SignalKind::ArrayPulse);
            assert!((packet.arrival_year-packet.launch_year-
                (s.systems[packet.source].pos-s.systems[5].pos).len()).abs()<1e-8);
            assert!((packet.arrival_year-arrival).abs()<1e-8);
        }
    }
    #[test] fn kugelblitz_attack_consumes_stock_and_arrives_at_sub_c(){
        let mut s=StrategicState::new(100_000.0);
        s.set_owner(1,Faction::Blight);
        s.systems[3].kugelblitz=1;
        let distance=(s.systems[1].pos-s.systems[3].pos).len();
        let arrival=s.launch_halo_strike(1).unwrap();
        assert_eq!(s.systems[3].kugelblitz,0);
        assert!((arrival-s.year-distance/0.76).abs()<1e-8);
        assert_eq!(s.systems[1].owner,Faction::Blight);
        s.advance_to(arrival-1.0);
        assert_eq!(s.systems[1].owner,Faction::Blight);
        s.advance_to(arrival);
        assert_eq!(s.systems[1].owner,Faction::Drummers);
    }
    #[test] fn autonomous_blight_schedules_only_causal_packets(){
        let mut s=StrategicState::default();
        s.advance_to(102_500.0);
        assert!(s.signals.iter().any(|p|p.kind==SignalKind::BlightSeed));
        for packet in &s.signals {
            if packet.kind!=SignalKind::BlightSeed{continue;}
            let distance=(s.systems[packet.source].pos-s.systems[packet.target].pos).len();
            assert!(packet.arrival_year>=packet.launch_year+distance/0.10);
            assert!((s.systems[packet.target].pos-ORION_OUTBREAK).len()
                <=Front::new(packet.arrival_year).reach_ly());
        }
    }
    #[test] fn convoy_escrow_and_arrival_reconcile_resources(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.systems[3].energy=300.0;
        s.systems[3].material=350.0;
        s.systems[5].energy=100.0;
        s.systems[5].material=100.0;
        let id=s.send_convoy(3,5).unwrap();
        assert_eq!(id,1);
        assert_eq!(s.systems[3].energy,235.0);
        assert_eq!(s.systems[3].material,290.0);
        let convoy=s.convoys[0].clone();
        assert!((convoy.arrival_year-convoy.launch_year-
            (s.systems[3].pos-s.systems[5].pos).len()/convoy.beta).abs()<1e-8);
        assert!(convoy.beta<=0.75);
        // Advance directly to just before and just after arrival;
        // production is zeroed to isolate cargo transfer.
        s.systems.iter_mut().for_each(|sys|sys.industry=0.0);
        s.advance_to(convoy.arrival_year-1.0);
        assert_eq!(s.convoys.len(),1);
        assert_eq!(s.systems[5].material,100.0);
        s.advance_to(convoy.arrival_year);
        assert_eq!(s.convoys.len(),0);
        assert_eq!(s.systems[5].material,160.0);
        assert_eq!(s.systems[5].energy,140.0);
        s.advance_to(convoy.arrival_year+500.0);
        assert_eq!(s.systems[5].material,160.0);
    }
    #[test] fn convoy_ownership_risk_keeps_cargo_out_of_enemy_hands(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.systems[3].energy=200.0;s.systems[3].material=200.0;
        s.send_convoy(3,5).unwrap();
        s.set_owner(5,Faction::Weena);
        s.systems.iter_mut().for_each(|sys|sys.industry=0.0);
        let eta=s.convoys[0].arrival_year;
        let before=s.systems[5].material;
        s.advance_to(eta);
        assert!(s.convoys.is_empty());
        assert_eq!(s.systems[5].material,before);
    }
    #[test] fn blight_holds_the_convoy_cordon_until_passage(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.set_owner(0,Faction::Drummers);
        let baseline=s.systems[0].material;
        s.send_convoy(1,0).unwrap();
        let eta=s.convoys[0].arrival_year;
        s.systems.iter_mut().for_each(|system|system.industry=0.0);
        s.advance_to(eta);
        assert!(s.convoys.is_empty());
        assert_eq!(s.systems[0].material,baseline);
        assert!(s.events.iter().any(|event|event.contains("interdicted")));
    }
    #[test] fn losing_a_guf_interception_leaves_beam_energy_intact(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        for i in [2,3,6] {
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        ready_array(&mut s,5);
        s.fire_guf(5).unwrap();
        assert_eq!(s.incoming_guf(5).unwrap().2,12.0);
        s.resolve_tactical(5,TacticalKind::GufDefense,
            TacticalOutcome::Defeat).unwrap();
        assert_eq!(s.systems[5].guf_ward,0.0);
        assert_eq!(s.incoming_guf(5).unwrap().2,12.0);
    }
    #[test] fn timely_guf_interception_absorbs_later_beams(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        for i in [2,3,6] {
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        ready_array(&mut s,5);
        s.systems[5].energy=200.0;s.systems[5].material=200.0;
        let eta=s.fire_guf(5).unwrap();
        assert_eq!(s.incoming_guf(5).unwrap().1,3);
        let initial_defense=s.systems[5].defense;
        s.resolve_tactical(5,TacticalKind::GufDefense,TacticalOutcome::Victory).unwrap();
        assert!(s.systems[5].guf_ward>0.0);
        s.systems.iter_mut().for_each(|sys|sys.industry=0.0);
        s.advance_to(eta);
        assert!(s.systems[5].defense>initial_defense-12.0);
        assert!(s.incoming_guf(5).is_none());
        assert_eq!(s.systems[5].guf_ward,0.0);
    }
    #[test] fn guf_intercept_requires_a_real_future_attack(){
        let mut s=StrategicState::default();
        assert!(s.resolve_tactical(5,TacticalKind::GufDefense,
            TacticalOutcome::Victory).is_err());
        assert_eq!(s.systems[5].guf_ward,0.0);
    }
    #[test] fn beam_aim_must_finish_before_guf_can_fire(){
        let mut s=StrategicState::default();s.blight_autonomy=false;
        for i in [2,3,6] {
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        s.array_charged=100.0;
        assert!(s.fire_guf(5).is_err());
        let ready=s.aim_guf(5).unwrap();
        assert!(ready>s.year);
        assert!(s.fire_guf(5).is_err());
        s.advance_to(ready-1.0);
        assert!(s.fire_guf(5).is_err());
        s.advance_to(ready);
        s.array_charged=100.0;
        assert!(s.fire_guf(5).is_ok());
    }
    #[test] fn displaced_target_escapes_locked_guf_beam(){
        let mut s=StrategicState::default();
        for i in [2,3,6]{
            s.set_owner(i,Faction::Weena);
            s.systems[i].array_node=true;s.systems[i].energy=500.0;
        }
        ready_array(&mut s,5);
        let eta=s.fire_guf(5).unwrap();
        let defense=s.systems[5].defense;
        // Move the star after launch, beyond narrow beam aperture.
        s.systems[5].pos.y+=300.0;
        s.advance_to(eta);
        assert_eq!(s.systems[5].defense,defense);
        assert!(s.events.iter().any(|e|e.contains("off-axis")));
    }
    #[test] fn insufficient_nodes_cannot_fire(){
        let mut s=StrategicState::default();s.array_charged=100.0;
        s.systems[6].array_node=true;
        assert!(s.fire_guf(5).is_err());
    }
    #[test] fn tactical_cooldown_prevents_instant_friendly_resource_farming(){
        let mut s=StrategicState::default();
        let before=s.systems[5].energy;
        assert!(s.resolve_tactical(5,TacticalKind::Intercept,
            TacticalOutcome::Victory).is_ok());
        let after=s.systems[5].energy;
        assert!(after>before);
        assert!(!s.tactical_ready(5));
        assert!(s.resolve_tactical(5,TacticalKind::Intercept,
            TacticalOutcome::Victory).is_err());
        assert_eq!(after,s.systems[5].energy);
        s.advance_to(100_250.0);
        assert!(s.tactical_ready(5));
    }
    #[test] fn weena_claims_are_delayed(){
        let mut s=StrategicState::new(100_000.0);
        let eta=s.claim_for_weena(5).unwrap();
        assert_eq!(s.systems[5].owner,Faction::Drummers);
        s.advance_to(eta);
        assert_eq!(s.systems[5].owner,Faction::Weena);
    }
    #[test] fn expansion_respects_orion_light_cone(){
        let mut s=StrategicState::new(100_000.0);
        s.set_owner(0,Faction::Blight);
        assert!(s.schedule_blight_seed(9,1).is_ok());
        let arrival=s.signals[0].arrival_year;
        assert!(arrival>=s.year+(s.systems[1].pos-s.systems[9].pos).len());
    }
    #[test] fn convoy_escort_victory_shields_convoy_and_defeat_destroys_it(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.set_owner(0,Faction::Drummers);
        s.send_convoy(1,0).unwrap();
        assert!(!s.convoys[0].escorted);
        let energy_before=s.systems[1].energy;
        s.resolve_tactical(1,TacticalKind::ConvoyEscort,TacticalOutcome::Victory).unwrap();
        assert!(s.convoys[0].escorted);
        assert!(s.systems[1].energy>energy_before);

        // Escort defeats the first interdiction check cleanly.
        let mut s2=StrategicState::default();
        s2.blight_autonomy=false;
        s2.systems[3].energy=300.0;s2.systems[3].material=300.0;
        s2.send_convoy(3,5).unwrap();
        s2.resolve_tactical(3,TacticalKind::ConvoyEscort,TacticalOutcome::Withdrawn).unwrap();
        assert_eq!(s2.convoys.len(),1);
        assert!(s2.tactical_ready(3));
        s2.resolve_tactical(3,TacticalKind::ConvoyEscort,TacticalOutcome::Defeat).unwrap();
        assert!(s2.convoys.is_empty());
        assert!(!s2.tactical_ready(3));
    }
    #[test] fn tactical_halo_strike_consumes_kugelblitz_once_and_never_on_withdrawal(){
        let mut s=StrategicState::default();
        s.blight_autonomy=false;
        s.systems[5].kugelblitz=1;
        s.systems[6].defense=6.0;
        s.systems[6].industry=3.0;
        // Withdrawal must NOT consume the kugelblitz or lock cooldown.
        s.resolve_tactical(6,TacticalKind::HaloStrike,TacticalOutcome::Withdrawn).unwrap();
        assert_eq!(s.systems[5].kugelblitz,1);
        assert_eq!(s.systems[6].defense,6.0);
        assert!(s.tactical_ready(6));
        // Victory consumes 1 kugelblitz and deals 4.5 defense + industry strike damage.
        s.resolve_tactical(6,TacticalKind::HaloStrike,TacticalOutcome::Victory).unwrap();
        assert_eq!(s.systems[5].kugelblitz,0);
        assert!((s.systems[6].defense-1.5).abs()<1e-9);
        assert!((s.systems[6].industry-2.25).abs()<1e-9);
        assert!(!s.tactical_ready(6));
    }
}
