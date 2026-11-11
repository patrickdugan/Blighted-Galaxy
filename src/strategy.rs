//! BG-STR-01: deterministic thousand-year strategic scenario layer.
//! Abstract resource units are game balance parameters, not astrophysical yields.
//! All messages and attacks propagate no faster than light in galactic coordinates.
use crate::overworld::{Front,ORION_OUTBREAK,WAYPOINTS,SectorStatus};
use crate::V2;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Faction { Blight, Drummers, Weena, Independent }
impl Faction {
    pub fn label(self)->&'static str {
        match self { Self::Blight=>"JENKINS BLIGHT",Self::Drummers=>"DRUMMERS / REMNANT",
            Self::Weena=>"WEENA / FORMOI",Self::Independent=>"INDEPENDENT" }
    }
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Project { Industry, Kugelblitz, GufNode }
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum TacticalKind { Intercept, Tidal }
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum TacticalOutcome { Victory, Defeat, Withdrawn }
#[derive(Clone,Debug)]
pub struct System {
    pub name:&'static str, pub pos:V2, pub owner:Faction,
    pub industry:f64, pub energy:f64, pub material:f64,
    pub kugelblitz:u32, pub array_node:bool, pub defense:f64,
}
#[derive(Clone,Debug)]
pub struct Signal {
    pub kind:SignalKind, pub source:usize, pub target:usize,
    pub launch_year:f64, pub arrival_year:f64, pub power:f64,
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum SignalKind { BlightSeed, ArrayPulse, WeenaClaim }
#[derive(Clone,Debug)]
pub struct StrategicState {
    pub year:f64, pub systems:Vec<System>,
    pub signals:Vec<Signal>,pub events:Vec<String>,
    pub array_charged:f64,pub pulses:u32,
    pub historical:bool,
    /// Incremented whenever a world changes its territorial footprint.
    pub territory_revision:u64,
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
            energy:100.0,material:150.0,kugelblitz:0,
            array_node:false, defense:if i==5{4.0}else{2.0},
        }).collect();
        systems.push(System{name:"ORION SEED",pos:ORION_OUTBREAK,
            owner:Faction::Blight,industry:4.0,energy:250.0,material:250.0,
            kugelblitz:0,array_node:false,defense:25.0});
        Self{year,systems,signals:vec![],events:vec![
            "STRATEGIC SCENARIO: industrial networks online".into()],
            array_charged:0.0,pulses:0,historical:false,
            territory_revision:0}
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
    pub fn resolve_tactical(&mut self,site:usize,kind:TacticalKind,
        outcome:TacticalOutcome)->Result<String,&'static str>{
        let s=self.systems.get_mut(site).ok_or("UNKNOWN TACTICAL SYSTEM")?;
        let name=s.name;
        let descriptor=match kind {
            TacticalKind::Intercept=>"RKM INTERCEPT",
            TacticalKind::Tidal=>"WHITE-DWARF INTERDICTION",
        };
        let message=match outcome {
            TacticalOutcome::Withdrawn=>format!("{name}: {descriptor} disengaged / no strategic change"),
            TacticalOutcome::Victory=>{
                if matches!(s.owner,Faction::Blight|Faction::Weena){
                    let damage=if kind==TacticalKind::Intercept{3.0}else{2.0};
                    s.defense-=damage;
                    s.material=(s.material-25.0).max(0.0);
                    if s.defense<=0.0{
                        s.owner=Faction::Drummers;s.array_node=false;
                        s.industry=(s.industry*0.6).max(0.25);
                        s.defense=1.0;
                        format!("{name}: {descriptor} victory / hostile control broken")
                    }else{
                        format!("{name}: {descriptor} victory / enemy defense {:.1}",s.defense)
                    }
                }else{
                    s.defense=(s.defense+1.25).min(18.0);
                    s.energy=(s.energy+20.0).min(1000.0);
                    format!("{name}: {descriptor} victory / defenses reinforced")
                }
            },
            TacticalOutcome::Defeat=>{
                s.defense=(s.defense-1.5).max(0.0);
                s.material=(s.material-18.0).max(0.0);
                if s.defense<=0.0 && s.owner==Faction::Drummers &&
                    (s.pos-ORION_OUTBREAK).len()<=Front::new(self.year).reach_ly(){
                    s.owner=Faction::Blight;s.array_node=false;
                    s.defense=1.0;
                    format!("{name}: {descriptor} defeat / site lost to Blight")
                }else{
                    format!("{name}: {descriptor} defeat / defenses damaged")
                }
            },
        };
        if outcome!=TacticalOutcome::Withdrawn{self.territory_revision+=1;}
        self.events.push(message.clone());
        Ok(message)
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
        }.map(|msg|{self.territory_revision+=1;self.events.push(msg.clone());msg})
    }
    pub fn set_owner(&mut self,site:usize,owner:Faction){
        if let Some(s)=self.systems.get_mut(site){s.owner=owner;self.territory_revision+=1;}
    }
    pub fn schedule_blight_seed(&mut self,source:usize,target:usize)->Result<(),&'static str>{
        let (a,b)=(self.systems.get(source).ok_or("INVALID SOURCE")?,
                   self.systems.get(target).ok_or("INVALID TARGET")?);
        if a.owner!=Faction::Blight{return Err("SOURCE IS NOT BLIGHT");}
        // Blight tendrils have effective expansion 0.10c; extra delay
        // encodes resistance. This also enforces the causal front geometry.
        let delay=(b.pos-a.pos).len()/0.10+1_500.0;
        let arrival=self.year+delay;
        if (b.pos-ORION_OUTBREAK).len()>Front::new(arrival).reach_ly(){
            return Err("OUTSIDE ORION CAUSAL ENVELOPE");
        }
        self.signals.push(Signal{kind:SignalKind::BlightSeed,source,target,
            launch_year:self.year,arrival_year:arrival,power:1.0});
        Ok(())
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
        self.systems[source].energy-=80.0;
        self.signals.push(Signal{kind:SignalKind::WeenaClaim,source,target,
            arrival_year:arrival,power:1.0});
        self.events.push(format!("WEENA influence relay sent to {}",self.systems[target].name));
        Ok(arrival)
    }
    /// Requires three independent industrial array nodes; power travels at c.
    /// The remaining beam is game-scale and may be intercepted on arrival.
    pub fn fire_guf(&mut self,target:usize)->Result<f64,&'static str>{
        let destination=self.systems.get(target).ok_or("INVALID TARGET")?;
        let emitters:Vec<usize>=self.systems.iter().enumerate()
            .filter(|(_,s)|s.owner==Faction::Weena&&s.array_node&&s.energy>=50.0)
            .map(|(i,_)|i).collect();
        if emitters.len()<3{return Err("THREE ACTIVE GUF NODES REQUIRED");}
        if self.array_charged<100.0{return Err("ARRAY CHARGING");}
        let target_pos=destination.pos;
        // Each emitter's light cone is independent. A distributed strike
        // reaches the destination at the latest packet arrival.
        let arrival=emitters.iter().map(|&i|
            self.year+(self.systems[i].pos-target_pos).len())
            .fold(self.year,f64::max);
        for &i in &emitters {self.systems[i].energy-=50.0;}
        self.array_charged=0.0;self.pulses+=1;
        // Delayed launches align all beam arrivals at the same target time.
        // Each separate flight is visible in the overworld and resolves once.
        for &source in &emitters {
            let light_time=(self.systems[source].pos-target_pos).len();
            self.signals.push(Signal{kind:SignalKind::ArrayPulse,
                source,target,launch_year:arrival-light_time,
                arrival_year:arrival,power:4.0});
        }
        self.events.push(format!("CHAMBER OF GUF fired; target arrival T+{arrival:.0} yr"));
        Ok(arrival)
    }
    pub fn advance_to(&mut self,year:f64){
        if !year.is_finite()||year<=self.year{return;}
        let end=year.min(200_000.0);
        while self.year<end {
            let next=(self.year+250.0).min(end);
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
                if let Some(s)=self.systems.get_mut(packet.target){
                    match packet.kind {
                        SignalKind::WeenaClaim=>{
                            if s.owner!=Faction::Blight {
                                s.owner=Faction::Weena;s.array_node=false;
                                self.territory_revision+=1;
                                self.events.push(format!("WEENA recruited {}",s.name));
                            }
                        }
                        SignalKind::BlightSeed=>{
                            if (s.pos-ORION_OUTBREAK).len()<=Front::new(self.year).reach_ly(){
                                s.defense-=packet.power*2.0;
                                if s.defense<=0.0{s.owner=Faction::Blight;s.array_node=false;}
                                self.territory_revision+=1;
                                self.events.push(format!("BLIGHT wave reached {}",s.name));
                            }
                        }
                        SignalKind::ArrayPulse=>{
                            s.defense-=packet.power;
                            self.territory_revision+=1;
                            if s.defense<=0.0 {
                                s.industry=(s.industry*0.25).max(0.0);
                                s.array_node=false;
                            }
                            self.events.push(format!("GUF beam struck {} / defense {:.1}",s.name,s.defense));
                        }
                    }
                }
            }
            self.signals=pending;
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
        s.array_charged=100.0;
        let latest=[2,3,6].iter().map(|&i|(s.systems[i].pos-s.systems[5].pos).len())
            .fold(0.0_f64,f64::max);
        let arrival=s.fire_guf(5).unwrap();
        assert!((arrival-100_000.0-latest).abs()<1e-6);
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
        s.array_charged=100.0;
        let arrival=s.fire_guf(5).unwrap();
        assert_eq!(s.signals.len(),3);
        for packet in &s.signals {
            assert_eq!(packet.kind,SignalKind::ArrayPulse);
            assert!((packet.arrival_year-packet.launch_year-
                (s.systems[packet.source].pos-s.systems[5].pos).len()).abs()<1e-8);
            assert!((packet.arrival_year-arrival).abs()<1e-8);
        }
    }
    #[test] fn insufficient_nodes_cannot_fire(){
        let mut s=StrategicState::default();s.array_charged=100.0;
        s.systems[6].array_node=true;
        assert!(s.fire_guf(5).is_err());
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
}
