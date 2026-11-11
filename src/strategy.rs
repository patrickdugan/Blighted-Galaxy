//! BG-STR-01: deterministic thousand-year strategic scenario layer.
//! Abstract resource units are game balance parameters, not astrophysical yields.
//! All messages and attacks propagate no faster than light in galactic coordinates.
use crate::overworld::{Front,ORION_OUTBREAK,WAYPOINTS};
use crate::V2;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Faction { Blight, Drummers, Weena, Independent }
impl Faction {
    pub fn label(self)->&'static str {
        match self { Self::Blight=>"JENKINS BLIGHT",Self::Drummers=>"DRUMMERS / REMNANT",
            Self::Weena=>"WEEN A / FORMOI",Self::Independent=>"INDEPENDENT" }
    }
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Project { Industry, Kugelblitz, GufNode }
#[derive(Clone,Debug)]
pub struct System {
    pub name:&'static str, pub pos:V2, pub owner:Faction,
    pub industry:f64, pub energy:f64, pub material:f64,
    pub kugelblitz:u32, pub array_node:bool, pub defense:f64,
}
#[derive(Clone,Debug)]
pub struct Signal {
    pub kind:SignalKind, pub source:usize, pub target:usize,
    pub arrival_year:f64, pub power:f64,
}
#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum SignalKind { BlightSeed, ArrayPulse }
#[derive(Clone,Debug)]
pub struct StrategicState {
    pub year:f64, pub systems:Vec<System>,
    pub signals:Vec<Signal>,pub events:Vec<String>,
    pub array_charged:f64,pub pulses:u32,
    pub historical:bool,
}
impl Default for StrategicState { fn default()->Self {Self::new(100_000.0)} }
impl StrategicState {
    pub fn new(year:f64)->Self {
        // Waypoint positions and project ownership are schematic campaign setup.
        let owners=[Faction::Independent,Faction::Drummers,Faction::Independent,
            Faction::Drummers,Faction::Drummers,Faction::Drummers,Faction::Weena,
            Faction::Independent,Faction::Independent];
        let systems=WAYPOINTS.iter().enumerate().map(|(i,w)|System {
            name:w.name,pos:w.pos,owner:owners[i],
            industry:if i==6 {3.0} else if matches!(i,3|4|5){2.0}else{0.5},
            energy:100.0,material:150.0,kugelblitz:0,
            array_node:false, defense:if i==5{4.0}else{2.0},
        }).collect();
        Self{year,systems,signals:vec![],events:vec![
            "STRATEGIC SCENARIO: industrial networks online".into()],
            array_charged:0.0,pulses:0,historical:false}
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
                if s.owner!=Faction::Weena{return Err("WEEN A CONTROL REQUIRED");}
                if s.industry<2.0||s.material<65.0||s.energy<75.0 {
                    return Err("ARRAY NODE NEEDS INDUSTRY / STORES");
                }
                s.material-=65.0;s.energy-=75.0;s.array_node=true;
                Ok(format!("{name}: synchronized Guf relay prepared"))
            },
        }.inspect(|_|{}).map(|msg|{self.events.push(msg.clone());msg})
    }
    pub fn set_owner(&mut self,site:usize,owner:Faction){
        if let Some(s)=self.systems.get_mut(site){s.owner=owner;}
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
            arrival_year:arrival,power:1.0});
        Ok(())
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
        self.signals.push(Signal{kind:SignalKind::ArrayPulse,
            source:emitters[0],target,arrival_year:arrival,
            power:emitters.len() as f64*4.0});
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
                        SignalKind::BlightSeed=>{
                            if (s.pos-ORION_OUTBREAK).len()<=Front::new(self.year).reach_ly(){
                                s.defense-=packet.power*2.0;
                                if s.defense<=0.0{s.owner=Faction::Blight;s.array_node=false;}
                                self.events.push(format!("BLIGHT wave reached {}",s.name));
                            }
                        }
                        SignalKind::ArrayPulse=>{
                            s.defense-=packet.power;
                            if s.defense<=0.0 {
                                s.industry=(s.industry*0.25).max(0.0);
                                s.array_node=false;
                            }
                            self.events.push(format!("GU F beam struck {} / defense {:.1}",s.name,s.defense));
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
    #[test] fn insufficient_nodes_cannot_fire(){
        let mut s=StrategicState::default();s.array_charged=100.0;
        s.systems[6].array_node=true;
        assert!(s.fire_guf(5).is_err());
    }
    #[test] fn expansion_respects_orion_light_cone(){
        let mut s=StrategicState::new(100_000.0);
        s.set_owner(0,Faction::Blight);
        assert!(s.schedule_blight_seed(0,5).is_ok());
        let arrival=s.signals[0].arrival_year;
        assert!(arrival>=s.year+(s.systems[5].pos-s.systems[0].pos).len());
    }
}
