//! SI-to-gameplay bridge. Scene state is metres, seconds, kg and joules.
 //! Physics: Newtonian stellar fields and relativistic thrust kinematics.
 //! Fictional systems: shields, high-power engines, and breakup impulses.
use crate::physics::{
    kinetic_energy_relativistic_j, lorentz_gamma, near_side_tidal_ratio_si,
    velocity_verlet_si, SiBody, AU_M, C_SI, G_SI, SOLAR_MASS_KG,
};
use crate::{swept_circle_hit, V2};
use crate::causal::{WorldlineHistory,WorldlineSample};
use crate::debris::{self, DebrisConfig};
use std::f64::consts::TAU;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Metres(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Seconds(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kilograms(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Joules(pub f64);
impl Metres {
    pub fn from_au(au:f64)->Self { assert!(au.is_finite()); Self(au*AU_M) }
    pub fn as_au(self)->f64 { self.0/AU_M }
}
impl Seconds {
    pub fn new(s:f64)->Self { assert!(s.is_finite()&&s>=0.0); Self(s) }
}
impl Kilograms {
    pub fn new(m:f64)->Self { assert!(m.is_finite()&&m>=0.0); Self(m) }
}
impl Joules {
    pub fn new(e:f64)->Self { assert!(e.is_finite()&&e>=0.0); Self(e) }
}

#[derive(Debug, Clone, Copy)]
pub struct CameraMeters {
    pub center_m:V2,
    pub metres_per_pixel:Metres,
}
impl CameraMeters {
    pub fn new(center_m:V2, metres_per_pixel:Metres)->Self {
        assert!(metres_per_pixel.0.is_finite()&&metres_per_pixel.0>0.0);
        Self{center_m,metres_per_pixel}
    }
    pub fn screen_of(&self,world_m:V2,width:f64,height:f64)->V2 {
        (world_m-self.center_m)/self.metres_per_pixel.0+V2::new(width*0.5,height*0.5)
    }
    pub fn world_of(&self,pixel:V2,width:f64,height:f64)->V2 {
        self.center_m+(pixel-V2::new(width*0.5,height*0.5))*self.metres_per_pixel.0
    }
    pub fn zoom(&mut self,factor:f64) {
        if factor.is_finite()&&factor>0.0 {
            self.metres_per_pixel.0=(self.metres_per_pixel.0/factor).clamp(1e-3,1e13);
        }
    }
}
pub fn gravitational_field_at(point:V2,fields:&[SiBody])->V2 {
    let mut a=V2::default();
    for body in fields {
        let delta=body.pos_m-point;
        let r2=delta.len_sq().max(body.radius_m.max(1.0).powi(2));
        a+=delta*(G_SI*body.mass_kg/(r2*r2.sqrt()));
    }
    a
}

/// Specific relativistic momentum u=gamma*v, converted back without an
/// unphysical velocity clamp. All components are in m/s.
pub fn velocity_from_specific_momentum(p:V2)->V2 {
    p/(p.len()/C_SI).hypot(1.0)
}
pub fn boost_velocity_m_s(ship:V2,local:V2)->Option<V2> {
    let speed=ship.len();
    let gamma=lorentz_gamma(speed)?;
    lorentz_gamma(local.len())?;
    if speed<=1e-9 {return Some(local);}
    let n=ship/speed;
    let parallel=n*local.dot(n);
    let transverse=local-parallel;
    let d=1.0+ship.dot(local)/C_SI.powi(2);
    if d<=0.0 {return None;}
    let combined=(parallel+ship+transverse/gamma)/d;
    if combined.len()>=C_SI {
        Some(combined.normalized()*(C_SI*(1.0-1e-13)))
    } else {Some(combined)}
}

#[derive(Debug, Clone, Copy)]
pub struct Craft {
    pub pos_m:V2,
    pub vel_m_s:V2,
    pub facing_rad:f64,
    pub dry_mass_kg:f64,
    pub cross_section_radius_m:f64,
    pub shield_energy_j:f64,
    pub shield_capacity_j:f64,
    pub proper_time_s:f64,
}
impl Craft {
    pub fn new(pos_m:V2,vel_m_s:V2,facing_rad:f64)->Self {
        assert!(lorentz_gamma(vel_m_s.len()).is_some());
        Self{pos_m,vel_m_s,facing_rad,
            dry_mass_kg:2.0e5,cross_section_radius_m:4.0e6,
            shield_energy_j:7.0e17,shield_capacity_j:7.0e17,proper_time_s:0.0}
    }
    pub fn heading(&self)->V2 {V2::from_angle(self.facing_rad)}
    pub fn alive(&self)->bool {self.shield_energy_j>0.0}
    pub fn gamma(&self)->f64 {
        lorentz_gamma(self.vel_m_s.len()).unwrap_or(f64::INFINITY)
    }
    pub fn tick(&mut self,dt:Seconds,thrust_m_s2:f64,fields:&[SiBody]) {
        if dt.0==0.0 {return;}
        assert!(thrust_m_s2.is_finite());
        let old=self.vel_m_s;
        let acc=self.heading()*thrust_m_s2+gravitational_field_at(self.pos_m,fields);
        let momentum=old*self.gamma()+acc*dt.0;
        self.vel_m_s=velocity_from_specific_momentum(momentum);
        let mid=(old+self.vel_m_s)*0.5;
        self.pos_m+=mid*dt.0;
        self.proper_time_s+=dt.0/
            lorentz_gamma(mid.len()).expect("subluminal midpoint");
    }
}
#[derive(Debug, Clone, Copy)]
pub struct KineticRound {
    pub pos_m:V2,
    pub vel_m_s:V2,
    pub mass_kg:f64,
    pub radius_m:f64,
    pub is_player:bool,
    pub lifetime_s:f64,
}
impl KineticRound {
    pub fn from_craft(craft:&Craft,is_player:bool)->Self {
        let muzzle=craft.heading()*(0.72*C_SI);
        Self{pos_m:craft.pos_m+craft.heading()*120.0,
            vel_m_s:boost_velocity_m_s(craft.vel_m_s,muzzle).unwrap(),
            mass_kg:2.0,radius_m:0.2,is_player,lifetime_s:12.0}
    }
}
pub fn relative_impact_energy_j(round:&KineticRound,target:&Craft)->Joules {
    let gp=lorentz_gamma(round.vel_m_s.len()).unwrap();
    let g_rel=(gp*target.gamma()*
        (1.0-round.vel_m_s.dot(target.vel_m_s)/C_SI.powi(2))).max(1.0);
    Joules::new((g_rel-1.0)*round.mass_kg*C_SI.powi(2))
}
#[derive(Debug, Clone, Copy)]
pub struct ImpactEvent {
    pub pos_m:V2,
    pub energy_j:Joules,
    pub target_is_player:bool,
    pub time_s:f64,
}
#[derive(Debug, Clone, Copy, Default)]
pub struct Control {
    pub player_thrust_m_s2:f64,
    pub player_facing_rad:f64,
    pub fire_player:bool,
}
pub struct SiJoust {
    pub sun:SiBody,
    pub player:Craft,
    pub opponent:Craft,
    pub rounds:Vec<KineticRound>,
    pub events:Vec<ImpactEvent>,
    pub time_s:f64,
    pub next_player_shot_s:f64,
    pub next_enemy_shot_s:f64,
    pub shots_fired:u64,
    pub player_history:WorldlineHistory,
    pub opponent_history:WorldlineHistory,
}
impl Default for SiJoust {fn default()->Self {Self::new()}}
impl SiJoust {
    pub fn new()->Self {
        let center=V2::new(AU_M,0.0);
        let player=Craft::new(center+V2::new(-250e6,-20e6),V2::default(),0.0);
        let opponent=Craft::new(center+V2::new(250e6,20e6),
            V2::new(-0.012*C_SI,0.0),std::f64::consts::PI);
        let mut player_history=WorldlineHistory::new(8192);
        let mut opponent_history=WorldlineHistory::new(8192);
        player_history.push(WorldlineSample{
            time_s:0.0,pos_m:player.pos_m,vel_m_s:player.vel_m_s,
        }).unwrap();
        opponent_history.push(WorldlineSample{
            time_s:0.0,pos_m:opponent.pos_m,vel_m_s:opponent.vel_m_s,
        }).unwrap();
        Self{
            sun:SiBody::new(SOLAR_MASS_KG,6.957e8,V2::default(),V2::default()),
            player,opponent,player_history,opponent_history,
            rounds:vec![],events:vec![],time_s:0.0,
            next_player_shot_s:0.0,next_enemy_shot_s:1.2,shots_fired:0,
        }
    }
    /// Optical position from a strictly recorded, bounded worldline.
    /// Returns None until a physical light signal has arrived from the
    /// recorded time window, rather than extrapolating nonexistent history.
    pub fn apparent_opponent(&self)->Option<(V2,f64)> {
        let e=self.opponent_history.observe(self.player.pos_m,self.time_s)?;
        Some((e.pos_m,e.delay_s))
    }
    pub fn last_impact(&self)->Option<&ImpactEvent> {self.events.last()}
    pub fn step(&mut self,dt:Seconds,input:Control) {
        assert!(dt.0>=0.0 && dt.0<=0.1);
        if dt.0==0.0 {return;}
        self.player.facing_rad=input.player_facing_rad;
        let player_before=self.player.pos_m;
        let opponent_before=self.opponent.pos_m;
        self.player.tick(dt,input.player_thrust_m_s2,&[self.sun]);
        let facing=(self.player.pos_m-self.opponent.pos_m).normalized();
        self.opponent.facing_rad=facing.y.atan2(facing.x);
        self.opponent.tick(dt,0.0,&[self.sun]);
        if input.fire_player && self.time_s>=self.next_player_shot_s && self.player.alive() {
            self.rounds.push(KineticRound::from_craft(&self.player,true));
            self.next_player_shot_s=self.time_s+0.25;
            self.shots_fired+=1;
        }
        if self.time_s>=self.next_enemy_shot_s && self.opponent.alive() {
            // Enemy fire control is limited to signals contained in its
            // recorded causal past. No target data means no targeted shot.
            if let Some(e)=self.player_history.observe(
                self.opponent.pos_m,self.time_s) {
                let aim=(e.pos_m-self.opponent.pos_m).normalized();
                self.opponent.facing_rad=aim.y.atan2(aim.x);
                self.rounds.push(KineticRound::from_craft(&self.opponent,false));
                self.next_enemy_shot_s=self.time_s+2.5;
            }
        }
        let mut impacts=Vec::new();
        for shot in &mut self.rounds {
            let before=shot.pos_m;
            shot.pos_m+=shot.vel_m_s*dt.0;
            shot.lifetime_s-=dt.0;
            let (b0,b1,r)=if shot.is_player {
                (opponent_before,self.opponent.pos_m,self.opponent.cross_section_radius_m)
            } else {(player_before,self.player.pos_m,self.player.cross_section_radius_m)};
            if let Some(f)=swept_circle_hit(before,shot.pos_m,shot.radius_m,b0,b1,r) {
                impacts.push((before+(shot.pos_m-before)*f,*shot));
                shot.lifetime_s=0.0;
            }
        }
        self.rounds.retain(|shot|shot.lifetime_s>0.0);
        for (pos,shot) in impacts {
            let target=if shot.is_player {&mut self.opponent} else {&mut self.player};
            let energy=relative_impact_energy_j(&shot,target);
            // 90% energy-to-shield coupling is a fictional game tuning value.
            target.shield_energy_j=(target.shield_energy_j-energy.0*0.9).max(0.0);
            self.events.push(ImpactEvent{pos_m:pos,energy_j:energy,
                target_is_player:!shot.is_player,time_s:self.time_s+dt.0*0.5});
        }
        if self.events.len()>80 {
            self.events.drain(0..self.events.len()-80);
        }
        self.time_s+=dt.0;
        self.player_history.push(WorldlineSample{
            time_s:self.time_s,pos_m:self.player.pos_m,
            vel_m_s:self.player.vel_m_s,
        }).expect("player worldline must remain causal and monotonic");
        self.opponent_history.push(WorldlineSample{
            time_s:self.time_s,pos_m:self.opponent.pos_m,
            vel_m_s:self.opponent.vel_m_s,
        }).expect("opponent worldline must remain causal and monotonic");
    }
}

/// Historical alias retained for renderer and encounter compatibility.
pub type TidalFragment = debris::FluidParcel;
#[derive(Debug, Clone, Copy)]
pub struct TidalDisruption {
    pub time_s:f64,
    pub ratio:f64,
    pub position_m:V2,
    pub original_mass_kg:f64,
}
pub struct SiTidalEncounter {
    pub black_hole:SiBody,
    pub white_dwarf:SiBody,
    pub fragments:Vec<TidalFragment>,
    pub disruption:Option<TidalDisruption>,
    pub time_s:f64,
    pub peak_tidal_ratio:f64,
}
impl Default for SiTidalEncounter {fn default()->Self {Self::new()}}
impl SiTidalEncounter {
    pub fn new()->Self {
        let m=10.0*SOLAR_MASS_KG;
        Self{
            black_hole:SiBody::new(m,2.0*G_SI*m/C_SI.powi(2),V2::default(),V2::default()),
            white_dwarf:SiBody::new(0.65*SOLAR_MASS_KG,9e6,
                V2::new(-1.45e8,2e7),V2::new(6e5,-8e4)),
            fragments:vec![],disruption:None,time_s:0.0,peak_tidal_ratio:0.0,
        }
    }
    pub fn total_visible_mass_kg(&self)->Kilograms {
        let debris=self.fragments.iter().map(|p|p.mass_kg).sum::<f64>();
        if self.disruption.is_some() {Kilograms(self.black_hole.mass_kg+debris)}
        else {Kilograms(self.black_hole.mass_kg+self.white_dwarf.mass_kg)}
    }
    pub fn step(&mut self,dt:Seconds) {
        assert!(dt.0>=0.0&&dt.0<=0.25);
        if self.disruption.is_none() {
            let mut pair=[self.black_hole,self.white_dwarf];
            velocity_verlet_si(&mut pair,dt.0,0.0);
            self.black_hole=pair[0];self.white_dwarf=pair[1];
            if let Some(q)=near_side_tidal_ratio_si(&self.black_hole,&self.white_dwarf) {
                self.peak_tidal_ratio=self.peak_tidal_ratio.max(q);
                if q>=1.0 {
                    self.disruption=Some(TidalDisruption{
                        time_s:self.time_s+dt.0,ratio:q,
                        position_m:self.white_dwarf.pos_m,
                        original_mass_kg:self.white_dwarf.mass_kg,
                    });
                    self.fragments=split_star_into_tracers(self.white_dwarf,96);
                }
            }
        } else {
            // Pairwise self-gravity + symmetric 2D pressure and black-hole
            // recoil; sound-speed closure approximates fluid behavior.
            debris::advance(&mut self.black_hole, &mut self.fragments,
                DebrisConfig::white_dwarf_proxy(),dt.0);
        }
        self.time_s+=dt.0;
    }
}
pub fn split_star_into_tracers(star:SiBody,count:usize)->Vec<TidalFragment> {
    assert!(count>=2&&count%2==0);
    let half=count/2;
    let mut out=Vec::with_capacity(count);
    for i in 0..half {
        let dir=V2::from_angle(TAU*(i as f64+0.5)/half as f64);
        let tangent=V2::new(-dir.y,dir.x);
        let offset=dir*star.radius_m*(0.15+0.8*(i%9) as f64/8.0);
        let kick=dir*25_000.0+tangent*9_000.0;
        out.push(TidalFragment{pos_m:star.pos_m+offset,
            vel_m_s:star.vel_m_s+kick,mass_kg:star.mass_kg/count as f64});
        out.push(TidalFragment{pos_m:star.pos_m-offset,
            vel_m_s:star.vel_m_s-kick,mass_kg:star.mass_kg/count as f64});
    }
    out
}
pub fn wicker_impact_energy_j()->Joules {
    Joules::new(kinetic_energy_relativistic_j(2.0,0.9993*C_SI).unwrap())
}
#[cfg(test)]
mod tests {
    use super::*;
    fn approx(a:f64,b:f64,rel:f64) {
        assert!((a-b).abs() <= rel*b.abs().max(1e-30),"{a:.12e} != {b:.12e}");
    }
    #[test]
    fn au_roundtrip() {
        approx(Metres::from_au(3.0).as_au(),3.0,1e-15);
    }
    #[test]
    fn camera_preserves_si_precision_at_one_au() {
        let c=V2::new(AU_M,AU_M);
        let v=CameraMeters::new(c,Metres(10.0));
        let p=v.screen_of(c+V2::new(120.0,-80.0),800.0,600.0);
        approx(p.x,412.0,1e-10);approx(p.y,292.0,1e-10);
        assert!((v.world_of(p,800.0,600.0)-c-V2::new(120.0,-80.0)).len()<1e-4);
    }
    #[test]
    fn boost_collinear() {
        let v=boost_velocity_m_s(V2::new(0.8*C_SI,0.0),
            V2::new(0.5*C_SI,0.0)).unwrap();
        approx(v.x/C_SI,1.3/1.4,1e-13);
    }
    #[test]
    fn boost_perpendicular() {
        let v=boost_velocity_m_s(V2::new(0.8*C_SI,0.0),
            V2::new(0.0,0.5*C_SI)).unwrap();
        approx(v.x/C_SI,0.8,1e-13);
        approx(v.y/C_SI,0.3,1e-13);
    }
    #[test]
    fn momentum_inversion_is_subluminal() {
        let v=velocity_from_specific_momentum(V2::new(1e13,-1e13));
        assert!(v.len()<C_SI);
        approx(v.x,1e13/(2.0_f64.sqrt()*1e13/C_SI).hypot(1.0),1e-14);
    }
    #[test]
    fn constant_speed_proper_time() {
        let mut ship=Craft::new(V2::default(),V2::new(0.8*C_SI,0.0),0.0);
        ship.tick(Seconds(10.0),0.0,&[]);
        approx(ship.proper_time_s,6.0,1e-14);
        approx(ship.pos_m.x,8.0*C_SI,1e-14);
    }
    #[test]
    fn high_thrust_remains_causal() {
        let mut ship=Craft::new(V2::default(),V2::default(),0.0);
        for _ in 0..300 {
            ship.tick(Seconds(0.1),100.0*C_SI,&[]);
            assert!(ship.vel_m_s.len()<C_SI);
        }
        assert!(ship.proper_time_s<30.0);
    }
    #[test]
    fn boosted_projectile_remains_causal() {
        let ship=Craft::new(V2::default(),V2::new(0.8*C_SI,0.0),0.0);
        let shot=KineticRound::from_craft(&ship,true);
        assert!(shot.vel_m_s.len()<C_SI);
        assert!(shot.vel_m_s.len()>ship.vel_m_s.len());
    }
    #[test]
    fn impact_energy_matches_stationary_target() {
        let target=Craft::new(V2::default(),V2::default(),0.0);
        let shot=KineticRound{pos_m:V2::default(),
            vel_m_s:V2::new(0.9993*C_SI,0.0),
            mass_kg:2.0,radius_m:0.2,is_player:true,lifetime_s:1.0};
        approx(relative_impact_energy_j(&shot,&target).0,
            wicker_impact_energy_j().0,1e-12);
    }
    #[test]
    fn fast_projectile_swept_collision() {
        let a=V2::new(-1.5e6,0.0);
        let b=a+V2::new(0.8*C_SI,0.0)*(1.0/60.0);
        assert!(swept_circle_hit(a,b,0.2,V2::default(),
            V2::default(),20.0).is_some());
        assert!(swept_circle_hit(V2::new(-1e6,100.0),
            V2::new(1e6,100.0),0.2,V2::default(),
            V2::default(),20.0).is_none());
    }
    #[test]
    fn rest_ship_without_forces_does_not_drift() {
        let mut ship=Craft::new(V2::default(),V2::default(),0.0);
        ship.tick(Seconds(0.05),0.0,&[]);
        assert_eq!(ship.pos_m,V2::default());
    }
    #[test]
    fn combat_damage_uses_si_impact_energy() {
        let mut game=SiJoust::new();
        game.opponent.pos_m=game.player.pos_m+V2::new(3e6,0.0);
        game.opponent.vel_m_s=V2::default();
        // Reset recorded past after test-fixture placement. Real gameplay
        // disallows discontinuous teleports between causal samples.
        game.opponent_history=WorldlineHistory::new(8192);
        game.opponent_history.push(WorldlineSample {
            time_s:0.0,pos_m:game.opponent.pos_m,
            vel_m_s:game.opponent.vel_m_s,
        }).unwrap();
        let shield=game.opponent.shield_energy_j;
        for _ in 0..20 {
            game.step(Seconds(1.0/120.0),Control{
                player_thrust_m_s2:0.0,
                player_facing_rad:0.0,fire_player:true,
            });
            if game.opponent.shield_energy_j<shield {break;}
        }
        assert!(game.opponent.shield_energy_j<shield);
        assert!(game.last_impact().unwrap().energy_j.0>0.0);
    }
    #[test]
    fn tracer_split_preserves_initial_mass_com_and_momentum() {
        let star=SiBody::new(0.65*SOLAR_MASS_KG,9e6,
            V2::new(1e8,2e8),V2::new(1.3e6,-7e5));
        let debris=split_star_into_tracers(star,96);
        let m=debris.iter().map(|x|x.mass_kg).sum::<f64>();
        let p=debris.iter().fold(V2::default(),|s,x|s+x.pos_m*x.mass_kg)/m;
        let v=debris.iter().fold(V2::default(),|s,x|s+x.vel_m_s*x.mass_kg)/m;
        approx(m,star.mass_kg,1e-14);
        assert!((p-star.pos_m).len()<1e-5);
        assert!((v-star.vel_m_s).len()<1e-6);
    }
    #[test]
    fn target_observation_uses_recorded_past_not_instant_prediction() {
        let mut game=SiJoust::new();
        assert!(game.apparent_opponent().is_none());
        for _ in 0..350 {
            game.step(Seconds(1.0/120.0),Control::default());
        }
        let (observed,delay)=game.apparent_opponent().unwrap();
        assert!(delay>1.0 && delay<3.0);
        assert!((observed-game.opponent.pos_m).len()>1.0e5);
        assert!(game.player_history.validate_segments());
    }
    #[test]
    fn debris_after_breakup_preserves_global_linear_momentum() {
        let mut game=SiTidalEncounter::new();
        for _ in 0..2500 {
            game.step(Seconds(0.1));
            if game.disruption.is_some() {break;}
        }
        assert!(game.disruption.is_some());
        let p0=debris::system_momentum_kg_m_s(&game.black_hole,&game.fragments);
        let before=game.black_hole.vel_m_s;
        for _ in 0..50 {game.step(Seconds(0.02));}
        let p1=debris::system_momentum_kg_m_s(&game.black_hole,&game.fragments);
        let scale=game.fragments[0].mass_kg*game.fragments[0].vel_m_s.len();
        assert!((p1-p0).len()/scale<1e-8);
        assert!((game.black_hole.vel_m_s-before).len()>0.0);
    }
    #[test]
    fn tidal_breakup_trigger_and_mass_accounting() {
        let mut game=SiTidalEncounter::new();
        let starting_mass=game.total_visible_mass_kg().0;
        for _ in 0..2500 {
            game.step(Seconds(0.1));
            if game.disruption.is_some() {break;}
        }
        let hit=game.disruption.expect("tidal threshold must trigger");
        assert!(hit.ratio>=1.0&&hit.time_s>0.0);
        assert_eq!(game.fragments.len(),96);
        approx(game.total_visible_mass_kg().0,starting_mass,1e-14);
        assert!((game.white_dwarf.pos_m-game.black_hole.pos_m).len()>
            game.black_hole.radius_m);
    }
}
