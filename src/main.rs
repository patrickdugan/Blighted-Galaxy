//! Main executable: SI-driven combat and tidal encounter lab.
 //! Legacy normalized prototype remains available via --bin arcade_legacy.
use blighted_galaxy::physics::{AU_M, C_SI, SOLAR_MASS_KG};
use blighted_galaxy::si_bridge::{
    CameraMeters, Control, Metres, Seconds, SiJoust, SiTidalEncounter,
};
use blighted_galaxy::V2;
use blighted_galaxy::overworld_nav::OverworldNav;
use macroquad::prelude::*;

fn config()->Conf {
    Conf {
        window_title:"Blighted Galaxy • Galactic Overworld".into(),
        window_width:1600,window_height:950,high_dpi:true,
        ..Default::default()
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode { Intercept, Tidal }
struct Effect {
    pos_m:V2,
    start_s:f64,
    magnitude:f64,
    tone:Color,
}
struct Game {
    mode:Mode,
    intercept:SiJoust,
    tidal:SiTidalEncounter,
    camera:CameraMeters,
    sim_rate:f64,
    accumulator:f64,
    manual_aim:bool,
    facing_rad:f64,
    effects:Vec<Effect>,
    displayed_latest_time_s:f64,
}
const DELTA:f64=1.0/120.0;
const TIDAL_STEP:f64=0.10;
impl Game {
    fn new(mode:Mode)->Self {
        let (center,scale)=match mode {
            Mode::Intercept=>(V2::new(AU_M,0.0),Metres(550_000.0)),
            Mode::Tidal=>(V2::default(),Metres(350_000.0)),
        };
        Self {mode,intercept:SiJoust::new(),tidal:SiTidalEncounter::new(),
            camera:CameraMeters::new(center,scale),
            sim_rate:1.0,accumulator:0.0,manual_aim:false,
            facing_rad:0.0,effects:Vec::new(),displayed_latest_time_s:-1.0}
    }
    fn screen(&self,pos:V2)->Vec2 {
        let p=self.camera.screen_of(pos,screen_width() as f64,screen_height() as f64);
        vec2(p.x as f32,p.y as f32)
    }
    fn current_time(&self)->f64 {
        match self.mode {
            Mode::Intercept=>self.intercept.time_s,
            Mode::Tidal=>self.tidal.time_s,
        }
    }
    fn input(&mut self) {
        if is_key_pressed(KeyCode::Key1) { *self=Self::new(Mode::Intercept); }
        if is_key_pressed(KeyCode::Key2) { *self=Self::new(Mode::Tidal); }
        if is_key_pressed(KeyCode::R) {
            let mode=self.mode;*self=Self::new(mode);
        }
        if is_key_pressed(KeyCode::Q) {self.sim_rate=(self.sim_rate/2.0).max(0.125);}
        if is_key_pressed(KeyCode::E) {self.sim_rate=(self.sim_rate*2.0).min(64.0);}
        if is_key_pressed(KeyCode::M) {self.manual_aim=!self.manual_aim;}
        let (_,scroll)=mouse_wheel();
        if scroll!=0.0 {self.camera.zoom((1.0+scroll as f64*0.1).clamp(0.2,3.0));}
        if is_key_down(KeyCode::X) {self.camera.zoom(1.03);}
        if is_key_down(KeyCode::Z) {self.camera.zoom(1.0/1.03);}
        if self.mode==Mode::Intercept {
            if self.manual_aim {
                if is_key_down(KeyCode::A) {self.facing_rad-=0.045;}
                if is_key_down(KeyCode::D) {self.facing_rad+=0.045;}
            } else {
                let (mx,my)=mouse_position();
                let cursor=self.camera.world_of(V2::new(mx as f64,my as f64),
                    screen_width() as f64,screen_height() as f64);
                let toward=cursor-self.intercept.player.pos_m;
                if toward.len()>0.0 {self.facing_rad=toward.y.atan2(toward.x);}
            }
        }
    }
    fn advance(&mut self,real_dt:f64) {
        let step=if self.mode==Mode::Intercept {DELTA} else {TIDAL_STEP};
        self.accumulator=(self.accumulator+real_dt.min(0.05)*self.sim_rate).min(8.0);
        let mut cycles=0;
        while self.accumulator>=step && cycles<400 {
            match self.mode {
                Mode::Intercept=>{
                    let thrust=if is_key_down(KeyCode::W) {0.02*C_SI}
                        else if is_key_down(KeyCode::S) {-0.01*C_SI} else {0.0};
                    self.intercept.step(Seconds(step),Control {
                        player_thrust_m_s2:thrust,
                        player_facing_rad:self.facing_rad,
                        fire_player:is_key_down(KeyCode::Space),
                    });
                }
                Mode::Tidal=>self.tidal.step(Seconds(step)),
            }
            self.accumulator-=step;
            cycles+=1;
        }
        if self.mode==Mode::Intercept {
            for evt in self.intercept.events.iter().filter(|x|x.time_s>self.displayed_latest_time_s) {
                self.effects.push(Effect{
                    pos_m:evt.pos_m,start_s:evt.time_s,
                    magnitude:evt.energy_j.0/1e17,
                    tone:if evt.target_is_player {RED} else {SKYBLUE},
                });
            }
            if let Some(last)=self.intercept.events.last() {
                self.displayed_latest_time_s=last.time_s;
            }
            let camera_target=(self.intercept.player.pos_m+self.intercept.opponent.pos_m)*0.5;
            self.camera.center_m+=(camera_target-self.camera.center_m)*0.08;
        } else {
            self.camera.center_m+=(self.tidal.black_hole.pos_m-self.camera.center_m)*0.06;
        }
        let now=self.current_time();
        self.effects.retain(|e|now-e.start_s<2.4);
    }
    fn stars(&self) {
        for i in 0..370_u32 {
            let x=((i as i64*2207+3181)%8000) as f32/8000.0*screen_width();
            let y=((i as i64*3911+701)%8000) as f32/8000.0*screen_height();
            let sz=if i%31==0 {1.7} else {0.8};
            draw_circle(x,y,sz,Color::new(0.45,0.58,0.79,0.2+(i%5) as f32*0.06));
        }
    }
    fn glow(&self,p:Vec2,r:f32,c:Color) {
        for layer in (1..=5).rev() {
            let k=layer as f32;
            draw_circle(p.x,p.y,r*(1.0+k*0.48),
                Color::new(c.r,c.g,c.b,0.032));
        }
    }
    fn actor(&self,pos_m:V2,theta:f64,c:Color,size:f32) {
        let p=self.screen(pos_m);
        let h=vec2(theta.cos() as f32,theta.sin() as f32);
        let w=vec2(-h.y,h.x);
        self.glow(p,size,c);
        draw_triangle(p+h*size*1.4,p-h*size+w*size*0.80,
            p-h*size-w*size*0.80,c);
        draw_circle_lines(p.x,p.y,size*1.5,1.0,c);
    }
    fn labeled(&self,title:&str,value:&str,y:f32,c:Color) {
        draw_text(title,24.0,y,20.0,LIGHTGRAY);
        draw_text(value,245.0,y,20.0,c);
    }
    fn paint_intercept(&self) {
        let player=&self.intercept.player;
        let foe=&self.intercept.opponent;
        let a=self.screen(player.pos_m);
        let b=self.screen(foe.pos_m);
        let shield_player=(player.cross_section_radius_m/
            self.camera.metres_per_pixel.0) as f32;
        let shield_foe=(foe.cross_section_radius_m/
            self.camera.metres_per_pixel.0) as f32;
        draw_circle_lines(a.x,a.y,shield_player.max(9.0),1.4,SKYBLUE);
        draw_circle_lines(b.x,b.y,shield_foe.max(9.0),1.4,ORANGE);
        self.actor(player.pos_m,player.facing_rad,SKYBLUE,16.0);
        self.actor(foe.pos_m,foe.facing_rad,ORANGE,19.0);
        if let Some((apparent,delay))=self.intercept.apparent_opponent() {
            let ghost=self.screen(apparent);
            draw_circle_lines(ghost.x,ghost.y,12.0,1.8,GOLD);
            draw_line(ghost.x,ghost.y,b.x,b.y,1.2,
                Color::new(1.0,0.86,0.38,0.54));
            draw_text("LIGHT-CONE OBSERVATION",ghost.x+17.0,ghost.y+5.0,15.0,GOLD);
            draw_text(&format!("T - {:.2} s",delay),ghost.x+17.0,ghost.y+24.0,15.0,YELLOW);
        }
        for round in &self.intercept.rounds {
            let p=self.screen(round.pos_m);
            let tail=self.screen(round.pos_m-round.vel_m_s*0.055);
            let tone=if round.is_player {SKYBLUE} else {RED};
            draw_line(tail.x,tail.y,p.x,p.y,2.5,tone);
            draw_circle(p.x,p.y,2.4,WHITE);
        }
        let now=self.current_time();
        for effect in &self.effects {
            let p=self.screen(effect.pos_m);
            let age=((now-effect.start_s)/2.4).clamp(0.0,1.0) as f32;
            let radius=(30.0+effect.magnitude.sqrt() as f32*34.0)*(0.3+age*2.2);
            draw_circle_lines(p.x,p.y,radius,2.3,
                Color::new(effect.tone.r,effect.tone.g,effect.tone.b,1.0-age));
        }
        let px=self.camera.metres_per_pixel.0;
        let distance=(foe.pos_m-player.pos_m).len();
        let light_delay=distance/C_SI;
        let last=self.intercept.last_impact();
        self.labeled("RANGE", &format!("{:.0} Mm",distance/1e6),245.0,WHITE);
        self.labeled("LIGHT TIME", &format!("{:.3} s",light_delay),274.0,GOLD);
        self.labeled("SHIP SPEED", &format!("{:.5} c",player.vel_m_s.len()/C_SI),303.0,SKYBLUE);
        self.labeled("LORENTZ", &format!("gamma = {:.5}",player.gamma()),332.0,SKYBLUE);
        self.labeled("PROPER TIME", &format!("{:.3} s",player.proper_time_s),361.0,WHITE);
        self.labeled("PLAYER SHIELD", &format!("{:.1}%",
            100.0*player.shield_energy_j/player.shield_capacity_j),390.0,SKYBLUE);
        self.labeled("ENEMY SHIELD", &format!("{:.1}%",
            100.0*foe.shield_energy_j/foe.shield_capacity_j),419.0,ORANGE);
        self.labeled("RKM SHOTS", &self.intercept.shots_fired.to_string(),448.0,LIGHTGRAY);
        if let Some(evt)=last {
            self.labeled("IMPACT ENERGY", &format!("{:.3e} J",evt.energy_j.0),
                477.0,YELLOW);
        }
        self.labeled("SCALE", &format!("{:.3} km/px",px/1000.0),
            506.0,LIGHTGRAY);
        if !foe.alive() {
            draw_text("TARGET SHIELD DISABLED",screen_width()*0.32,screen_height()*0.8,
                39.0,SKYBLUE);
        }
        if !player.alive() {
            draw_text("PLAYER SHIELD DISABLED",screen_width()*0.32,screen_height()*0.8,
                39.0,RED);
        }
        draw_text("SHIELD RADII ARE FICTIONAL DEFENSE ENVELOPES; SHIP MARKERS EXAGGERATED",
            25.0,screen_height()-17.0,17.0,GRAY);
    }
    fn paint_tidal(&self) {
        let b=&self.tidal.black_hole;
        let star=&self.tidal.white_dwarf;
        let bh=self.screen(b.pos_m);
        self.glow(bh,19.0,ORANGE);
        draw_circle(bh.x,bh.y,21.0,ORANGE);
        draw_circle(bh.x,bh.y,16.0,BLACK);
        draw_circle_lines(bh.x,bh.y,24.0,2.0,YELLOW);
        draw_text("10 M_SUN BLACK HOLE",bh.x+32.0,bh.y+7.0,18.0,ORANGE);
        let white=self.screen(star.pos_m);
        if self.tidal.disruption.is_none() {
            self.glow(white,13.0,SKYBLUE);
            draw_circle(white.x,white.y,14.0,LIGHTGRAY);
            draw_text("0.65 M_SUN WHITE DWARF",white.x+19.0,white.y+7.0,17.0,WHITE);
        } else {
            let event=self.tidal.disruption.unwrap();
            let locus=self.screen(event.position_m);
            draw_circle_lines(locus.x,locus.y,
                (self.tidal.time_s-event.time_s).max(0.0) as f32*18.0+15.0,
                1.5,Color::new(1.0,0.45,0.1,0.6));
            for fragment in &self.tidal.fragments {
                let p=self.screen(fragment.pos_m);
                if p.x>0.0&&p.y>0.0&&p.x<screen_width()&&p.y<screen_height() {
                    draw_circle(p.x,p.y,2.3,Color::new(0.56,0.82,1.0,0.86));
                }
            }
        }
        let dist=(b.pos_m-star.pos_m).len();
        let q=blighted_galaxy::physics::near_side_tidal_ratio_si(b,star)
            .unwrap_or(f64::NAN);
        self.labeled("SEPARATION",&format!("{:.2} Mm",dist/1e6),245.0,WHITE);
        self.labeled("ORBIT SPEED",&format!("{:.1} km/s",star.vel_m_s.len()/1000.0),
            274.0,WHITE);
        self.labeled("TIDAL / SELF G",&format!("{:.4}",q),303.0,YELLOW);
        self.labeled("MAX TIDAL RATIO",&format!("{:.4}",self.tidal.peak_tidal_ratio),
            332.0,YELLOW);
        self.labeled("ELAPSED",&format!("{:.2} s",self.tidal.time_s),361.0,LIGHTGRAY);
        self.labeled("MASS TRACED",&format!("{:.4} M_SUN",
            self.tidal.total_visible_mass_kg().0/SOLAR_MASS_KG),
            390.0,SKYBLUE);
        if let Some(event)=self.tidal.disruption {
            self.labeled("DISRUPTION",&format!("at {:.2} s; q = {:.3}",
                event.time_s,event.ratio),419.0,ORANGE);
            self.labeled("DEBRIS", &format!("{} tracers",self.tidal.fragments.len()),
                448.0,WHITE);
            draw_text("PROXY BREAKUP: DEBRIS IS BALLISTIC, WITHOUT STELLAR HYDRODYNAMICS",
                20.0,screen_height()-17.0,16.0,GRAY);
        } else {
            draw_text("NEAR-SIDE TIDAL ACCELERATION / SELF-GRAVITY; FAILURE GATE AT 1",
                20.0,screen_height()-17.0,17.0,GRAY);
        }
    }
    fn draw(&self) {
        clear_background(Color::new(0.010,0.015,0.042,1.0));
        self.stars();
        match self.mode {
            Mode::Intercept=>self.paint_intercept(),
            Mode::Tidal=>self.paint_tidal(),
        }
        draw_rectangle(13.0,12.0,800.0,207.0,
            Color::new(0.022,0.033,0.058,0.92));
        draw_rectangle_lines(13.0,12.0,800.0,207.0,1.0,DARKGRAY);
        draw_text("BLIGHTED GALAXY / SI GAMEPLAY BRIDGE",23.0,45.0,27.0,WHITE);
        let title=match self.mode {
            Mode::Intercept=>"1  RELATIVISTIC RKM INTERCEPT • 2kg PROJECTILES",
            Mode::Tidal=>"2  WHITE-DWARF TIDAL BREAKUP • 10 SOLAR MASS HOLE",
        };
        draw_text(title,23.0,80.0,22.0,SKYBLUE);
        draw_text("1 / 2 SCENARIOS   3 / TAB GALAXY   Q / E TIME   Z / X / WHEEL ZOOM",
            23.0,113.0,19.0,LIGHTGRAY);
        if self.mode==Mode::Intercept {
            draw_text("MOUSE aim • M keyboard aim • A/D rotate • W/S thrust • SPACE fire",
                23.0,145.0,18.0,LIGHTGRAY);
        } else {
            draw_text("Observe calibrated freefall and physically scaled tidal ratio",
                23.0,145.0,18.0,LIGHTGRAY);
        }
        draw_text(&format!("SIM RATE {:5.3}x | COORDINATE TIME {:8.3} s | R RESET | ESC QUIT",
            self.sim_rate,self.current_time()),23.0,179.0,18.0,GOLD);
        draw_text("SI kinematics + Newtonian gravity • exotic shields and breakup are game rules",
            23.0,204.0,16.0,GRAY);
    }
}

#[macroquad::main(config)]
async fn main() {
    let mut game=Game::new(Mode::Intercept);
    let mut overworld=OverworldNav::new();
    // The campaign opens on its galactic navigation layer. Existing
    // physics encounters remain directly reachable with keys 1 and 2.
    // --encounter keeps the previous startup for dedicated combat tests.
    let encounter_start=std::env::args().any(|a|a=="--encounter");
    let mut galaxy_mode=!encounter_start;
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        let dt=get_frame_time() as f64;
        if galaxy_mode {
            if is_key_pressed(KeyCode::Key1){
                game=Game::new(Mode::Intercept);
                galaxy_mode=false;
            }else if is_key_pressed(KeyCode::Key2){
                game=Game::new(Mode::Tidal);
                galaxy_mode=false;
            }else{
                overworld.update(dt);
            }
        }else if is_key_pressed(KeyCode::Key3)||
                is_key_pressed(KeyCode::Tab) {
            galaxy_mode=true;
        }else{
            game.input();
            game.advance(dt);
        }
        if galaxy_mode {overworld.draw();}
        else {game.draw();}
        next_frame().await;
    }
}
