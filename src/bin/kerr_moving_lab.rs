//! BG-PHY-08: moving equatorial Kerr observer and metric-aware target link.
 //! cargo run --release --bin kerr_moving_lab
 //!
 //! W toggles ZAMO / circular geodesic camera, F auto target acquisition.
 //! The photon is integrated backwards in Kerr coordinate time. Local
 //! ZAMO aberration and emitter-rest-frame target envelopes are included.
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::{KerrCamera,EquatorialNullPath};
use blighted_galaxy::kerr_moving_observer::{
    KerrMovingObserver,ObserverMotion,MovingOpticalObservation,
    observe_recorded_source,source_aligned_to_observer_ray,
};
use blighted_galaxy::kerr_optics::OrbitingCombatant;
use blighted_galaxy::causal::WorldlineHistory;
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::V2;
use macroquad::prelude::*;
use std::f64::consts::FRAC_PI_2;

fn conf()->Conf {
    Conf{window_title:"Blighted Galaxy • Moving Kerr Optical Observer".into(),
        window_width:1600,window_height:950,high_dpi:true,
        ..Default::default()}
}
struct Lab {
    observer:KerrMovingObserver,
    source:OrbitingCombatant,
    source_history:WorldlineHistory,
    observation:Option<MovingOpticalObservation>,
    manual_ray:Option<EquatorialNullPath>,
    time_s:f64,
    time_scale:f64,
    local_alpha_rg:f64,
    auto:bool,
    is_zamo:bool,
    scale:f32,
    accumulator:f64,
}
impl Lab {
    fn new(spin:f64,is_zamo:bool)->Self {
        let geometry=KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap();
        let camera=KerrCamera::new(geometry,85.0,FRAC_PI_2,36.0).unwrap();
        let motion=if is_zamo {ObserverMotion::Zamo}
            else {ObserverMotion::CircularGeodesic{direction:1}};
        let now=0.06;
        let temp=KerrMovingObserver::new(camera,motion,0.0).unwrap();
        // Place the moving observer near +x at the start of the lab,
        // so the same camera framing works for all spin configurations.
        let phi0=-temp.omega_rad_s*now;
        let observer=KerrMovingObserver::new(camera,motion,phi0).unwrap();
        let source=source_aligned_to_observer_ray(&observer,
            0.0,20.0,now,1,0.85).unwrap();
        let source_history=source.record(
            &camera,0.01,0.16,0.00008).unwrap();
        let mut lab=Self{
            observer,source,source_history,observation:None,
            manual_ray:None,time_s:now,time_scale:1.0,
            local_alpha_rg:0.0,auto:true,is_zamo,
            scale:7.1,accumulator:0.0,
        };
        lab.refresh();lab
    }
    fn refresh(&mut self) {
        if self.auto {
            self.observation=observe_recorded_source(
                &self.observer,&self.source_history,self.time_s,
                self.source.envelope_rg,4.0,81,17.0);
        } else {
            self.observation=None;
        }
        let alpha=self.observation.as_ref()
            .map(|x|x.alpha_camera_rg)
            .unwrap_or(self.local_alpha_rg);
        self.manual_ray=self.observer.ray_at(
            self.time_s,alpha,17.0,3600);
    }
    fn tick(&mut self,dt:f64) {
        let new_spin=if is_key_pressed(KeyCode::Key1){Some(0.0)}
            else if is_key_pressed(KeyCode::Key2){Some(0.75)}
            else if is_key_pressed(KeyCode::Key3){Some(-0.75)}
            else {None};
        if let Some(s)=new_spin {*self=Self::new(s,self.is_zamo);}
        if is_key_pressed(KeyCode::W) {
            *self=Self::new(self.observer.camera.geometry.spin,!self.is_zamo);
        }
        if is_key_pressed(KeyCode::R){
            *self=Self::new(self.observer.camera.geometry.spin,self.is_zamo);
        }
        if is_key_pressed(KeyCode::F) {
            self.auto=!self.auto;self.accumulator=10.0;
        }
        if is_key_pressed(KeyCode::Q) {
            self.time_scale=(self.time_scale/2.0).max(0.125);
        }
        if is_key_pressed(KeyCode::E) {
            self.time_scale=(self.time_scale*2.0).min(8.0);
        }
        if is_key_down(KeyCode::A) {
            self.local_alpha_rg=(self.local_alpha_rg-3.5*dt).max(-15.0);
            self.auto=false;self.accumulator=10.0;
        }
        if is_key_down(KeyCode::D) {
            self.local_alpha_rg=(self.local_alpha_rg+3.5*dt).min(15.0);
            self.auto=false;self.accumulator=10.0;
        }
        if is_key_down(KeyCode::Z) {self.scale=(self.scale/1.02).max(1.0);}
        if is_key_down(KeyCode::X) {self.scale=(self.scale*1.02).min(18.0);}
        let (_,scroll)=mouse_wheel();
        if scroll!=0.0 {
            self.scale=(self.scale*(1.0+scroll*0.08)).clamp(1.0,18.0);
        }
        self.time_s+=dt.min(0.05)*0.001*self.time_scale;
        if self.time_s>0.142 {
            self.time_s=0.04;
        }
        self.accumulator+=dt;
        if self.accumulator>0.11 {
            self.refresh();self.accumulator=0.0;
        }
    }
    fn pixel(&self,p:V2)->Vec2 {
        let rg=self.observer.camera.geometry.r_g_m;
        vec2(screen_width()*0.48,screen_height()*0.51)+
            vec2((p.x/rg) as f32,(p.y/rg) as f32)*self.scale
    }
    fn background(&self) {
        clear_background(Color::new(0.01,0.018,0.048,1.0));
        for i in 0..460_i64 {
            let x=((i*2833+1723)%5419) as f32/5419.0*screen_width();
            let y=((i*3617+2843)%5519) as f32/5519.0*screen_height();
            draw_circle(x,y,if i%29==0 {1.6} else {0.75},
                Color::new(0.72,0.8,1.0,0.17+0.04*(i%4) as f32));
        }
    }
    fn render(&self) {
        self.background();
        let geom=self.observer.camera.geometry;
        let centre=self.pixel(V2::default());
        let horizon=geom.outer_horizon_rg() as f32*self.scale;
        draw_circle(centre.x,centre.y,horizon,BLACK);
        draw_circle_lines(centre.x,centre.y,horizon,2.0,ORANGE);
        let isco=geom.isco_rg(1).unwrap();
        draw_circle_lines(centre.x,centre.y,
            isco as f32*self.scale,1.1,
            Color::new(0.65,0.4,1.0,0.65));
        draw_circle_lines(centre.x,centre.y,
            self.source.orbital_radius_rg as f32*self.scale,1.2,
            Color::new(0.8,0.55,0.25,0.22));
        let observer_pos=self.observer.position_m(self.time_s);
        let obs=self.pixel(observer_pos);
        draw_circle(obs.x,obs.y,9.0,WHITE);
        draw_circle_lines(obs.x,obs.y,17.0,1.2,WHITE);
        draw_text("ORBITING OBSERVER",obs.x-82.0,obs.y-28.0,16.0,WHITE);

        if let Some(ray)=self.manual_ray.as_ref() {
            for pair in ray.samples.windows(2) {
                let a=self.pixel(pair[0].chart_xy_m);
                let b=self.pixel(pair[1].chart_xy_m);
                draw_line(a.x,a.y,b.x,b.y,2.1,SKYBLUE);
            }
        }
        if let Some(sample)=self.source_history.sample_at(self.time_s) {
            let p=self.pixel(sample.pos_m);
            draw_circle(p.x,p.y,9.0,RED);
            draw_circle_lines(p.x,p.y,15.0,1.2,RED);
            draw_text("PRESENT SOURCE (DEBUG)",p.x+13.0,p.y-8.0,16.0,RED);
        }
        if let Some(acq)=self.observation.as_ref() {
            let p=self.pixel(acq.apparent_pos_m);
            draw_circle_lines(p.x,p.y,19.0,2.2,GOLD);
            draw_circle(p.x,p.y,5.0,GOLD);
            draw_text("RECEIVED EMISSION",p.x+16.0,p.y+4.0,16.0,YELLOW);
            let ray=self.pixel(acq.intersection_chart_m);
            draw_circle_lines(ray.x,ray.y,24.0,1.3,SKYBLUE);
        }
        draw_rectangle(13.0,13.0,940.0,210.0,
            Color::new(0.015,0.028,0.052,0.94));
        draw_text("BLIGHTED GALAXY / MOVING KERR OBSERVER",
            27.0,46.0,27.0,WHITE);
        draw_text("1 zero spin  2 +0.75  3 -0.75  |  W ZAMO/orbit  F auto acquire",
            27.0,79.0,19.0,SKYBLUE);
        draw_text("A/D manual photon   Q/E time rate   Z/X zoom   R reset",
            27.0,109.0,19.0,LIGHTGRAY);
        draw_text(&format!("SPIN a={:+.2}   β relative to ZAMO={:+.5}   dτ/dt={:.6}",
            geom.spin,self.observer.beta_zamo,
            self.observer.dproper_dcoordinate()),
            27.0,139.0,18.0,GOLD);
        draw_text(&format!("COORD TIME {:.6}s   OBSERVER φ={:.4}rad   MODE {}",
            self.time_s,self.observer.phi_at(self.time_s),
            if self.is_zamo {"ZAMO"} else {"CIRCULAR GEODESIC"}),
            27.0,169.0,18.0,LIGHTGRAY);
        draw_text(if self.auto {"LOCAL REST-FRAME WORLD TUBE  /  AUTO LOCK"}
            else {"MANUAL LOCAL CAMERA RAY  /  NO SENSOR LOCK"},
            27.0,200.0,17.0,GOLD);

        draw_rectangle(14.0,screen_height()-195.0,
            715.0,168.0,Color::new(0.015,0.028,0.052,0.94));
        let metrics=if let Some(acq)=&self.observation {
            [
                format!("EMISSION t     {:.9} s",acq.emitter_time_s),
                format!("CURVED DELAY   {:.6} ms",acq.delay_s*1000.0),
                format!("APPARENT SHIFT {:.3} r_g",
                    (acq.apparent_pos_m-acq.current_pos_m).len()/geom.r_g_m),
                format!("CAMERA α       {:+.5} r_g",acq.alpha_camera_rg),
                format!("RAY CHORDS     {}",acq.ray.samples.len()),
            ]
        } else {
            ["NO CURVED PHOTON/SOURCE MATCH".to_string(),
            "Historical source record is independently bounded".to_string(),
            "Manual probe: F enables source acquisition".to_string(),
            "Observer photons aberrated into the ZAMO frame".to_string(),
            "Metric-aware local target envelope".to_string()]
        };
        for (i,line) in metrics.iter().enumerate() {
            draw_text(line,28.0,screen_height()-168.0+i as f32*27.0,
                18.0,if i==1 {GOLD} else {LIGHTGRAY});
        }
        draw_text("Equatorial Kerr photons; rest-frame distance is local tangent metric, source history is causal.",
            20.0,screen_height()-12.0,16.0,GRAY);
    }
}
#[macroquad::main(conf)]
async fn main() {
    let mut lab=Lab::new(0.75,false);
    loop {
        if is_key_pressed(KeyCode::Escape){break;}
        lab.tick(get_frame_time() as f64);
        lab.render();
        next_frame().await;
    }
}
