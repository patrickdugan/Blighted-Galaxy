//! Live curved-Kerr optical acquisition of a moving posthuman target.
 //! Run: cargo run --release --bin kerr_optical_lab
 //!
 //! Equatorial ZAMO observer and a small source on a circular Kerr timelike
 //! geodesic. A reverse-propagated null ray intersects an interpolated
 //! *recorded* source worldline. This is a sensor/aiming display, not a
 //! physically correct outgoing plasma weapon or arbitrary worldtube metric.
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::{KerrCamera,EquatorialNullPath};
use blighted_galaxy::kerr_optics::{
    OrbitingCombatant,OpticalAcquisition,
    observe_combatant,source_aligned_to_past_ray,
};
use blighted_galaxy::causal::WorldlineHistory;
use blighted_galaxy::physics::{SOLAR_MASS_KG,C_SI};
use blighted_galaxy::V2;
use macroquad::prelude::*;
use std::f64::consts::FRAC_PI_2;

fn settings()->Conf {
    Conf{
        window_title:"Blighted Galaxy • Kerr Optical Intercept".into(),
        window_width:1600,window_height:920,high_dpi:true,
        ..Default::default()
    }
}
struct Lab {
    camera:KerrCamera,
    combatant:OrbitingCombatant,
    history:WorldlineHistory,
    observer_time_s:f64,
    sim_speed:f64,
    scale_px_rg:f32,
    probe_alpha_rg:f64,
    autolock:bool,
    acquisition:Option<OpticalAcquisition>,
    probe:Option<EquatorialNullPath>,
    refresh_s:f64,
}
impl Lab {
    fn new(spin:f64)->Self {
        let geometry=KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap();
        let camera=KerrCamera::new(geometry,85.0,FRAC_PI_2,36.0).unwrap();
        let combatant=source_aligned_to_past_ray(
            &camera,0.0,20.0,0.06,1,1.0).unwrap();
        let history=combatant.record(
            &camera,-0.01,0.15,0.00008).unwrap();
        let mut lab=Self{
            camera,combatant,history,observer_time_s:0.06,
            sim_speed:1.0,scale_px_rg:5.9,probe_alpha_rg:0.0,
            autolock:true,acquisition:None,probe:None,
            refresh_s:10.0,
        };
        lab.refresh();
        lab
    }
    fn refresh(&mut self) {
        self.acquisition=if self.autolock {
            observe_combatant(&self.camera,&self.history,
                self.observer_time_s,self.combatant.envelope_rg,
                6.0,97,17.0)
        } else {None};
        let alpha=self.acquisition.as_ref()
            .map(|p|p.alpha_rg)
            .unwrap_or(self.probe_alpha_rg);
        self.probe=self.camera.equatorial_null_path(
            alpha,17.0,3600);
    }
    fn controls(&mut self,dt:f64) {
        let new_spin=if is_key_pressed(KeyCode::Key1) {Some(0.0)}
            else if is_key_pressed(KeyCode::Key2) {Some(0.75)}
            else if is_key_pressed(KeyCode::Key3) {Some(-0.75)}
            else {None};
        if let Some(s)=new_spin {*self=Self::new(s);}
        if is_key_pressed(KeyCode::R) {
            let spin=self.camera.geometry.spin;*self=Self::new(spin);
        }
        if is_key_pressed(KeyCode::F) {self.autolock=!self.autolock;self.refresh_s=10.0;}
        if is_key_pressed(KeyCode::Q) {self.sim_speed=(self.sim_speed*0.5).max(0.125);}
        if is_key_pressed(KeyCode::E) {self.sim_speed=(self.sim_speed*2.0).min(16.0);}
        if is_key_down(KeyCode::A) {
            self.probe_alpha_rg=(self.probe_alpha_rg-4.0*dt).max(-22.0);
            self.autolock=false;self.refresh_s=10.0;
        }
        if is_key_down(KeyCode::D) {
            self.probe_alpha_rg=(self.probe_alpha_rg+4.0*dt).min(22.0);
            self.autolock=false;self.refresh_s=10.0;
        }
        if is_key_down(KeyCode::Z) {self.scale_px_rg=(self.scale_px_rg/1.03).max(1.5);}
        if is_key_down(KeyCode::X) {self.scale_px_rg=(self.scale_px_rg*1.03).min(12.0);}
        let (_,wheel)=mouse_wheel();
        if wheel!=0.0 {
            self.scale_px_rg=(self.scale_px_rg*(1.0+0.08*wheel))
                .clamp(1.5,12.0);
        }
        self.observer_time_s+=dt*0.002*self.sim_speed;
        if self.observer_time_s>=0.145 {self.observer_time_s=0.02;}
        self.refresh_s+=dt;
        if self.refresh_s>=0.16 {
            self.refresh();self.refresh_s=0.0;
        }
    }
    fn pixel(&self,p:V2)->Vec2 {
        let rg=self.camera.geometry.r_g_m;
        let center=vec2(screen_width()*0.44,screen_height()*0.54);
        center+vec2((p.x/rg) as f32,(p.y/rg) as f32)*
            self.scale_px_rg
    }
    fn world_ring(&self,r:f64,color:Color,label:&str) {
        let center=self.pixel(V2::default());
        let rpx=r as f32*self.scale_px_rg;
        draw_circle_lines(center.x,center.y,rpx,1.4,
            Color::new(color.r,color.g,color.b,0.6));
        if rpx<125.0 {
            draw_text(label,center.x+rpx+5.0,
                center.y+12.0,16.0,color);
        }
    }
    fn draw_ray(&self,path:&EquatorialNullPath) {
        for pair in path.samples.windows(2) {
            let a=self.pixel(pair[0].chart_xy_m);
            let b=self.pixel(pair[1].chart_xy_m);
            draw_line(a.x,a.y,b.x,b.y,2.3,SKYBLUE);
        }
    }
    fn hud(&self) {
        draw_rectangle(14.0,14.0,920.0,195.0,
            Color::new(0.015,0.027,0.054,0.94));
        draw_text("BLIGHTED GALAXY / KERR CAUSAL OPTICAL LINK",
            27.0,48.0,29.0,WHITE);
        draw_text("1 Schwarzschild   2 Kerr +0.75   3 Kerr -0.75",
            27.0,82.0,22.0,SKYBLUE);
        draw_text("F auto-acquire/probe   A/D steer probe   Q/E time rate",
            27.0,114.0,19.0,LIGHTGRAY);
        draw_text("R restart   Z/X or wheel zoom   Escape exit",
            27.0,143.0,19.0,LIGHTGRAY);
        draw_text(&format!("SPIN a={:+.2}  OBSERVER t={:.4}s  PROBE alpha={:+.3} rg",
            self.camera.geometry.spin,self.observer_time_s,
            self.acquisition.as_ref().map(|p|p.alpha_rg)
                .unwrap_or(self.probe_alpha_rg)),
            27.0,174.0,19.0,GOLD);
        draw_text(if self.autolock {"AUTO SENSOR LOCK"}
            else {"MANUAL PHOTON PROBE"},
            27.0,200.0,17.0,if self.autolock {GREEN}else{YELLOW});

        draw_rectangle(14.0,screen_height()-223.0,
            650.0,190.0,Color::new(0.015,0.027,0.054,0.94));
        let rg=self.camera.geometry.r_g_m;
        let rows=if let Some(acq)=&self.acquisition {
            vec![
                format!("ARRIVAL TIME       {:9.6} s",acq.observer_time_s),
                format!("EMISSION TIME      {:9.6} s",acq.emission_time_s),
                format!("CURVED LIGHT DELAY {:9.3} ms",acq.light_delay_s*1e3),
                format!("PAST/PRESENT DELTA {:9.4} rg",
                    (acq.observed_pos_m-acq.present_pos_m).len()/rg),
                format!("RADIAL/POLAR ERR   {:.2e} / {:.2e}",
                    acq.radial_error,acq.polar_error),
                format!("RAY CHORD SAMPLES  {:9}",acq.ray.samples.len()),
            ]
        } else {
            vec![
                "CAUSAL SENSOR        NO MATCH".into(),
                "RECORDED WORLDLINE   AVAILABLE".into(),
                "Probe manually or toggle F auto-acquire".into(),
                "Source moves on a circular Kerr test geodesic".into(),
                "Collision envelope defined in BL coordinate chart".into(),
                "Waiting for optical path intersection".into(),
            ]
        };
        for (i,row) in rows.iter().enumerate() {
            draw_text(row,28.0,screen_height()-194.0+i as f32*27.0,
                18.0,if i==2 {GOLD}else{LIGHTGRAY});
        }
        draw_text("Photon path is curved Kerr null geodesic; target is historical orbital worldline.",
            24.0,screen_height()-14.0,17.0,GRAY);
    }
    fn draw(&self) {
        clear_background(Color::new(0.008,0.014,0.040,1.0));
        for i in 0..410_i32 {
            let px=((i as i64*2357+2713)%5003) as f32/5003.0*screen_width();
            let py=((i as i64*4121+977)%5009) as f32/5009.0*screen_height();
            draw_circle(px,py,if i%39==0 {1.8}else{0.8},
                Color::new(0.6,0.7,1.0,0.24));
        }
        let geo=self.camera.geometry;
        self.world_ring(geo.outer_horizon_rg(),ORANGE,"HORIZON");
        self.world_ring(geo.photon_circular_radius_rg(1).unwrap(),GOLD,"+L PHOTON");
        self.world_ring(geo.isco_rg(1).unwrap(),VIOLET,"+L ISCO");
        self.world_ring(self.combatant.orbital_radius_rg,
            Color::new(0.6,0.7,0.95,0.24),"");
        let origin=self.pixel(V2::default());
        let h_radius=geo.outer_horizon_rg() as f32*self.scale_px_rg;
        draw_circle(origin.x,origin.y,h_radius,BLACK);
        draw_circle_lines(origin.x,origin.y,h_radius,2.0,ORANGE);
        if let Some(path)=self.probe.as_ref() {
            self.draw_ray(path);
        }
        let observer_pos=V2::new(self.camera.observer_rg*geo.r_g_m,0.0);
        let observer=self.pixel(observer_pos);
        draw_circle(observer.x,observer.y,10.0,WHITE);
        draw_circle_lines(observer.x,observer.y,18.0,1.0,WHITE);
        draw_text("ZAMO OBSERVER",observer.x-50.0,observer.y-30.0,18.0,WHITE);
        if let Some(source)=self.history.sample_at(self.observer_time_s) {
            let pos=self.pixel(source.pos_m);
            draw_circle(pos.x,pos.y,12.0,RED);
            draw_circle_lines(pos.x,pos.y,17.0,2.0,RED);
            draw_text("PRESENT TARGET (DEBUG)",pos.x+16.0,
                pos.y-11.0,17.0,RED);
        }
        if let Some(acq)=&self.acquisition {
            let ghost=self.pixel(acq.observed_pos_m);
            draw_circle_lines(ghost.x,ghost.y,16.0,2.0,YELLOW);
            draw_circle(ghost.x,ghost.y,5.0,YELLOW);
            draw_text("RECEIVED SOURCE LIGHT",ghost.x+17.0,
                ghost.y+5.0,17.0,GOLD);
            let hit=self.pixel(acq.photon_chart_pos_m);
            draw_circle_lines(hit.x,hit.y,22.0,2.4,SKYBLUE);
        }
        self.hud();
    }
}
#[macroquad::main(settings)]
async fn main() {
    let mut lab=Lab::new(0.75);
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        lab.controls(get_frame_time() as f64);
        lab.draw();
        next_frame().await;
    }
}
