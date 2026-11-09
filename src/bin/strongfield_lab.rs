//! Interactive strong-field and 2-D self-gravitating debris laboratory.
 //!
 //! cargo run --release --bin strongfield_lab
 //! 1 stable circular geodesic, 2 precessing eccentric, 3 plunge,
 //! 4 black-hole / white-dwarf SPH-lite tidal debris.
use blighted_galaxy::debris::{self,DebrisConfig};
use blighted_galaxy::physics::{C_SI,G_SI,SOLAR_MASS_KG};
use blighted_galaxy::schwarzschild::{
    gravitational_radius_m,isco_radius_m,OrbitStatus,SchwarzschildOrbit,
};
use blighted_galaxy::si_bridge::{Seconds,SiTidalEncounter};
use blighted_galaxy::V2;
use macroquad::prelude::*;
use std::f64::consts::PI;

fn window_conf()->Conf {
    Conf {
        window_title:"Blighted Galaxy • Strong-field & self-gravitating debris".into(),
        window_width:1600,window_height:950,high_dpi:true,
        ..Default::default()
    }
}
#[derive(Clone,Copy,PartialEq,Eq)]
enum Mode { Circular,Precessing,Plunge,Debris }
impl Mode {
    fn title(self)->&'static str {
        match self {
            Self::Circular=>"1   STABLE GEODESIC • r = 10 GM/c²",
            Self::Precessing=>"2   SCHWARZSCHILD PERIAPSIS PRECESSION",
            Self::Plunge=>"3   TIMELIKE RADIAL PLUNGE • HORIZON CHART CUTOFF",
            Self::Debris=>"4   WHITE-DWARF BREAKUP • SELF-GRAVITATING 2D SPH",
        }
    }
    fn is_geodesic(self)->bool {self!=Self::Debris}
}
struct Lab {
    mode:Mode,
    orbit:Option<SchwarzschildOrbit>,
    tidal:Option<SiTidalEncounter>,
    trail:Vec<V2>,
    speed:f64,
    scale:f32,
    paused:bool,
    accumulated:f64,
    original_momentum:V2,
}
impl Lab {
    fn new(mode:Mode)->Self {
        let mass=10.0*SOLAR_MASS_KG;
        let rg=gravitational_radius_m(mass).unwrap();
        let orbit=match mode {
            Mode::Circular=>SchwarzschildOrbit::circular(mass,10.0*rg),
            Mode::Precessing=>SchwarzschildOrbit::bound_from_p_e(mass,30.0,0.2),
            Mode::Plunge=>SchwarzschildOrbit::from_initial(
                mass,12.0*rg,-0.22*C_SI,0.0),
            Mode::Debris=>None,
        };
        let tidal=if mode==Mode::Debris {
            Some(SiTidalEncounter::new())
        } else {None};
        let original_momentum=tidal.as_ref().map(|t|
            t.black_hole.vel_m_s*t.black_hole.mass_kg+
            t.white_dwarf.vel_m_s*t.white_dwarf.mass_kg).unwrap_or_default();
        Self{
            mode,
            trail:orbit.as_ref().map(|o|vec![o.xy_m()]).unwrap_or_default(),
            orbit,tidal,original_momentum,
            speed:if mode==Mode::Debris {8.0} else {30.0},
            scale:match mode {
                Mode::Circular=>32.0,Mode::Precessing=>11.0,
                Mode::Plunge=>35.0,Mode::Debris=>2.5,
            },
            paused:false,accumulated:0.0,
        }
    }
    fn world_to_screen(&self,p:V2)->Vec2 {
        let normalized=if self.mode.is_geodesic() {
            let r=self.orbit.as_ref().unwrap().r_g_m;
            p/r
        } else {p/1.0e6};
        vec2(
            screen_width()*0.53+normalized.x as f32*self.scale,
            screen_height()*0.50+normalized.y as f32*self.scale,
        )
    }
    fn controls(&mut self) {
        let new_mode=if is_key_pressed(KeyCode::Key1) {Some(Mode::Circular)}
            else if is_key_pressed(KeyCode::Key2) {Some(Mode::Precessing)}
            else if is_key_pressed(KeyCode::Key3) {Some(Mode::Plunge)}
            else if is_key_pressed(KeyCode::Key4) {Some(Mode::Debris)}
            else {None};
        if let Some(m)=new_mode {*self=Lab::new(m);}
        if is_key_pressed(KeyCode::R) {*self=Lab::new(self.mode);}
        if is_key_pressed(KeyCode::Space) {self.paused=!self.paused;}
        if is_key_pressed(KeyCode::Q) {self.speed=(self.speed*0.5).max(0.05);}
        if is_key_pressed(KeyCode::E) {self.speed=(self.speed*2.0).min(300.0);}
        let (_,wheel)=mouse_wheel();
        if wheel!=0.0 {self.scale*=(1.0+wheel*0.075).clamp(0.2,4.0);}
        if is_key_down(KeyCode::Z) {self.scale*=0.98;}
        if is_key_down(KeyCode::X) {self.scale*=1.02;}
        self.scale=self.scale.clamp(0.6,160.0);
    }
    fn step(&mut self,dt_real:f64) {
        if self.paused {return;}
        self.accumulated=(self.accumulated+dt_real.min(0.05)*self.speed).min(5.0);
        let mut steps=0;
        let increment=if self.mode.is_geodesic(){0.035}else{0.10};
        while self.accumulated>=increment && steps<500 {
            if let Some(o)=&mut self.orbit {
                if o.status!=OrbitStatus::Active {self.paused=true;break;}
                let dt_s=increment*o.r_g_m/C_SI;
                o.advance_proper(dt_s).expect("valid orbit integration");
                self.trail.push(o.xy_m());
                if self.trail.len()>18_000 {self.trail.drain(..100);}
            }
            if let Some(t)=&mut self.tidal {t.step(Seconds(increment));}
            self.accumulated-=increment;
            steps+=1;
        }
    }
    fn draw_background(&self) {
        clear_background(Color::new(0.012,0.017,0.044,1.0));
        for i in 0..400_i32 {
            let x=((i as i64*2731+3721)%7309) as f32/7309.0*screen_width();
            let y=((i as i64*4523+1207)%7297) as f32/7297.0*screen_height();
            draw_circle(x,y,if i%19==0 {1.6}else{0.7},
                Color::new(0.75,0.8,1.0,0.13+(i%6) as f32*0.08));
        }
    }
    fn disk(&self,p:Vec2,r:f32,color:Color) {
        for i in (1..=6).rev() {
            draw_circle(p.x,p.y,r*(1.0+i as f32*0.33),
                Color::new(color.r,color.g,color.b,0.03));
        }
        draw_circle(p.x,p.y,r,color);
    }
    fn schwarzschild(&self) {
        let o=self.orbit.as_ref().unwrap();
        let origin=self.world_to_screen(V2::default());
        let ring=|r:f32,color:Color|{
            draw_circle_lines(origin.x,origin.y,r*self.scale,1.2,color);
        };
        ring(6.0,Color::new(0.12,0.59,1.0,0.25));
        ring(3.0,Color::new(1.0,0.8,0.3,0.55));
        ring(2.0,Color::new(1.0,0.32,0.21,0.7));
        draw_circle(origin.x,origin.y,self.scale*2.0,BLACK);
        draw_circle_lines(origin.x,origin.y,self.scale*2.0,2.0,ORANGE);
        draw_text("ISCO 6 r_g",origin.x+6.0*self.scale+8.0,
            origin.y+12.0,17.0,SKYBLUE);
        draw_text("PHOTON 3 r_g",origin.x+3.0*self.scale+6.0,
            origin.y-20.0,17.0,GOLD);
        draw_text("HORIZON 2 r_g",origin.x+2.0*self.scale+7.0,
            origin.y+40.0,17.0,ORANGE);

        for points in self.trail.windows(2) {
            let a=self.world_to_screen(points[0]);
            let b=self.world_to_screen(points[1]);
            draw_line(a.x,a.y,b.x,b.y,2.0,Color::new(0.38,0.86,1.0,0.8));
        }
        let particle=self.world_to_screen(o.xy_m());
        self.disk(particle,6.0,SKYBLUE);
        draw_text("TIMELIKE TEST PARTICLE",particle.x+12.0,particle.y-12.0,17.0,SKYBLUE);

        let lines=[
            ("AREAL RADIUS",format!("{:.7} r_g",o.x)),
            ("PHI UNWRAPPED",format!("{:.6} radians",o.phi_rad)),
            ("SPECIFIC ENERGY",format!("{:.12}",o.epsilon)),
            ("SPECIFIC L",format!("{:.7} (c r_g)",o.ell)),
            ("PROPER TIME",format!("{:.8} seconds",o.proper_time_s)),
            ("COORDINATE TIME",format!("{:.8} seconds",o.coordinate_time_s)),
            ("ENERGY CONSTRAINT",format!("{:.3e}",o.current_energy_constraint_error())),
            ("LAST STATE",format!("{:?}",o.status)),
        ];
        self.metrics(&lines);
    }
    fn tidal(&self) {
        let scene=self.tidal.as_ref().unwrap();
        let black=self.world_to_screen(scene.black_hole.pos_m);
        self.disk(black,17.0,ORANGE);
        draw_circle(black.x,black.y,15.0,BLACK);
        draw_circle_lines(black.x,black.y,17.5,2.0,GOLD);
        draw_text("BLACK HOLE",black.x+25.0,black.y-15.0,18.0,GOLD);
        if scene.disruption.is_none() {
            let star=self.world_to_screen(scene.white_dwarf.pos_m);
            self.disk(star,12.0,SKYBLUE);
            draw_text("WHITE DWARF",star.x+21.0,star.y-14.0,18.0,SKYBLUE);
        } else {
            let config=DebrisConfig::white_dwarf_proxy();
            let densities=debris::surface_densities_kg_m2(
                &scene.fragments,config.h_m);
            let peak=densities.iter().copied().fold(0.0_f64,f64::max);
            for (p,rho) in scene.fragments.iter().zip(&densities) {
                let pixel=self.world_to_screen(p.pos_m);
                if pixel.x<0.0||pixel.x>screen_width()||
                    pixel.y<0.0||pixel.y>screen_height() {continue;}
                let heat=(*rho/peak).clamp(0.0,1.0) as f32;
                draw_circle(pixel.x,pixel.y,3.2,
                    Color::new(0.28+heat*0.65,0.65-heat*0.4,1.0-heat*0.1,0.9));
            }
        }
        let total=scene.total_visible_mass_kg().0;
        let momentum=if scene.disruption.is_some() {
            debris::system_momentum_kg_m_s(&scene.black_hole,&scene.fragments)
        } else {
            scene.black_hole.vel_m_s*scene.black_hole.mass_kg+
                scene.white_dwarf.vel_m_s*scene.white_dwarf.mass_kg
        };
        let drift=(momentum-self.original_momentum).len()/
            (0.65*SOLAR_MASS_KG*6e5);
        let distance=(scene.white_dwarf.pos_m-scene.black_hole.pos_m).len();
        let lines=[
            ("TIME",format!("{:.2} seconds",scene.time_s)),
            ("SEPARATION",format!("{:.3} Mm",distance/1e6)),
            ("PEAK TIDAL RATIO",format!("{:.5}",scene.peak_tidal_ratio)),
            ("FRAGMENTS",format!("{}",scene.fragments.len())),
            ("MASS TRACKED",format!("{:.6} M_sun",total/SOLAR_MASS_KG)),
            ("RELATIVE P DRIFT",format!("{:.3e}",drift)),
            ("DYNAMICS",if scene.disruption.is_some() {
                "2D SELF-GRAVITY + PRESSURE".to_string()
            } else {"BINARY VERLET".to_string()}),
            ("EOS", "ISOTHERMAL 2-D PROXY".to_string()),
        ];
        self.metrics(&lines);
    }
    fn metrics(&self,lines:&[(&str,String)]) {
        let y=screen_height()-288.0;
        draw_rectangle(16.0,y-29.0,610.0,280.0,
            Color::new(0.01,0.023,0.047,0.91));
        draw_text("PHYSICAL INVARIANTS",28.0,y-4.0,20.0,WHITE);
        for (i,(name,value)) in lines.iter().enumerate() {
            draw_text(name,28.0,y+28.0+i as f32*28.0,18.0,LIGHTGRAY);
            draw_text(value,243.0,y+28.0+i as f32*28.0,18.0,
                if i==6 {GOLD}else{SKYBLUE});
        }
    }
    fn draw(&self) {
        self.draw_background();
        if self.mode.is_geodesic() {self.schwarzschild();}
        else {self.tidal();}
        draw_rectangle(16.0,14.0,1020.0,162.0,
            Color::new(0.01,0.025,0.048,0.90));
        draw_text("BLIGHTED GALAXY • STRONG-FIELD / DEBRIS LAB",
            28.0,46.0,29.0,WHITE);
        draw_text(self.mode.title(),28.0,77.0,21.0,SKYBLUE);
        draw_text("1–4 scenarios  |  SPACE pause  |  Q/E time scale  |  Z/X or wheel zoom",
            28.0,108.0,18.0,LIGHTGRAY);
        draw_text(&format!("R reset  |  ESC exit  |  rate {:.2} {} / real second  |  {}",
            self.speed,
            if self.mode.is_geodesic() {"gravitational time units"} else {"SI seconds"},
            if self.paused {"PAUSED"} else {"RUNNING"}),
            28.0,135.0,17.0,GOLD);
        draw_text("SCHWARZSCHILD TEST PARTICLE / ISOTHERMAL 2D FLUID PROXY",
            28.0,161.0,16.0,GRAY);
        draw_text("GR: EQUATORIAL, NONSPINNING, TEST BODY. DEBRIS: SYMMETRIC 2D SPH FORCES.",
            24.0,screen_height()-15.0,17.0,GRAY);
    }
}
#[macroquad::main(window_conf)]
async fn main() {
    let mut lab=Lab::new(Mode::Precessing);
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        lab.controls();
        lab.step(get_frame_time() as f64);
        lab.draw();
        next_frame().await;
    }
}
