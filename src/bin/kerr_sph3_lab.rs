//! Kerr equatorial photon lensing and three-dimensional stellar debris lab.
 //!
 //! cargo run --release --bin kerr_sph3_lab
 //! 1 Schwarzschild light paths, 2 Kerr a*=0.85 light paths,
 //! 3 a*=0.85 reversed spin, 4 SPH3 debris plus batched GPU billboards.
use blighted_galaxy::kerr::{KerrGeometry,KerrPhoton,PhotonStatus};
use blighted_galaxy::particle_batch::{
    build_meshes,count_indexed_quads,draw_billboards,ParticleBillboard,
};
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::si_bridge::{Seconds,SiTidalEncounter};
use blighted_galaxy::sph3::{
    self,Fluid3,Massive3,Sph3Config,V3,
};
use blighted_galaxy::V2;
use macroquad::prelude::*;

fn window_conf()->Conf {
    Conf {window_title:"Blighted Galaxy: Kerr photons + 3D debris".to_string(),
        window_width:1600,window_height:950,high_dpi:true,
        ..Default::default()}
}
#[derive(Clone,Copy,PartialEq,Eq)]
enum Scene {Schwarzschild,Spinning,Reversed,Debris3}
impl Scene {
    fn title(self)->&'static str {
        match self {
            Self::Schwarzschild=>"1  NONROTATING BH / EQUATORIAL PHOTONS",
            Self::Spinning=>"2  KERR a*=+0.85 / FRAME DRAGGING",
            Self::Reversed=>"3  KERR a*=-0.85 / REVERSED FRAME DRAGGING",
            Self::Debris3=>"4  WHITE-DWARF TIDAL DEBRIS / 3D SPH POLYTROPE",
        }
    }
    fn spin(self)->f64 {
        match self {
            Self::Schwarzschild|Self::Debris3=>0.0,
            Self::Spinning=>0.85,
            Self::Reversed=>-0.85,
        }
    }
}
struct Trace {
    points_rg:Vec<V2>,
    impact:f64,
    status:PhotonStatus,
    radial_error:f64,
}
fn trace_bundle(spin:f64)->(KerrGeometry,Vec<Trace>) {
    let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap();
    let mut traces=Vec::new();
    for k in -11..=11 {
        let b=k as f64*0.85;
        let mut ray=KerrPhoton::inbound(bh,46.0,b).unwrap();
        let mut pts=vec![ray.pos_rg()];
        for step in 0..3200 {
            ray.advance(0.055).unwrap();
            if step%4==0 || ray.status!=PhotonStatus::Active {
                pts.push(ray.pos_rg());
            }
            if ray.status!=PhotonStatus::Active {break;}
        }
        traces.push(Trace{points_rg:pts,impact:b,
            status:ray.status,radial_error:ray.radial_constraint_residual()});
    }
    (bh,traces)
}
struct StarDebris3 {
    progenitor:SiTidalEncounter,
    central:Option<Massive3>,
    particles:Vec<Fluid3>,
    initial_momentum:V3,
    initial_mass:f64,
    initial_energy:f64,
    now_s:f64,
    breakup_s:Option<f64>,
}
impl StarDebris3 {
    fn new()->Self {
        Self {progenitor:SiTidalEncounter::new(),
            central:None,particles:vec![],
            initial_momentum:V3::default(),
            initial_mass:0.0,initial_energy:0.0,now_s:0.0,
            breakup_s:None}
    }
    fn tick(&mut self,dt:f64) {
        if let Some(ref mut bh)=self.central {
            sph3::advance(bh,&mut self.particles,
                Sph3Config::white_dwarf_proxy(),dt);
        } else {
            self.progenitor.step(Seconds(dt));
            if self.progenitor.disruption.is_some() {
                let (particles,_)=sph3::seed_from_progenitor(
                    self.progenitor.white_dwarf,0.0,0.0,128);
                let hole=self.progenitor.black_hole;
                let bh=Massive3 {pos_m:V3::new(hole.pos_m.x,hole.pos_m.y,0.0),
                    vel_m_s:V3::new(hole.vel_m_s.x,hole.vel_m_s.y,0.0),
                    mass_kg:hole.mass_kg,radius_m:hole.radius_m};
                self.initial_momentum=sph3::total_momentum(&bh,&particles);
                self.initial_mass=sph3::total_mass_kg(&bh,&particles);
                self.initial_energy=sph3::energy_j(
                    &bh,&particles,Sph3Config::white_dwarf_proxy());
                self.central=Some(bh);
                self.particles=particles;
                self.breakup_s=Some(self.progenitor.time_s);
            }
        }
        self.now_s+=dt;
    }
    fn bh_position(&self)->V3 {
        self.central.as_ref().map(|b|b.pos_m).unwrap_or_else(||{
            let p=self.progenitor.black_hole.pos_m;
            V3::new(p.x,p.y,0.0)
        })
    }
}
struct Game {
    scene:Scene,
    geometry:KerrGeometry,
    traces:Vec<Trace>,
    debris:StarDebris3,
    gpu_stars:Vec<macroquad::models::Mesh>,
    scale:f32,
    angle:f32,
    rate:f64,
    paused:bool,
    accumulator:f64,
}
impl Game {
    fn new(scene:Scene)->Self {
        let (geometry,traces)=trace_bundle(scene.spin());
        let mut sprites=Vec::with_capacity(9_500);
        // The 9,500 deterministic background particles are a GPU batch /
        // draw-call stress specimen, not an astrophysical simulation.
        let mut seed=0xBE_11_47_59_u64;
        for _ in 0..9_500 {
            seed=seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let x=((seed>>32) as f64/u32::MAX as f64) as f32*1600.0;
            seed=seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let y=((seed>>32) as f64/u32::MAX as f64) as f32*950.0;
            sprites.push(ParticleBillboard{center_px:vec2(x,y),
                radius_px:if seed%11==0 {1.4} else {0.65},
                tint:Color::new(0.50,0.66,0.85,0.28),
                depth:0.0});
        }
        let gpu_stars=build_meshes(&sprites);
        Self {scene,geometry,traces,debris:StarDebris3::new(),
            gpu_stars,scale:if scene==Scene::Debris3 {3.0}else{10.0},
            angle:0.30,rate:15.0,paused:false,accumulator:0.0}
    }
    fn controls(&mut self) {
        let next=if is_key_pressed(KeyCode::Key1) {Some(Scene::Schwarzschild)}
            else if is_key_pressed(KeyCode::Key2) {Some(Scene::Spinning)}
            else if is_key_pressed(KeyCode::Key3) {Some(Scene::Reversed)}
            else if is_key_pressed(KeyCode::Key4) {Some(Scene::Debris3)}
            else {None};
        if let Some(scene)=next {*self=Self::new(scene);}
        if is_key_pressed(KeyCode::R) {*self=Self::new(self.scene);}
        if is_key_pressed(KeyCode::Space) {self.paused=!self.paused;}
        if is_key_pressed(KeyCode::Q) {self.rate=(self.rate*0.5).max(0.2);}
        if is_key_pressed(KeyCode::E) {self.rate=(self.rate*2.0).min(64.0);}
        let (_,scroll)=mouse_wheel();
        if scroll!=0.0 {self.scale*=(1.0+scroll*0.07).clamp(0.4,2.0);}
        if is_key_down(KeyCode::Z) {self.scale*=0.98;}
        if is_key_down(KeyCode::X) {self.scale*=1.02;}
        if is_key_down(KeyCode::A) {self.angle-=0.014;}
        if is_key_down(KeyCode::D) {self.angle+=0.014;}
        self.scale=self.scale.clamp(0.4,100.0);
    }
    fn update(&mut self,dt_real:f64) {
        if self.paused||self.scene!=Scene::Debris3 {return;}
        self.accumulator=(self.accumulator+dt_real.min(0.05)*self.rate).min(3.0);
        let mut loops=0;
        while loops<140 {
            let dt=if self.debris.central.is_some() {0.004} else {0.1};
            if self.accumulator<dt {break;}
            self.debris.tick(dt);
            self.accumulator-=dt;
            loops+=1;
        }
    }
    fn center()->Vec2 {
        vec2(screen_width()*0.55,screen_height()*0.51)
    }
    fn pixel_kerr(&self,position_rg:V2)->Vec2 {
        Self::center()+vec2(position_rg.x as f32,position_rg.y as f32)*self.scale
    }
    fn pixel_debris(&self,position_m:V3)->Vec2 {
        let p=(position_m-self.debris.bh_position())/1e6;
        let c=self.angle.cos() as f64;
        let s=self.angle.sin() as f64;
        Self::center()+vec2(p.x as f32,(p.y*c-p.z*s) as f32)*self.scale
    }
    fn hud(&self) {
        draw_rectangle(12.0,12.0,990.0,183.0,
            Color::new(0.016,0.026,0.05,0.92));
        draw_text("BLIGHTED GALAXY • KERR PHOTONS / 3D STELLAR FLUID",
            24.0,47.0,27.0,WHITE);
        draw_text(self.scene.title(),24.0,80.0,21.0,SKYBLUE);
        draw_text("1–4 SCENES   |   Q/E TIME RATE   |   Z/X OR WHEEL ZOOM",
            24.0,112.0,19.0,LIGHTGRAY);
        draw_text("SPACE PAUSE  |  R RESET  |  A/D ROTATE FLUID  |  ESC EXIT",
            24.0,144.0,19.0,LIGHTGRAY);
        draw_text(&format!("RUNNING: {}    DRAW CALLS: {} FOR 9,500 BACKGROUND QUADS",
            if self.paused {"PAUSED"} else {"ACTIVE"},self.gpu_stars.len()),
            24.0,177.0,17.0,GOLD);
    }
    fn metrics(&self,lines:&[(&str,String)]) {
        let y=screen_height()-265.0;
        draw_rectangle(14.0,y-33.0,640.0,254.0,
            Color::new(0.012,0.025,0.048,0.90));
        draw_text("NUMERICAL DIAGNOSTICS",25.0,y-7.0,19.0,SKYBLUE);
        for (i,(name,value)) in lines.iter().enumerate() {
            draw_text(name,25.0,y+24.0+i as f32*28.0,18.0,LIGHTGRAY);
            draw_text(value,275.0,y+24.0+i as f32*28.0,18.0,WHITE);
        }
    }
    fn paint_kerr(&self) {
        let bh=self.geometry;
        let center=Self::center();
        for (radius,color,label) in [
            (bh.outer_horizon_rg(),ORANGE,"HORIZON"),
            (bh.ergosphere_equatorial_rg(),VIOLET,"EQUATORIAL ERGOSPHERE"),
            (bh.isco_rg(1).unwrap(),SKYBLUE,"POSITIVE-L ISCO"),
            (bh.isco_rg(-1).unwrap(),DARKBLUE,"NEGATIVE-L ISCO"),
        ] {
            let r=radius as f32*self.scale;
            if r>0.0 {
                draw_circle_lines(center.x,center.y,r,1.0,
                    Color::new(color.r,color.g,color.b,0.6));
                draw_text(label,center.x+r+5.0,center.y+13.0,14.0,color);
            }
        }
        draw_circle(center.x,center.y,
            (bh.outer_horizon_rg() as f32*self.scale).max(2.0),BLACK);
        draw_circle_lines(center.x,center.y,
            (bh.outer_horizon_rg() as f32*self.scale).max(2.0),2.0,ORANGE);
        let mut captured=0;
        let mut escaped=0;
        let mut residual=0.0_f64;
        for trace in &self.traces {
            let color=if trace.impact>0.0 {SKYBLUE}
                else if trace.impact<0.0 {VIOLET}else{YELLOW};
            for pair in trace.points_rg.windows(2) {
                let a=self.pixel_kerr(pair[0]);let b=self.pixel_kerr(pair[1]);
                draw_line(a.x,a.y,b.x,b.y,1.25,
                    Color::new(color.r,color.g,color.b,0.70));
            }
            match trace.status {
                PhotonStatus::Escaped=>escaped+=1,
                PhotonStatus::HorizonChartCutoff=>captured+=1,
                PhotonStatus::Active=>{},
            }
            residual=residual.max(trace.radial_error.abs());
        }
        let lines=[
            ("SPIN a*",format!("{:.3}",bh.spin)),
            ("OUTER HORIZON",format!("{:.4} GM/c²",bh.outer_horizon_rg())),
            ("ISCO +L / -L",format!("{:.3} / {:.3} GM/c²",
                bh.isco_rg(1).unwrap(),bh.isco_rg(-1).unwrap())),
            ("CRITICAL b +L",format!("{:.5} r_g",bh.photon_critical_b_rg(1).unwrap())),
            ("CRITICAL b -L",format!("{:.5} r_g",bh.photon_critical_b_rg(-1).unwrap())),
            ("FRAME DRAG at 5rg",format!("{:.4} rad/s",
                bh.zamo_omega_rad_s(5.0).unwrap())),
            ("CAPTURED / ESCAPED",format!("{} / {}",captured,escaped)),
        ];
        self.metrics(&lines);
        draw_text(&format!("MAX NULL RADIAL RESIDUAL {:.2e}",residual),
            25.0,screen_height()-22.0,18.0,GRAY);
    }
    fn paint_debris(&self) {
        let scene=&self.debris;
        let c=self.pixel_debris(scene.bh_position());
        draw_circle(c.x,c.y,19.0,ORANGE);
        draw_circle(c.x,c.y,14.0,BLACK);
        draw_circle_lines(c.x,c.y,21.0,2.0,GOLD);
        let mut particles=Vec::new();
        if let Some(bh)=&scene.central {
            let densities=sph3::densities_kg_m3(&scene.particles,
                Sph3Config::white_dwarf_proxy().smoothing_m);
            let max_density=densities.iter().copied().fold(0.0_f64,f64::max);
            for (p,rho) in scene.particles.iter().zip(&densities) {
                let shade=(*rho/max_density).clamp(0.0,1.0) as f32;
                particles.push(ParticleBillboard {
                    center_px:self.pixel_debris(p.pos_m),
                    radius_px:3.2,
                    tint:Color::new(0.3+shade*0.7,0.8-shade*0.45,1.0,0.8),
                    depth:0.0,
                });
            }
            let momentum=sph3::total_momentum(bh,&scene.particles);
            let mass=sph3::total_mass_kg(bh,&scene.particles);
            let e=sph3::energy_j(bh,&scene.particles,
                Sph3Config::white_dwarf_proxy());
            let p_scale=scene.initial_mass*1.0e6;
            let lines=[
                ("BREAKUP",format!("{:.3} s",scene.breakup_s.unwrap_or_default())),
                ("FLUID TIME",format!("{:.3} s",
                    scene.now_s-scene.breakup_s.unwrap_or_default())),
                ("PARCEL COUNT",format!("{}",scene.particles.len())),
                ("TRACED MASS",format!("{:.6} M_sun",mass/SOLAR_MASS_KG)),
                ("MOMENTUM DRIFT",format!("{:.3e}",
                    (momentum-scene.initial_momentum).len()/p_scale)),
                ("ENERGY DRIFT",format!("{:.3e}",
                    (e-scene.initial_energy).abs()/scene.initial_energy.abs())),
                ("POLYTROPIC gamma","5/3, fixed entropy K".into()),
            ];
            self.metrics(&lines);
        } else {
            let star=scene.progenitor.white_dwarf;
            let sprite=self.pixel_debris(V3::new(star.pos_m.x,star.pos_m.y,0.0));
            draw_circle(sprite.x,sprite.y,14.0,SKYBLUE);
            let lines=[
                ("BINARY TIME",format!("{:.2} s",scene.now_s)),
                ("DISTANCE",format!("{:.2} Mm",
                    (scene.progenitor.white_dwarf.pos_m-
                     scene.progenitor.black_hole.pos_m).len()/1e6)),
                ("TIDAL / SELF G",format!("{:.4}",scene.progenitor.peak_tidal_ratio)),
                ("DEBRIS","awaiting q >= 1".into()),
                ("SAMPLE COUNT","128 after break".into()),
                ("POLYTROPIC EOS","gamma=5/3 proxy".into()),
                ("PIPELINE","CPU SPH -> GPU mesh".into()),
            ];
            self.metrics(&lines);
        }
        let meshes=build_meshes(&particles);
        draw_billboards(&meshes);
        draw_text(&format!("DYNAMIC PARTICLES {} / GPU DRAW CALLS {}",
            count_indexed_quads(&meshes),meshes.len()),
            25.0,screen_height()-22.0,19.0,GOLD);
    }
    fn draw(&self) {
        clear_background(Color::new(0.007,0.015,0.035,1.0));
        draw_billboards(&self.gpu_stars);
        if self.scene==Scene::Debris3 {self.paint_debris();}
        else {self.paint_kerr();}
        self.hud();
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut scene=Game::new(Scene::Spinning);
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        scene.controls();
        scene.update(get_frame_time() as f64);
        scene.draw();
        next_frame().await;
    }
}
