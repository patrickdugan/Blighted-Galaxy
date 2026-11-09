//! Cold hydrostatic white-dwarf profile + 3D degenerate SPH demonstration.
 //! cargo run --release --bin white_dwarf_lab
use blighted_galaxy::physics::SOLAR_MASS_KG;
use blighted_galaxy::sph3::{
    self,Fluid3,Massive3,Sph3Config,V3,
};
use blighted_galaxy::white_dwarf::{
    DegenerateElectrons,WhiteDwarfProfile,profile_for_mass,
    solve_hydrostatic,
};
use blighted_galaxy::particle_batch::{
    ParticleBillboard,build_meshes,draw_billboards,
};
use blighted_galaxy::stellar_relax::{
    relax_isolated_white_dwarf,RelaxParameters,RelaxReport,
};
use macroquad::prelude::*;

fn config()->Conf {
    Conf {window_title:"Blighted Galaxy • Cold White-Dwarf Structure".into(),
        window_width:1600,window_height:950,high_dpi:true,
        ..Default::default()}
}
struct WdLab {
    profile:WhiteDwarfProfile,
    mass_radius:Vec<(f64,f64)>,
    particles:Vec<Fluid3>,
    central:Massive3,
    config:Sph3Config,
    energy_initial:f64,
    momentum_initial:V3,
    mass_initial:f64,
    elapsed_s:f64,
    accumulator_s:f64,
    time_rate:f64,
    paused:bool,
    rotation:f32,
    scale:f32,
    relaxation:Option<RelaxReport>,
}
impl WdLab {
    fn new(target_solar:f64)->Self {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let profile=profile_for_mass(target_solar*SOLAR_MASS_KG,eos,10_000.0)
            .expect("reference white dwarf must converge");
        let mut curve=Vec::new();
        for i in 0..38 {
            let rho=10_f64.powf(7.5+i as f64*4.5/37.0);
            let star=solve_hydrostatic(rho,eos,20_000.0).unwrap();
            curve.push((star.mass_solar(),star.radius_m/1e6));
        }
        let center=V3::new(1.2e8,0.0,0.0);
        let bulk=V3::new(-3.0e5,2.0e4,0.0);
        let particles=sph3::seed_from_hydrostatic(
            &profile,center,bulk,192);
        let central=Massive3 {
            pos_m:V3::default(),vel_m_s:V3::default(),
            mass_kg:10.0*SOLAR_MASS_KG,radius_m:3.0e4,
        };
        let config=Sph3Config::cold_white_dwarf(2.0).unwrap();
        let energy_initial=sph3::energy_j(&central,&particles,config);
        let momentum_initial=sph3::total_momentum(&central,&particles);
        let mass_initial=sph3::total_mass_kg(&central,&particles);
        Self{profile,mass_radius:curve,particles,central,config,
            energy_initial,momentum_initial,mass_initial,
            elapsed_s:0.0,accumulator_s:0.0,time_rate:1.0,
            paused:true,rotation:0.36,scale:2.7,relaxation:None}
    }
    fn controls(&mut self) {
        let target=if is_key_pressed(KeyCode::Key1) {Some(0.35)}
            else if is_key_pressed(KeyCode::Key2) {Some(0.65)}
            else if is_key_pressed(KeyCode::Key3) {Some(0.95)}
            else {None};
        if let Some(m)=target {*self=Self::new(m);}
        if is_key_pressed(KeyCode::R) {
            *self=Self::new(self.profile.mass_solar());
        }
        if is_key_pressed(KeyCode::Space) {self.paused=!self.paused;}
        if is_key_pressed(KeyCode::H) {
            // Reset to a hydrostatic stratified seed before relaxing the
            // star in isolation, then introduce it into the BH encounter.
            self.paused=true;
            self.particles=sph3::seed_from_hydrostatic(
                &self.profile,V3::new(1.2e8,0.0,0.0),
                V3::new(-3.0e5,2.0e4,0.0),192);
            self.central=Massive3{
                pos_m:V3::default(),vel_m_s:V3::default(),
                mass_kg:10.0*SOLAR_MASS_KG,radius_m:3.0e4
            };
            self.relaxation=relax_isolated_white_dwarf(
                &mut self.particles,self.config,&self.profile,
                RelaxParameters {
                    iterations:1500,dt_s:0.01,damping_time_s:0.7,
                }).ok();
            self.energy_initial=sph3::energy_j(
                &self.central,&self.particles,self.config);
            self.momentum_initial=sph3::total_momentum(
                &self.central,&self.particles);
            self.mass_initial=sph3::total_mass_kg(
                &self.central,&self.particles);
            self.elapsed_s=0.0;self.accumulator_s=0.0;
        }
        if is_key_pressed(KeyCode::Q) {self.time_rate=(self.time_rate/2.0).max(0.1);}
        if is_key_pressed(KeyCode::E) {self.time_rate=(self.time_rate*2.0).min(12.0);}
        if is_key_down(KeyCode::A) {self.rotation-=0.015;}
        if is_key_down(KeyCode::D) {self.rotation+=0.015;}
        if is_key_down(KeyCode::Z) {self.scale=(self.scale/1.02).max(0.5);}
        if is_key_down(KeyCode::X) {self.scale=(self.scale*1.02).min(30.0);}
    }
    fn step(&mut self,dt_real:f64) {
        if self.paused {return;}
        self.accumulator_s=(self.accumulator_s+dt_real.min(0.05)*
            self.time_rate).min(0.12);
        let mut count=0;
        let dt=0.002;
        while self.accumulator_s>=dt && count<90 {
            sph3::advance(&mut self.central,&mut self.particles,
                self.config,dt);
            self.elapsed_s+=dt;
            self.accumulator_s-=dt;
            count+=1;
        }
    }
    fn project(&self,p:V3)->Vec2 {
        let relative=(p-self.central.pos_m)/1e6;
        let theta=self.rotation as f64;
        let y=relative.y*theta.cos()-relative.z*theta.sin();
        vec2(screen_width()*0.68+relative.x as f32*self.scale,
            screen_height()*0.51+y as f32*self.scale)
    }
    fn draw_curve(&self) {
        let base=vec2(35.0,screen_height()-280.0);
        let size=vec2(560.0,188.0);
        draw_rectangle(base.x-10.0,base.y-32.0,size.x+35.0,size.y+72.0,
            Color::new(0.01,0.025,0.055,0.9));
        draw_text("COLD DEGENERACY MASS–RADIUS CURVE",
            base.x,base.y-9.0,19.0,SKYBLUE);
        for pair in self.mass_radius.windows(2) {
            let pos=|v:(f64,f64)| {
                vec2(base.x+(v.0/1.5) as f32*size.x,
                    base.y+size.y-(v.1/25.0) as f32*size.y)
            };
            let a=pos(pair[0]);let b=pos(pair[1]);
            draw_line(a.x,a.y,b.x,b.y,2.0,GOLD);
        }
        let marker=vec2(
            base.x+(self.profile.mass_solar()/1.5) as f32*size.x,
            base.y+size.y-(self.profile.radius_m/25.0e6) as f32*size.y);
        draw_circle(marker.x,marker.y,7.0,SKYBLUE);
        draw_text("0 Msun",base.x,base.y+size.y+21.0,16.0,LIGHTGRAY);
        draw_text("1.5 Msun",base.x+size.x-77.0,
            base.y+size.y+21.0,16.0,LIGHTGRAY);
        draw_text("25 Mm",base.x+2.0,base.y+18.0,16.0,LIGHTGRAY);
    }
    fn draw(&self) {
        clear_background(Color::new(0.012,0.016,0.035,1.0));
        for i in 0..380 {
            let x=(i*2741%1423) as f32/1423.0*screen_width();
            let y=(i*3877%1487) as f32/1487.0*screen_height();
            draw_circle(x,y,0.8,Color::new(0.65,0.72,1.0,0.25));
        }
        let bh=self.project(self.central.pos_m);
        draw_circle(bh.x,bh.y,20.0,GOLD);
        draw_circle(bh.x,bh.y,16.0,BLACK);
        draw_circle_lines(bh.x,bh.y,22.0,2.0,ORANGE);
        let density=sph3::densities_kg_m3(
            &self.particles,self.config.smoothing_m);
        let peak=density.iter().copied().fold(0.0_f64,f64::max);
        let sprites=self.particles.iter().zip(&density)
            .map(|(p,rho)| {
                let heat=(*rho/peak).clamp(0.0,1.0) as f32;
                ParticleBillboard{
                    center_px:self.project(p.pos_m),radius_px:3.6,
                    tint:Color::new(0.35+heat*0.63,0.8-0.3*heat,1.0,0.85),
                    depth:0.0
                }
            }).collect::<Vec<_>>();
        draw_billboards(&build_meshes(&sprites));
        let momentum=sph3::total_momentum(&self.central,&self.particles);
        let energy=sph3::energy_j(&self.central,&self.particles,self.config);
        let p_drift=(momentum-self.momentum_initial).len()/
            (self.mass_initial*3.0e5);
        let energy_drift=(energy-self.energy_initial).abs()/
            self.energy_initial.abs().max(1.0);
        draw_rectangle(14.0,14.0,920.0,385.0,
            Color::new(0.012,0.024,0.047,0.93));
        draw_text("BLIGHTED GALAXY / WHITE-DWARF STRUCTURE",
            27.0,47.0,26.0,WHITE);
        draw_text("1: 0.35 Msun   2: 0.65 Msun   3: 0.95 Msun",
            27.0,78.0,19.0,SKYBLUE);
        draw_text("SPACE pause | H 15-second damped relaxation | A/D rotate | Q/E rate | Z/X zoom",
            27.0,105.0,18.0,LIGHTGRAY);
        let rows=[
            format!("Mass            {:9.6} solar masses",self.profile.mass_solar()),
            format!("Radius          {:9.2} km",self.profile.radius_m/1e3),
            format!("Central density {:9.3e} kg/m³",self.profile.rho_c_kg_m3),
            format!("Central pressure{:9.3e} Pa",
                self.profile.eos.pressure_pa(self.profile.rho_c_kg_m3)),
            format!("Elapsed         {:9.4} s       {}",
                self.elapsed_s,if self.paused {"PAUSED"} else {"RUNNING"}),
            format!("SPH ΔE/|E0|     {:9.3e}       ΔP/scale {:9.3e}",
                energy_drift,p_drift),
            if let Some(r)=self.relaxation {
                format!("RELAX FORCE RATIO {:9.5}    {} iterations",
                    r.force_reduction(),r.iterations)
            } else {"RELAX STATUS    H: prepare damped initial state".to_string()},
            if let Some(r)=self.relaxation {
                format!("LOG-RHO MISMATCH {:7.4}->{:.4} | max R {:.1} km",
                    r.initial_density_log_rms,r.final_density_log_rms,
                    r.max_radius_m/1000.0)
            } else {"SOURCE PROFILE   hydrostatic 1D, SPH discrete seed".to_string()},
        ];
        for (i,row) in rows.iter().enumerate() {
            draw_text(row,27.0,141.0+i as f32*29.0,19.0,
                if i==5 {GOLD}else{WHITE});
        }
        self.draw_curve();
        draw_text("Hydrostatic Newtonian reference; 192 equal-mass SPH parcels under a cold electron EOS.",
            25.0,screen_height()-20.0,18.0,LIGHTGRAY);
    }
}
#[macroquad::main(config)]
async fn main() {
    let mut scene=WdLab::new(0.65);
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        scene.controls();
        scene.step(get_frame_time() as f64);
        scene.draw();
        next_frame().await;
    }
}
