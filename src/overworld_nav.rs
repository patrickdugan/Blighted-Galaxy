//! Macroquad galactic navigation, built entirely from procedural pixels.
//!
//! Visual provenance: the user's original NASA-style face-on Milky Way
//! reference with fine radial/longitude stenciling, sparse blood-red
//! star speckles, and thick infrared-darkened Jenkins Swarm patches.
//! Neither the screenshot nor the user's artwork is rasterized into
//! the game. Both star field and moving blight texture regenerate from
//! fixed seeds and physically bounded world positions.
//!
//! Controls: mouse selects destination; Enter launches at 0.76c;
//! Space pause; Q/E calendar speed; ,/. jump +/- 1000 yr;
//! WASD/arrows pan, wheel/Z/X zoom, right drag pan, F focus wavefront,
//! C galactic center, O Orion; G galactic grid, L local sector stencils,
//! B toggle Blight; 1/2 encounters; Tab/3 return to overworld.
use crate::overworld::{
    self,Front,GalacticNavModel,MAP_EXTENT_LY,ORION_OUTBREAK,
    SectorStatus,SOL,WAYPOINTS,
};
use crate::V2;
use crate::strategy::{StrategicState,Project,Faction,TacticalKind,TacticalOutcome};
use macroquad::prelude::*;
use std::f64::consts::{PI,TAU};

const TEX_SIZE:u16=1152;
const GALAXY_DIAMETER_LY:f64=2.0*MAP_EXTENT_LY;
const PANEL_WIDTH:f32=338.0;
const BG:Color=Color::new(0.005,0.007,0.017,1.0);
const BLUE:Color=Color::new(0.37,0.65,0.88,1.0);
const CYAN:Color=Color::new(0.43,0.86,0.94,1.0);
const GOLD:Color=Color::new(0.95,0.79,0.48,1.0);
const REDDIM:Color=Color::new(0.78,0.17,0.17,1.0);
const BONE:Color=Color::new(0.83,0.84,0.80,1.0);

#[derive(Clone,Copy)]
struct Random(u64);
impl Random {
    fn new(seed:u64)->Self{Self(seed)}
    fn next(&mut self)->u64 {
        self.0=self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut v=self.0;
        v=(v^(v>>30)).wrapping_mul(0xbf58476d1ce4e5b9);
        v=(v^(v>>27)).wrapping_mul(0x94d049bb133111eb);
        v^(v>>31)
    }
    fn unit(&mut self)->f64 {
        ((self.next()>>11) as f64)/(1u64<<53) as f64
    }
    fn signed(&mut self)->f64{2.0*self.unit()-1.0}
    fn bell(&mut self)->f64{
        // Irwin-Hall, bounded and stable across platforms.
        (0..6).map(|_|self.unit()).sum::<f64>()-3.0
    }
}
fn tex_xy(p:V2)->(f64,f64){
    let n=TEX_SIZE as f64;
    (n*(p.x+MAP_EXTENT_LY)/GALAXY_DIAMETER_LY,
     n*(MAP_EXTENT_LY-p.y)/GALAXY_DIAMETER_LY)
}
fn set_rgb(image:&mut Image,x:i32,y:i32,rgb:[f32;3]){
    if x<0||y<0||x>=TEX_SIZE as i32||y>=TEX_SIZE as i32{return;}
    let i=(y as usize*TEX_SIZE as usize+x as usize)*4;
    for k in 0..3{
        image.bytes[i+k]=((image.bytes[i+k] as f32/255.0+
            rgb[k]).clamp(0.0,1.0)*255.0) as u8;
    }
}
fn point_glow(image:&mut Image,x:f64,y:f64,rgb:[f32;3],power:f32){
    let ix=x.round() as i32;
    let iy=y.round() as i32;
    set_rgb(image,ix,iy,[rgb[0]*power,rgb[1]*power,rgb[2]*power]);
    let fringe=0.16*power;
    for (dx,dy) in [(-1,0),(1,0),(0,-1),(0,1)]{
        set_rgb(image,ix+dx,iy+dy,[
            rgb[0]*fringe,rgb[1]*fringe,rgb[2]*fringe]);
    }
}
fn gaussian(rng:&mut Random)->f64{rng.bell()/1.15}
fn generate_galaxy()->Image {
    let size=TEX_SIZE as usize;
    let mut image=Image::gen_image_color(TEX_SIZE,TEX_SIZE,
        Color::new(0.003,0.004,0.010,1.0));
    // A diffuse, light-absorbing four-arm disk with a golden central
    // bulge and a tilted bar. No source image is embedded.
    for y in 0..size{
        let world_y=MAP_EXTENT_LY-
            y as f64*GALAXY_DIAMETER_LY/size as f64;
        for x in 0..size{
            let world_x=x as f64*GALAXY_DIAMETER_LY/size as f64-
                MAP_EXTENT_LY;
            let r=(world_x*world_x+world_y*world_y).sqrt();
            let theta=world_y.atan2(world_x);
            let disk=(-r/45_000.0).exp()*
                (1.0/(1.0+((r-52_000.0)/2600.0).exp()));
            let phase=theta-((r+1400.0)/6500.0).ln()/0.29_f64.tan()
                +0.11*(r/4200.0+0.4*theta).sin();
            // Broad, overlapping stellar clouds avoid the thin
            // vector-drawn spiral ribbons of the first prototype.
            let arm=(0.5+0.5*(4.0*phase).cos()).powi(2);
            let dust=(0.49*(world_x/2000.0).sin()+
                0.26*(world_y/1700.0).cos()+
                0.25*((world_x+world_y)/1100.0).sin()).abs();
            let barx=world_x*0.89+world_y*0.46;
            let bary=world_y*0.89-world_x*0.46;
            let bar=(-(barx/12_500.0).powi(2)-
                (bary/2900.0).powi(2)).exp();
            let core=(-(r/5000.0).powi(2)).exp();
            let dust_dark=1.0-0.46*dust*arm;
            let brightness=(disk*(0.70+1.18*arm)*dust_dark)
                .clamp(0.0,1.0) as f32;
            let warm=(0.45*core+0.23*bar).clamp(0.0,1.0) as f32;
            let idx=(y*size+x)*4;
            let channels=[
                0.009+0.18*brightness+0.56*warm,
                0.014+0.25*brightness+0.40*warm,
                0.027+0.34*brightness+0.23*warm,
            ];
            for k in 0..3{
                image.bytes[idx+k]=(255.0_f32*channels[k]
                    .min(1.0)) as u8;
            }
        }
    }
    let mut rng=Random::new(0xB11_9A1A_600D_1993);
    // Four logarithmic stellar arms and an independent diffuse population.
    for i in 0..70_000usize {
        let u=rng.unit().max(1e-9);
        let v=rng.unit().max(1e-9);
        let spiral=i%3!=0;
        let r=if i%13==0 {
            (-u.ln()).sqrt()*6100.0
        }else{
            ((-u.ln()-v.ln())*12_400.0).clamp(1000.0,52_000.0)
        };
        if r>DISK_VISUAL_RADIUS {continue;}
        let arm=(rng.next()%4) as f64;
        let theta=if spiral{
            ((r+1400.0)/6500.0).ln()/0.29_f64.tan()
                +arm*PI/2.0
                -0.11*(r/4200.0+arm*0.58).sin()
                +gaussian(&mut rng)*(0.18+460.0/r)
        }else{rng.unit()*TAU};
        let jitter=gaussian(&mut rng)*(150.0+0.006*r);
        let radius=(r+jitter).max(0.0);
        let pos=V2::from_angle(theta)*radius;
        let (tx,ty)=tex_xy(pos);
        let core=((-(r/9000.0).powi(2)).exp()) as f32;
        let bright=(rng.unit().powi(6)*0.41+0.08) as f32;
        let white=[0.72,0.81,1.0];
        let warm=[1.0,0.79,0.55];
        let tone=[
            white[0]*(1.0-core)+warm[0]*core,
            white[1]*(1.0-core)+warm[1]*core,
            white[2]*(1.0-core)+warm[2]*core,
        ];
        point_glow(&mut image,tx,ty,tone,bright);
        if i%97==0 {
            point_glow(&mut image,tx,ty,[0.74,0.89,1.0],0.45);
        }
    }
    image
}
const DISK_VISUAL_RADIUS:f64=54_000.0;
fn set_overlay(image:&mut Image,x:i32,y:i32,color:[u8;4]){
    if x<0||y<0||x>=TEX_SIZE as i32||y>=TEX_SIZE as i32{return;}
    let off=(y as usize*TEX_SIZE as usize+x as usize)*4;
    image.bytes[off..off+4].copy_from_slice(&color);
}
fn blight_image(front:Front,strategy:&StrategicState)->Image {
    let mut image=Image::gen_image_color(
        TEX_SIZE,TEX_SIZE,Color::new(0.0,0.0,0.0,0.0));
    if front.reach_ly()<5.0{return image;}
    let (cx,cy)=tex_xy(ORION_OUTBREAK);
    let pixel_radius=front.reach_ly()*
        TEX_SIZE as f64/GALAXY_DIAMETER_LY;
    let bound=pixel_radius.ceil() as i32+3;
    for py in (cy as i32-bound).max(0)..
        (cy as i32+bound+1).min(TEX_SIZE as i32) {
        for px in (cx as i32-bound).max(0)..
            (cx as i32+bound+1).min(TEX_SIZE as i32) {
            let w=V2::new(
                px as f64/TEX_SIZE as f64*GALAXY_DIAMETER_LY-
                    MAP_EXTENT_LY,
                MAP_EXTENT_LY-py as f64/TEX_SIZE as f64*
                    GALAXY_DIAMETER_LY);
            let delta=w-ORION_OUTBREAK;
            let r=delta.len();
            if r>=front.reach_ly() {continue;}
            let a=delta.y.atan2(delta.x);
            let local=front.sector_reach_ly(a);
            let grit=0.56*((px as f64*0.23).sin()*
                (py as f64*0.17).cos())+
                0.44*((px as f64*0.71+py as f64*0.47).sin());
            let edge=local-110.0*grit;
            let status=strategy.territory_status(w,front);
            // A captured Blight outpost may push a narrow occupied patch
            // beyond the historical resistance notch, within causality.
            let strategic_blight=strategy.systems.iter().any(|s|{
                s.owner==Faction::Blight &&
                (w-s.pos).len()<=(500.0+s.industry.max(0.0)*70.0+
                    s.defense.max(0.0)*50.0).clamp(300.0,1550.0)
            });
            if r>edge && !strategic_blight {continue;}
            if status==SectorStatus::Unaffected {continue;}
            if status==SectorStatus::Reclaimed{
                set_overlay(&mut image,px,py,[2,23,31,95]);
                continue;
            }
            if status==SectorStatus::Defended{
                set_overlay(&mut image,px,py,[9,10,24,110]);
                continue;
            }
            let dark=(210.0+27.0*grit as f32).clamp(130.0,246.0) as u8;
            let near=(r>local-450.0) as u8;
            let (red,green,blue)=if near>0{(25,3,6)}else{(4,3,7)};
            set_overlay(&mut image,px,py,[red,green,blue,dark]);
        }
    }
    let mut rng=Random::new(0xA11CE_313_1000);
    // Sparse red point sources suggest occupied stars; they are visibly
    // stellar speckles on closer inspection, never solid red paint.
    for _ in 0..1_050 {
        let theta=rng.unit()*TAU;
        let achieved=front.sector_reach_ly(theta);
        let radius=achieved*(rng.unit().sqrt()*0.99);
        let point=ORION_OUTBREAK+V2::from_angle(theta)*radius;
        let status=strategy.territory_status(point,front);
        if status!=SectorStatus::Held&&
            status!=SectorStatus::Contested{continue;}
        let (x,y)=tex_xy(point);
        let brightness=rng.unit();
        let tone=if brightness>0.965{
            [218,32,26,237]
        }else if brightness>0.62{
            [138,15,21,205]
        }else{
            [70,7,15,165]
        };
        set_overlay(&mut image,x.round() as i32,y.round() as i32,tone);
        if brightness>0.992 {
            for (dx,dy) in [(-1,0),(1,0),(0,-1),(0,1)]{
                set_overlay(&mut image,
                    x.round() as i32+dx,y.round() as i32+dy,
                    [90,13,20,160]);
            }
        }
    }
    // Directional Jenkins corridors end BEFORE their local front;
    // the 0.10c causal boundary is a hard spatial invariant.
    for lane in 0..15 {
        let theta=lane as f64*TAU/15.0+0.17;
        let local=front.sector_reach_ly(theta);
        let mut r=320.0;
        while r+180.0<local {
            let bend=0.055*(r/950.0+lane as f64).sin()*
                (r/local).clamp(0.0,1.0);
            let point=ORION_OUTBREAK+
                V2::from_angle(theta+bend)*r;
            if matches!(front.status(point),
                SectorStatus::Held|SectorStatus::Contested){
                let (x,y)=tex_xy(point);
                if rng.unit()>0.75{
                    set_overlay(&mut image,x as i32,y as i32,
                        [114,12,20,189]);
                }
            }
            r+=105.0;
        }
    }
    image
}
fn translucent(c:Color,alpha:f32)->Color{
    Color::new(c.r,c.g,c.b,alpha)
}
fn draw_dashed_circle(p:Vec2,r:f32,tone:Color,strokes:usize){
    if r<=2.0{return;}
    for i in 0..strokes{
        if i%4>=2{continue;}
        let a=i as f32*std::f32::consts::TAU/strokes as f32;
        let b=(i+1) as f32*std::f32::consts::TAU/strokes as f32;
        draw_line(p.x+a.cos()*r,p.y+a.sin()*r,
            p.x+b.cos()*r,p.y+b.sin()*r,
            0.85,tone);
    }
}
fn shadow_text(label:&str,x:f32,y:f32,size:f32,c:Color){
    draw_text(label,x+1.0,y+2.0,size,Color::new(0.0,0.0,0.0,0.94));
    draw_text(label,x,y,size,c);
}
pub struct OverworldNav {
    pub model:GalacticNavModel,
    pub strategy:StrategicState,
    strategy_visible:bool,
    center:V2,
    zoom:f32,
    base:Texture2D,
    blight:Texture2D,
    last_front_bucket:i64,
    last_territory_revision:u64,
    grid:bool,
    local_grid:bool,
    show_blight:bool,
    show_help:bool,
    selected:Option<V2>,
    selected_name:Option<&'static str>,
    dragging:Option<Vec2>,
    encounter_notice:Option<String>,
}
impl Default for OverworldNav{
    fn default()->Self{Self::new()}
}
impl OverworldNav{
    pub fn new()->Self{
        let base=Texture2D::from_image(&generate_galaxy());
        let strategy=StrategicState::default();
        let blight=Texture2D::from_image(
            &blight_image(Front::new(100_000.0),&strategy));
        base.set_filter(FilterMode::Linear);
        blight.set_filter(FilterMode::Linear);
        Self{
            model:GalacticNavModel::default(),
            strategy,strategy_visible:true,
            center:V2::new(-500.0,-8_000.0),
            zoom:1.0,base,blight,
            last_front_bucket:100,last_territory_revision:0,
            grid:true,local_grid:false,show_blight:true,
            show_help:false,selected:Some(SOL),
            selected_name:Some("SOL / ORION SPUR"),
            dragging:None,encounter_notice:None,
        }
    }
    /// Region detail view: the 10-kly ring spans most of the viewport.
    /// Also used by the actual-frame screenshot regression.
    pub fn focus_frontier(&mut self){
        self.center=ORION_OUTBREAK;
        self.zoom=3.7;
        self.local_grid=true;
    }
    fn map_width(&self)->f32{
        (screen_width()-PANEL_WIDTH).max(500.0)
    }
    fn viewport_center(&self)->Vec2{
        vec2(self.map_width()*0.5,screen_height()*0.5+10.0)
    }
    fn pixels_per_ly(&self)->f32 {
        let available=(self.map_width()*0.47)
            .min(screen_height()*0.47);
        (available/54_000.0)*self.zoom
    }
    fn world_to_screen(&self,w:V2)->Vec2{
        let p=self.viewport_center();
        let s=self.pixels_per_ly() as f64;
        vec2(p.x+((w.x-self.center.x)*s) as f32,
            p.y-((w.y-self.center.y)*s) as f32)
    }
    fn screen_to_world(&self,p:Vec2)->V2{
        let c=self.viewport_center();
        let s=self.pixels_per_ly() as f64;
        V2::new(self.center.x+(p.x-c.x) as f64/s,
            self.center.y-(p.y-c.y) as f64/s)
    }
    fn waypoint_at(&self,mouse:Vec2)->Option<usize>{
        WAYPOINTS.iter().enumerate().filter_map(|(i,w)|{
            let q=self.world_to_screen(w.pos);
            let d=(q-mouse).length();
            if d<14.0 {Some((i,d))}else{None}
        }).min_by(|a,b|a.1.total_cmp(&b.1)).map(|x|x.0)
    }
    pub fn selected_tactical_site(&self)->Option<usize>{
        // A committed campaign battle requires the player's ship to have
        // completed its relativistic navigation leg to the chosen system.
        // The 1/2 sandbox encounters remain playable anywhere.
        self.selected.and_then(|point|self.strategy.nearest(point,500.0))
            .filter(|&site|(self.model.position-self.strategy.systems[site].pos).len()<=500.0)
    }
    pub fn tactical_result(&mut self,site:usize,kind:TacticalKind,
        result:TacticalOutcome){
        self.encounter_notice=Some(match self.strategy.resolve_tactical(site,kind,result){
            Ok(message)=>message,Err(reason)=>reason.into(),
        });
        self.refresh_occupation();
    }
    fn refresh_occupation(&mut self){
        let bucket=(self.model.year/1_000.0).round() as i64;
        if bucket==self.last_front_bucket &&
            self.last_territory_revision==self.strategy.territory_revision{return;}
        let image=blight_image(self.model.front(),&self.strategy);
        self.blight.update(&image);
        self.last_front_bucket=bucket;
        self.last_territory_revision=self.strategy.territory_revision;
    }
    pub fn update(&mut self,seconds:f64){
        let dt=seconds.clamp(0.0,0.10);
        if is_key_pressed(KeyCode::T){self.strategy_visible=!self.strategy_visible;}
        if is_key_pressed(KeyCode::N){
            self.strategy.blight_autonomy=!self.strategy.blight_autonomy;
            self.encounter_notice=Some(format!("BLIGHT AUTONOMY: {}",
                if self.strategy.blight_autonomy{"ENABLED"}else{"PAUSED"}));
        }
        // Site-level industrial actions use the selected navigable waypoint.
        // I: capacity, K: kugelblitz, U: Guf node, P: fire at selected site.
        let action=if is_key_pressed(KeyCode::I){Some(Project::Industry)}
            else if is_key_pressed(KeyCode::K){Some(Project::Kugelblitz)}
            else if is_key_pressed(KeyCode::U){Some(Project::GufNode)}else{None};
        if let Some(project)=action {
            if let Some(site)=self.selected.and_then(|p|self.strategy.nearest(p,500.0)) {
                self.encounter_notice=Some(match self.strategy.project(site,project){
                    Ok(message)=>message,Err(reason)=>reason.into(),
                });
            }else{self.encounter_notice=Some("SELECT A NAMED INDUSTRIAL SYSTEM".into());}
        }
        if is_key_pressed(KeyCode::P){
            if let Some(site)=self.selected.and_then(|p|self.strategy.nearest(p,500.0)) {
                self.encounter_notice=Some(match self.strategy.fire_guf(site){
                    Ok(arrival)=>format!("GUF STRIKE ETA T+{arrival:.0}"),
                    Err(reason)=>reason.into(),
                });
            }
        }
        if is_key_pressed(KeyCode::Y){
            if let Some(site)=self.selected.and_then(|p|self.strategy.nearest(p,500.0)) {
                self.encounter_notice=Some(match self.strategy.launch_halo_strike(site){
                    Ok(arrival)=>format!("HALO STRIKE ETA T+{arrival:.0}"),
                    Err(reason)=>reason.into(),
                });
            }
        }
        if is_key_pressed(KeyCode::V){
            if let Some(site)=self.selected.and_then(|p|self.strategy.nearest(p,500.0)) {
                self.encounter_notice=Some(match self.strategy.claim_for_weena(site){
                    Ok(arrival)=>format!("WEENA INFLUENCE ETA T+{arrival:.0}"),
                    Err(reason)=>reason.into(),
                });
            }
        }
        if is_key_pressed(KeyCode::J){
            if let Some(site)=self.selected.and_then(|p|self.strategy.nearest(p,500.0)) {
                self.encounter_notice=Some(match self.strategy.schedule_blight_seed(
                        self.strategy.systems.iter().enumerate()
                            .filter(|(i,x)|*i!=site&&x.owner==Faction::Blight)
                            .min_by(|a,b|(a.1.pos-self.strategy.systems[site].pos).len()
                                .total_cmp(&(b.1.pos-self.strategy.systems[site].pos).len()))
                            .map(|(i,_)|i).unwrap_or(9),site){
                    Ok(())=>"BLIGHT EXPANSION ORDER IN FLIGHT".into(),
                    Err(reason)=>reason.into(),
                });
            }
        }
        if is_key_pressed(KeyCode::Space){
            self.model.paused=!self.model.paused;
        }
        if is_key_pressed(KeyCode::Q){
            self.model.rate_yr_per_second=
                (self.model.rate_yr_per_second/2.0).max(25.0);
        }
        if is_key_pressed(KeyCode::E){
            self.model.rate_yr_per_second=
                (self.model.rate_yr_per_second*2.0).min(20_000.0);
        }
        if is_key_pressed(KeyCode::Comma){
            self.model.seek(self.model.year-1_000.0);
        }
        if is_key_pressed(KeyCode::Period){
            self.model.seek(self.model.year+1_000.0);
        }
        if is_key_pressed(KeyCode::Home){self.model.seek(0.0);}
        if is_key_pressed(KeyCode::End){
            self.model.seek(100_000.0);
        }
        if is_key_pressed(KeyCode::G){self.grid=!self.grid;}
        if is_key_pressed(KeyCode::L){
            self.local_grid=!self.local_grid;
        }
        if is_key_pressed(KeyCode::B){
            self.show_blight=!self.show_blight;
        }
        if is_key_pressed(KeyCode::H){
            self.show_help=!self.show_help;
        }
        if is_key_pressed(KeyCode::C){
            self.center=V2::default();self.zoom=1.0;
        }
        if is_key_pressed(KeyCode::O){
            self.center=ORION_OUTBREAK;
            self.zoom=2.5;
        }
        if is_key_pressed(KeyCode::F){
            self.focus_frontier();
        }
        if is_key_pressed(KeyCode::Z){
            self.zoom=(self.zoom*0.85).max(0.50);
        }
        if is_key_pressed(KeyCode::X){
            self.zoom=(self.zoom*1.17).min(12.0);
        }
        let (_,wheel)=mouse_wheel();
        if wheel!=0.0 {
            self.zoom=(self.zoom*(1.0+wheel*0.12).max(0.2))
                .clamp(0.5,12.0);
        }
        let pan_speed=11_000.0*dt/
            self.zoom.max(0.5) as f64;
        if is_key_down(KeyCode::A)||is_key_down(KeyCode::Left){
            self.center.x-=pan_speed;
        }
        if is_key_down(KeyCode::D)||is_key_down(KeyCode::Right){
            self.center.x+=pan_speed;
        }
        if is_key_down(KeyCode::W)||is_key_down(KeyCode::Up){
            self.center.y+=pan_speed;
        }
        if is_key_down(KeyCode::S)||is_key_down(KeyCode::Down){
            self.center.y-=pan_speed;
        }
        let mouse=vec2(mouse_position().0,mouse_position().1);
        if is_mouse_button_down(MouseButton::Right) &&
            mouse.x<self.map_width() {
            if let Some(last)=self.dragging {
                let s=self.pixels_per_ly() as f64;
                self.center.x-=(mouse.x-last.x) as f64/s;
                self.center.y+=(mouse.y-last.y) as f64/s;
            }
            self.dragging=Some(mouse);
        }else{self.dragging=None;}
        if is_mouse_button_pressed(MouseButton::Left) &&
            mouse.x<self.map_width() {
            if let Some(index)=self.waypoint_at(mouse){
                self.selected=Some(WAYPOINTS[index].pos);
                self.selected_name=Some(WAYPOINTS[index].name);
            }else{
                self.selected=Some(self.screen_to_world(mouse));
                self.selected_name=None;
            }
            if let Some(selected)=self.selected{
                self.model.choose(selected);
            }
        }
        if is_key_pressed(KeyCode::Enter) {
            if self.model.depart(0.76){
                self.encounter_notice=None;
            }
        }
        if is_key_pressed(KeyCode::Backspace){
            self.model.destination=None;
            self.selected=None;
            self.selected_name=None;
        }
        self.model.tick(dt);
        if self.model.year<self.strategy.year {
            // Historical rewind clears orders but preserves the user's AI setting.
            let autonomy=self.strategy.blight_autonomy;
            self.strategy=StrategicState::new(self.model.year);
            self.strategy.blight_autonomy=autonomy;
        }else{self.strategy.advance_to(self.model.year);}
        if let Some(p)=self.model.last_arrival.take(){
            let near=WAYPOINTS.iter()
                .find(|site|(site.pos-p).len()<500.0);
            self.encounter_notice=Some(match near {
                Some(site)=>format!("ARRIVED: {} / encounter available",
                    site.name),
                None=>"ARRIVED: uncharted sector / survey ready".into(),
            });
        }
        self.refresh_occupation();
    }
    fn draw_textures(&self){
        let top=self.world_to_screen(
            V2::new(-MAP_EXTENT_LY,MAP_EXTENT_LY));
        let side=(GALAXY_DIAMETER_LY*self.pixels_per_ly() as f64) as f32;
        let params=DrawTextureParams{
            dest_size:Some(vec2(side,side)),
            ..Default::default()
        };
        draw_texture_ex(&self.base,top.x,top.y,WHITE,params.clone());
        if self.show_blight{
            draw_texture_ex(&self.blight,top.x,top.y,WHITE,params);
        }
    }
    fn draw_galactic_grid(&self){
        if !self.grid{return;}
        let origin=self.world_to_screen(V2::default());
        let scale=self.pixels_per_ly();
        for r in [15_000.0_f64,30_000.0,45_000.0,60_000.0]{
            let radius=r as f32*scale;
            draw_circle_lines(origin.x,origin.y,radius,0.65,
                Color::new(0.62,0.71,0.83,0.23));
            let label=self.world_to_screen(V2::new(-1_700.0,r));
            if label.x>10.0&&label.x<self.map_width()-60.0&&
                label.y>70.0&&label.y<screen_height()-25.0{
                shadow_text(&format!("{:.0} ly",r),
                    label.x+7.0,label.y,14.0,
                    Color::new(0.70,0.77,0.85,0.54));
            }
        }
        for k in 0..12{
            let angle=k as f64*PI/6.0;
            let direction=V2::new(angle.sin(),angle.cos());
            let endpoint=self.world_to_screen(direction*62_000.0);
            draw_line(origin.x,origin.y,endpoint.x,endpoint.y,
                0.65,Color::new(0.72,0.78,0.87,0.22));
            let p=self.world_to_screen(direction*67_000.0);
            if p.x>25.0&&p.x<self.map_width()-25.0&&
               p.y>87.0&&p.y<screen_height()-20.0{
                shadow_text(&format!("{}",k*30),p.x-13.0,p.y,
                    16.0,Color::new(0.75,0.78,0.84,0.65));
            }
        }
        let p=self.world_to_screen(V2::new(0.0,66_000.0));
        if p.x>30.0&&p.x<self.map_width()-45.0&&p.y>75.0{
            shadow_text("GALACTIC LONGITUDE",p.x-76.0,p.y-23.0,
                14.0,Color::new(0.68,0.75,0.90,0.53));
        }
    }
    fn draw_local_grid(&self){
        let origin=self.world_to_screen(ORION_OUTBREAK);
        let scale=self.pixels_per_ly();
        if self.local_grid{
            for r in (2_000..=12_000).step_by(2_000){
                let color=if r==10_000{
                    Color::new(0.94,0.77,0.63,0.60)
                }else{Color::new(0.71,0.54,0.49,0.22)};
                draw_dashed_circle(origin,r as f32*scale,
                    color,160);
            }
            for k in 0..24{
                let a=k as f64*TAU/24.0;
                let end=self.world_to_screen(
                    ORION_OUTBREAK+V2::from_angle(a)*12_000.0);
                draw_line(origin.x,origin.y,end.x,end.y,0.6,
                    Color::new(0.78,0.46,0.38,0.21));
            }
        }
        let reach=self.model.front().reach_ly();
        if reach>100.0{
            draw_dashed_circle(origin,reach as f32*scale,
                Color::new(0.92,0.79,0.68,0.69),144);
            let annotation=self.world_to_screen(
                ORION_OUTBREAK+V2::from_angle(-0.72)*reach);
            if annotation.x>80.0&&
                annotation.x<self.map_width()-210.0 &&
                annotation.y>120.0&&annotation.y<screen_height()-50.0{
                shadow_text(&format!("CAUSAL REACH / {:.0} ly",
                    reach),annotation.x+5.0,annotation.y-9.0,
                    14.0,BONE);
            }
        }
    }
    fn draw_waypoints(&self){
        for (idx,site) in WAYPOINTS.iter().enumerate(){
            let p=self.world_to_screen(site.pos);
            if p.x<5.0||p.x>self.map_width()-20.0||
                p.y<80.0||p.y>screen_height()-50.0{continue;}
            let selected=self.selected_name==Some(site.name);
            let tone=if idx==0{GOLD}else if idx>=4{CYAN}else{BLUE};
            draw_circle_lines(p.x,p.y,selected.then_some(7.2)
                .unwrap_or(4.8),1.15,tone);
            draw_circle(p.x,p.y,1.6,tone);
            if selected||self.zoom>1.3||idx==0||idx==5 {
                let short=match idx{
                    0=>"SOL",1=>"AMICUS",2=>"OLYMPUS",
                    3=>"GREAT LIBRARY",4=>"OUTER ARM",
                    5=>"BEL-AIR",6=>"NEBULAR WASTES",
                    7=>"PERSEUS FRONT",_=>"CANIS EXODUS",
                };
                shadow_text(short,p.x+9.0,p.y-7.0,14.0,
                    translucent(tone,0.83));
            }
        }
        if self.strategy_visible {
            for site in &self.strategy.systems {
                let p=self.world_to_screen(site.pos);
                if p.x<3.0||p.x>self.map_width()-3.0||p.y<85.0||p.y>screen_height()-45.0{continue;}
                let tone=match site.owner {
                    Faction::Blight=>REDDIM,Faction::Drummers=>CYAN,
                    Faction::Weena=>GOLD,Faction::Independent=>BLUE,
                };
                let radius=7.5+site.industry.min(6.0) as f32*1.1;
                draw_circle_lines(p.x,p.y,radius,1.4,tone);
                if site.array_node {
                    draw_circle_lines(p.x,p.y,radius+6.0,1.0,GOLD);
                }
                if site.kugelblitz>0 {draw_circle(p.x+radius,p.y-radius,3.0,CYAN);}
            }
            for signal in &self.strategy.signals {
                if self.strategy.year<signal.launch_year{continue;}
                let source=self.strategy.systems[signal.source].pos;
                let target=self.strategy.systems[signal.target].pos;
                let a=self.world_to_screen(source);
                let b=self.world_to_screen(target);
                let color=match signal.kind{
                    crate::strategy::SignalKind::ArrayPulse=>GOLD,
                    crate::strategy::SignalKind::HaloStrike=>CYAN,
                    crate::strategy::SignalKind::WeenaClaim=>GOLD,
                    crate::strategy::SignalKind::BlightSeed=>REDDIM,
                };
                draw_line(a.x,a.y,b.x,b.y,1.1,translucent(color,0.55));
                let transit=(self.strategy.year-signal.launch_year)/
                    (signal.arrival_year-signal.launch_year).max(1.0);
                let p=a+(b-a)*transit.clamp(0.0,1.0) as f32;
                draw_circle(p.x,p.y,3.0,color);
            }
        }
        let p=self.world_to_screen(ORION_OUTBREAK);
        draw_circle(p.x,p.y,3.1,Color::new(0.99,0.31,0.24,0.96));
        draw_circle_lines(p.x,p.y,10.0,1.0,
            Color::new(0.91,0.24,0.24,0.45));
    }
    fn draw_travel(&self){
        let pos=self.world_to_screen(self.model.position);
        if let Some(target)=self.model.destination {
            let end=self.world_to_screen(target);
            let delta=end-pos;
            let steps=(delta.length()/12.0).ceil().min(350.0) as usize;
            for i in 0..steps {
                if i%3==2{continue;}
                let a=i as f32/steps.max(1) as f32;
                let b=(i+1) as f32/steps.max(1) as f32;
                draw_line(pos.x+delta.x*a,pos.y+delta.y*a,
                    pos.x+delta.x*b,pos.y+delta.y*b,1.2,
                    Color::new(0.54,0.87,0.98,0.64));
            }
            draw_circle_lines(end.x,end.y,9.0,1.2,CYAN);
            draw_circle(end.x,end.y,2.0,CYAN);
        }
        let heading=if let Some(leg)=self.model.travel{
            (leg.destination-leg.position()).normalized()
        }else{V2::new(0.0,1.0)};
        let forward=vec2(heading.x as f32,-heading.y as f32);
        let normal=vec2(-forward.y,forward.x);
        draw_circle(pos.x,pos.y,13.0,
            Color::new(0.11,0.47,0.68,0.11));
        draw_triangle(pos+forward*11.0,
            pos-forward*6.0+normal*5.0,
            pos-forward*6.0-normal*5.0,CYAN);
        draw_circle_lines(pos.x,pos.y,13.0,1.0,
            Color::new(0.52,0.85,0.96,0.60));
    }
    fn draw_logo(&self){
        let x=24.0;let y=48.0;
        // Four offset passes evoke early-1990s metallic game typography.
        shadow_text("BLIGHTED GALAXY",x+3.0,y+4.0,42.0,
            Color::new(0.42,0.08,0.10,1.0));
        shadow_text("BLIGHTED GALAXY",x,y,42.0,
            Color::new(0.83,0.85,0.86,1.0));
        draw_line(x,y+8.0,x+368.0,y+8.0,2.1,
            Color::new(0.77,0.19,0.20,0.67));
        // Scattered crimson points descending from the title are
        // tiny blighted stars. They should read as blood spatter only
        // from a distance. Sparse on purpose.
        let mut rng=Random::new(0xB100D_5EED);
        for _ in 0..49 {
            let px=x+rng.unit() as f32*365.0;
            let py=y+9.0+rng.unit().powi(3) as f32*52.0;
            let power=0.28+0.55*rng.unit() as f32;
            draw_circle(px,py,0.45+rng.unit() as f32*0.42,
                Color::new(0.89,0.12,0.13,power));
        }
        shadow_text("THE ORION INFECTION / GALACTIC HISTORIOGRAPHY",
            x,y+39.0,13.0,
            Color::new(0.75,0.79,0.85,0.64));
    }
    fn panel_line(&self,text:&str,x:f32,y:f32,tone:Color){
        draw_text(text,x,y,17.0,tone);
    }
    fn draw_sidebar(&self){
        let right=self.map_width();
        let w=screen_width()-right;
        draw_rectangle(right,0.0,w,screen_height(),
            Color::new(0.009,0.016,0.030,0.97));
        draw_line(right,0.0,right,screen_height(),1.2,
            Color::new(0.38,0.51,0.64,0.69));
        let x=right+20.0;
        let mut y=36.0;
        shadow_text("GALACTIC NAVIGATION",x,y,23.0,BONE);
        y+=21.0;
        self.panel_line("JENKINS SWARMS / ORION FRONT",x,y,
            translucent(REDDIM,0.94));
        y+=37.0;
        self.panel_line(&format!("SCENARIO T+{:>7.0} yr",
            self.model.year),x,y,GOLD);
        y+=27.0;
        let radius=self.model.front().reach_ly();
        self.panel_line(&format!("CAUSAL LIMIT   {:>8.0} ly",
            radius),x,y,BONE);
        y+=26.0;
        self.panel_line("EFFECTIVE FRONT       0.10 c",x,y,CYAN);
        y+=25.0;
        self.panel_line(&format!("CLOCK      {:>7.0} yr/s",
            self.model.rate_yr_per_second),x,y,BLUE);
        y+=20.0;
        draw_line(x,y,right+w-16.0,y,0.7,
            Color::new(0.61,0.69,0.79,0.30));
        y+=28.0;
        shadow_text("SELECTED RADIAL SECTOR",x,y,19.0,GOLD);
        y+=26.0;
        let point=self.selected.unwrap_or(self.model.position);
        let sector=overworld::sector_of(point);
        let status=self.strategy.territory_status(point,self.model.front());
        self.panel_line(&format!("LONGITUDE  {:03}-{:03} deg",
            sector.galactic_angle_index*30,
            (sector.galactic_angle_index+1)*30),x,y,BONE);
        y+=25.0;
        self.panel_line(&format!("GC RING    {:02} / 15,000 ly",
            sector.galactic_ring_index),x,y,BONE);
        y+=25.0;
        self.panel_line(&format!("ORION RING {:02} / 2,000 ly",
            sector.local_ring_index),x,y,CYAN);
        y+=25.0;
        self.panel_line(&format!("ORION AZ   {:02} / 24",
            sector.local_angle_index),x,y,CYAN);
        y+=25.0;
        self.panel_line(&format!("FRONT      {}",
            status.name()),x,y,
            match status {
                SectorStatus::Held=>REDDIM,
                SectorStatus::Contested=>GOLD,
                SectorStatus::Defended=>BLUE,
                SectorStatus::Reclaimed=>CYAN,
                SectorStatus::Unaffected=>BONE,
            });
        y+=24.0;
        self.panel_line(&format!("FROM CORE  {:>8.0} ly",
            sector.from_core_ly),x,y,BONE);
        y+=23.0;
        self.panel_line(&format!("FROM ORION {:>8.0} ly",
            sector.from_outbreak_ly),x,y,BONE);
        y+=24.0;
        draw_line(x,y,right+w-16.0,y,0.7,
            Color::new(0.61,0.69,0.79,0.30));
        y+=27.0;
        shadow_text("FLIGHT PLANNER",x,y,19.0,GOLD);
        y+=26.0;
        if let Some(leg)=self.model.travel {
            self.panel_line(&format!("CRUISE  {:.2} c",leg.cruise_beta),
                x,y,CYAN);y+=24.0;
            self.panel_line(&format!("REMAIN  {:>8.0} yr",
                leg.remaining_yr()),x,y,BONE);y+=24.0;
            self.panel_line(&format!("PROPER  {:>8.0} yr",
                leg.proper_elapsed_yr),x,y,CYAN);y+=24.0;
        }else if let Some(dest)=self.model.destination{
            let distance=(dest-self.model.position).len();
            self.panel_line(&format!("DISTANCE {:>8.0} ly",distance),
                x,y,BONE);y+=24.0;
            self.panel_line(&format!("ETA @ .76c {:>8.0} yr",
                distance/0.76),x,y,CYAN);y+=24.0;
            self.panel_line("ENTER  COMMIT RELATIVISTIC LEG",
                x,y,GOLD);y+=24.0;
        }else{
            self.panel_line("CLICK GALAXY TO SELECT DESTINATION",
                x,y,BONE);y+=24.0;
            self.panel_line("ENTER DEPARTS AT 0.76 c",x,y,CYAN);y+=24.0;
        }
        self.panel_line(&format!("TOTAL PROPER {:>8.0} yr",
            self.model.proper_time_total_yr),x,y,BLUE);
        y+=30.0;
        if let Some(title)=self.selected_name {
            self.panel_line(title,x,y,GOLD);y+=21.0;
            if let Some(s)=WAYPOINTS.iter().find(|s|s.name==title){
                self.panel_line(s.source,x,y,BONE);y+=21.0;
            }
        }
        if self.strategy_visible {
            y+=16.0;
            shadow_text("STRATEGIC SYSTEM",x,y,19.0,GOLD);y+=25.0;
            if let Some(index)=self.strategy.nearest(point,500.0){
                let s=&self.strategy.systems[index];
                self.panel_line(s.owner.label(),x,y,CYAN);y+=23.0;
                self.panel_line(&format!("IND {:.1}   E {:.0}  M {:.0}",s.industry,s.energy,s.material),x,y,BONE);y+=23.0;
                self.panel_line(&format!("K-HOLES {}   GUF {}",s.kugelblitz,if s.array_node{"READY"}else{"OFF"}),x,y,GOLD);y+=23.0;
                self.panel_line(if (self.model.position-s.pos).len()<=500.0{
                    "1 / 2: CAMPAIGN BATTLE READY"
                }else{"TRAVEL HERE TO PLAY CAMPAIGN BATTLE"},x,y,CYAN);y+=23.0;
            }
            self.panel_line(&format!("GUF CHARGE {:.0}%  WAVES {}",self.strategy.array_charged,self.strategy.signals.len()),x,y,BONE);y+=22.0;
            self.panel_line(&format!("BLIGHT AUTONOMY {}",if self.strategy.blight_autonomy{"ON"}else{"OFF"}),x,y,REDDIM);
        }
        let footer=screen_height()-160.0;
        if let Some(msg)=&self.encounter_notice{
            draw_text(msg,x,footer-21.0,15.0,CYAN);
        }
        draw_line(x,footer,right+w-16.0,footer,1.0,
            Color::new(0.6,0.7,0.8,0.31));
        let hints=[
            "WASD / RIGHT-DRAG    PAN",
            "WHEEL / Z X        ZOOM",
            "F ORION   C CORE   O FRONT",
            "G GALACTIC / L LOCAL STENCIL",
            "B DARK OCCLUSION   H HELP",
            "SPACE PAUSE   , / . = 1 KYR",
            "T STRATEGY  I INDUSTRY  K KUGELBLITZ",
            "U GUF  V WEENA  J BLIGHT  P FIRE",
            "Y HALO STRIKE  N SWARM AUTONOMY",
        ];
        for (i,line) in hints.iter().enumerate(){
            draw_text(line,x,footer+18.0+i as f32*15.7,
                14.5,Color::new(0.58,0.67,0.78,0.89));
        }
    }
    fn timeline(&self){
        let x=23.0;
        let y=screen_height()-34.0;
        let length=(self.map_width()-65.0).max(200.0);
        draw_rectangle(x-4.0,y-20.0,length+12.0,32.0,
            Color::new(0.009,0.013,0.026,0.70));
        draw_line(x,y,x+length,y,1.0,
            Color::new(0.62,0.71,0.83,0.57));
        for i in 0..=5{
            let px=x+i as f32*length/5.0;
            draw_line(px,y-5.0,px,y+5.0,1.0,
                Color::new(0.65,0.69,0.74,0.50));
            shadow_text(&format!("{}k",i*20),px-11.0,y-8.0,
                14.0,BONE);
        }
        let t=(self.model.year/100_000.0).clamp(0.0,1.0) as f32;
        let px=x+t*length;
        draw_circle(px,y,4.7,Color::new(0.98,0.36,0.32,0.96));
        shadow_text("EPISODES VIII-XXII / DELAYED HISTORY",
            x,y-30.0,15.0,
            Color::new(0.68,0.79,0.89,0.72));
    }
    pub fn draw(&self){
        clear_background(BG);
        self.draw_textures();
        self.draw_galactic_grid();
        self.draw_local_grid();
        self.draw_waypoints();
        self.draw_travel();
        self.draw_logo();
        self.draw_sidebar();
        self.timeline();
        if self.show_help {
            draw_rectangle(42.0,122.0,self.map_width()-84.0,210.0,
                Color::new(0.006,0.016,0.030,0.94));
            draw_rectangle_lines(42.0,122.0,self.map_width()-84.0,
                210.0,1.0,BLUE);
            let text=[
                "THE BLIGHT IS AN INFRARED-DARK STELLAR OCCULTATION.",
                "RED POINTS TRACK ACTIVE ORION-BOUNDED OCCUPATION AND TENDRILS.",
                "THE GALACTIC GRID IS CENTERED ON THE CORE / 30 DEG, 15 KLY.",
                "THE LOCAL FRONTIER STENCIL IS CENTERED ON ORION / 15 DEG, 2 KLY.",
                "WAYPOINT COORDINATES ARE SCHEMATIC; SOURCE EPISODES ARE CANON.",
                "CRUISE TIME AND PROPER TIME OBEY SPECIAL RELATIVITY.",
                "SELECT STAR AND PRESS 1 / 2 TO START A STRATEGIC ENCOUNTER.",
            ];
            for (i,line) in text.iter().enumerate(){
                draw_text(line,61.0,154.0+i as f32*24.0,
                    17.0,if i==1 {GOLD}else{BONE});
            }
        }
        shadow_text("100,000-YEAR SCENARIO / 0.10c FRONT",
            22.0,screen_height()-72.0,15.0,
            Color::new(0.80,0.75,0.73,0.82));
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn world_texture_projection_preserves_orion_geometry(){
        let center=V2::default();
        let (x,y)=tex_xy(center);
        assert_eq!(x,TEX_SIZE as f64/2.0);
        assert_eq!(y,TEX_SIZE as f64/2.0);
        let (a,b)=tex_xy(ORION_OUTBREAK);
        assert!((a-x).abs()<15.0);
        assert!(b>y);
    }
    #[test]
    fn procedural_sector_limits_are_physical(){
        let f=Front::new(100_000.0);
        for i in 0..360{
            let t=i as f64*TAU/360.0;
            assert!(f.sector_reach_ly(t)<=10_000.0);
        }
    }
}
