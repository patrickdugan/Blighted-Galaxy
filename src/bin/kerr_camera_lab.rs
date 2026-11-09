//! Interactive off-equatorial Kerr photon camera.
 //! cargo run --release --bin kerr_camera_lab
 //! The source sky and thin disc colors are synthetic visualization rules.
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::{CameraImage,KerrCamera};
use blighted_galaxy::physics::SOLAR_MASS_KG;
use macroquad::prelude::*;
use std::f64::consts::PI;

fn config()->Conf {
    Conf {window_title:"Blighted Galaxy • Kerr 3D Photon Camera".into(),
        window_width:1600,window_height:950,
        high_dpi:true,..Default::default()}
}
struct ImageFrame {
    image:CameraImage,
    texture:Texture2D,
    spin:f64,
    inclination_deg:f64,
    span_rg:f64,
    seconds:f64,
}
#[macroquad::main(config)]
async fn main() {
    let spins=[0.0_f64,0.45,0.85,-0.85];
    let mut selected=2_usize;
    let mut inclination_deg=60.0_f64;
    let mut span_rg=36.0_f64;
    let mut frame:Option<ImageFrame>=None;
    loop {
        if is_key_pressed(KeyCode::Escape) {break;}
        let mut changed=frame.is_none();
        for (key,index) in [
            (KeyCode::Key1,0),(KeyCode::Key2,1),
            (KeyCode::Key3,2),(KeyCode::Key4,3)
        ] {
            if is_key_pressed(key) {selected=index;changed=true;}
        }
        if is_key_pressed(KeyCode::Up) {
            inclination_deg=(inclination_deg+5.0).min(175.0);
            changed=true;
        }
        if is_key_pressed(KeyCode::Down) {
            inclination_deg=(inclination_deg-5.0).max(5.0);
            changed=true;
        }
        if is_key_pressed(KeyCode::Z) {
            span_rg=(span_rg*1.2).min(75.0);
            changed=true;
        }
        if is_key_pressed(KeyCode::X) {
            span_rg=(span_rg/1.2).max(12.0);
            changed=true;
        }
        if is_key_pressed(KeyCode::R) {changed=true;}
        if changed {
            let spin=spins[selected];
            let bh=KerrGeometry::new(10.0*SOLAR_MASS_KG,spin).unwrap();
            let cam=KerrCamera::new(bh,85.0,
                inclination_deg*PI/180.0,span_rg).unwrap();
            let now=get_time();
            let image=cam.render_rgba(240,144);
            let elapsed=get_time()-now;
            let texture=Texture2D::from_rgba8(image.width,image.height,&image.pixels);
            texture.set_filter(FilterMode::Linear);
            frame=Some(ImageFrame {
                image,texture,spin,inclination_deg,span_rg,
                seconds:elapsed,
            });
        }
        clear_background(BLACK);
        if let Some(f)=frame.as_ref() {
            draw_texture_ex(&f.texture,0.0,0.0,WHITE,
                DrawTextureParams {dest_size:Some(vec2(screen_width(),screen_height())),
                    ..Default::default()});
            draw_rectangle(12.0,12.0,1150.0,181.0,
                Color::new(0.008,0.015,0.034,0.91));
            draw_text("BLIGHTED GALAXY / NON-EQUATORIAL KERR CAMERA",
                28.0,48.0,27.0,WHITE);
            draw_text("1 Schwarzschild  2 spin +0.45  3 spin +0.85  4 spin -0.85",
                28.0,78.0,20.0,SKYBLUE);
            draw_text("UP/DOWN inclination    Z/X field of view    R recompute    Esc exit",
                28.0,108.0,19.0,LIGHTGRAY);
            draw_text(&format!(
                "spin={:+.2}  inclination={:.1} deg  FOV span={:.1} rg  render={:.3}s",
                f.spin,f.inclination_deg,f.span_rg,f.seconds),
                28.0,140.0,19.0,GOLD);
            draw_text(&format!(
                "rays: horizon {} | disc {} | sky {} | unresolved {} | mean steps {:.1}",
                f.image.captured,f.image.disc,f.image.sky,
                f.image.unresolved,f.image.mean_steps_per_ray()),
                28.0,173.0,18.0,LIGHTGRAY);
            draw_rectangle(12.0,screen_height()-57.0,
                screen_width()-24.0,47.0,
                Color::new(0.008,0.015,0.034,0.86));
            draw_text("Past-directed ZAMO rays with Kerr Carter constant; synthetic sky & disc; no radiative transfer",
                26.0,screen_height()-27.0,18.0,WHITE);
        }
        next_frame().await;
    }
}
