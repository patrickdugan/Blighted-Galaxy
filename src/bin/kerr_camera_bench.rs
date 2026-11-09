//! Headless Kerr camera smoke test and optional binary PPM export.
 //! cargo run --release --bin kerr_camera_bench -- --output kerr-camera.png
use blighted_galaxy::kerr::KerrGeometry;
use blighted_galaxy::kerr_camera::KerrCamera;
use blighted_galaxy::physics::SOLAR_MASS_KG;
use std::f64::consts::PI;

fn main() {
    let output=std::env::args().skip(1).collect::<Vec<_>>();
    let path=output.windows(2).find(|x|x[0]=="--output")
        .map(|x|x[1].as_str());
    let mass=10.0*SOLAR_MASS_KG;
    let mut previous=None;
    for spin in [0.0,0.85] {
        let geometry=KerrGeometry::new(mass,spin).unwrap();
        let camera=KerrCamera::new(geometry,85.0,PI/3.0,36.0).unwrap();
        let image=camera.render_rgba(120,80);
        println!("Kerr camera a={spin:+.2}, inclination 60 degrees");
        println!("  pixels={} captured={} disc={} sky={} unresolved={}",
            image.width as usize*image.height as usize,
            image.captured,image.disc,image.sky,image.unresolved);
        println!("  average integration steps {:.2}; radial/polar worst {:.3e}/{:.3e}",
            image.mean_steps_per_ray(),
            image.worst_radial_residual,image.worst_polar_residual);
        assert!(image.captured>0 && image.disc>0 && image.sky>0);
        assert!(image.worst_radial_residual < 1e-4);
        assert!(image.worst_polar_residual < 1e-3);
        assert_eq!(image.pixels.len(),120*80*4);
        if let Some(prior)=previous.take() {
            assert_ne!(prior,image.pixels,"rotating Kerr must change sky image");
        }
        if spin>0.0 {
            if let Some(p)=path {
                if p.ends_with(".png") {
                    let file=std::fs::File::create(p).expect("create PNG");
                    let mut encoder=png::Encoder::new(
                        std::io::BufWriter::new(file),
                        image.width as u32,image.height as u32);
                    encoder.set_color(png::ColorType::Rgba);
                    encoder.set_depth(png::BitDepth::Eight);
                    let mut writer=encoder.write_header().expect("PNG header");
                    writer.write_image_data(&image.pixels).expect("PNG pixels");
                } else {
                    std::fs::write(p,image.ppm_bytes()).expect("write PPM");
                }
                println!("  wrote {p}");
            }
        } else {
            previous=Some(image.pixels);
        }
    }
}
