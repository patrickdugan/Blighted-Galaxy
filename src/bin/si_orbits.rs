//! SI-scale graphical orbital physics lab.
 //!
 //! Run: cargo run --release --bin si_orbits
 //! Real masses, distances and time steps. Body disks are schematic markers.
use blighted_galaxy::physics::{
    circular_binary_si, evolve_nbody_si, light_travel_time_s,
    near_side_tidal_ratio_si, newtonian_energy_j, retarded_emission_si,
    schwarzschild_radius_m, weak_lensing_deflection_rad,
    SiBody, AU_M, C_SI, EARTH_MASS_KG, G_SI, JULIAN_YEAR_S, SOLAR_MASS_KG,
};
use blighted_galaxy::V2;
use macroquad::prelude::*;

const SUN_RADIUS_M: f64 = 6.957e8;
const EARTH_RADIUS_M: f64 = 6.371e6;
const WHITE_DWARF_RADIUS_M: f64 = 9.0e6;
const MAX_TRAIL: usize = 720;

fn window_conf() -> Conf {
    Conf {
        window_title: "Blighted Galaxy: SI Orbital Lab".into(),
        window_width: 1500,
        window_height: 900,
        high_dpi: true,
        ..Default::default()
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Scenario {
    SunEarth,
    DwarfFlyby,
    BlackHoleBinary,
}
impl Scenario {
    fn heading(self) -> &'static str {
        match self {
            Self::SunEarth => "1 / SUN-EARTH: REGRESSION ORBIT",
            Self::DwarfFlyby => "2 / WHITE DWARF: CLOSE-ENCOUNTER PERTURBATION",
            Self::BlackHoleBinary => "3 / BLACK HOLE: BINARY AND INNER WORLD",
        }
    }
    fn notes(self) -> &'static str {
        match self {
            Self::SunEarth => "Reference: circular Newtonian binary, 1 AU, real masses",
            Self::DwarfFlyby => "Flyby: energy transfer and tidal stresses; material fracture not simulated",
            Self::BlackHoleBinary => "Weak-field orbit of an 8-solar-mass compact object; disk radius is schematic",
        }
    }
}
struct Lab {
    scenario: Scenario,
    bodies: Vec<SiBody>,
    names: Vec<&'static str>,
    tones: Vec<Color>,
    histories: Vec<Vec<V2>>,
    elapsed_s: f64,
    reference_energy_j: f64,
    paused: bool,
    days_per_real_second: f64,
    camera_au: V2,
    pixels_per_au: f32,
}
impl Lab {
    fn new(scenario: Scenario) -> Self {
        let (bodies, names, tones, scale) = match scenario {
            Scenario::SunEarth => {
                let mut binary = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
                binary[0].radius_m = SUN_RADIUS_M;
                binary[1].radius_m = EARTH_RADIUS_M;
                (binary.to_vec(), vec!["Sun", "Earth"], vec![GOLD, SKYBLUE], 270.0)
            }
            Scenario::DwarfFlyby => {
                let mut binary = circular_binary_si(SOLAR_MASS_KG, EARTH_MASS_KG, AU_M);
                binary[0].radius_m = SUN_RADIUS_M;
                binary[1].radius_m = EARTH_RADIUS_M;
                let white_dwarf = SiBody::new(
                    0.65 * SOLAR_MASS_KG, WHITE_DWARF_RADIUS_M,
                    V2::new(-4.8 * AU_M, 1.6 * AU_M),
                    V2::new(15_000.0, -8_000.0),
                );
                (vec![binary[0],binary[1],white_dwarf],
                    vec!["Sun", "Earth", "White dwarf"],
                    vec![GOLD, SKYBLUE, LIGHTGRAY], 110.0)
            }
            Scenario::BlackHoleBinary => {
                let mut binary = circular_binary_si(8.0 * SOLAR_MASS_KG, SOLAR_MASS_KG, 2.5 * AU_M);
                binary[0].radius_m = schwarzschild_radius_m(8.0 * SOLAR_MASS_KG).unwrap();
                binary[1].radius_m = SUN_RADIUS_M;
                let planet_offset = 0.3 * AU_M;
                let planet = SiBody::new(
                    EARTH_MASS_KG, EARTH_RADIUS_M,
                    binary[1].pos_m + V2::new(planet_offset, 0.0),
                    binary[1].vel_m_s + V2::new(0.0, (G_SI * SOLAR_MASS_KG / planet_offset).sqrt()),
                );
                (vec![binary[0], binary[1], planet],
                    vec!["8-Msun black hole", "Companion", "Inner planet"],
                    vec![ORANGE, GOLD, SKYBLUE], 135.0)
            }
        };
        let energy = newtonian_energy_j(&bodies, 0.0);
        let histories = bodies.iter().map(|body| vec![body.pos_m]).collect();
        Self {
            scenario, bodies, names, tones, histories,
            elapsed_s: 0.0, reference_energy_j: energy,
            paused: false, days_per_real_second: 22.0,
            camera_au: V2::default(),
            pixels_per_au: scale,
        }
    }
    fn frame_advance(&mut self, frame_seconds: f64) {
        if self.paused { return; }
        let elapsed = frame_seconds.min(0.05) * self.days_per_real_second * 86_400.0;
        // Fixed SI time resolution: 3 h on all bodies. A true close encounter
        // requires local adaptivity (later milestone).
        let max_step_s = 3.0 * 3600.0;
        evolve_nbody_si(&mut self.bodies, elapsed, max_step_s, 0.0);
        self.elapsed_s += elapsed;
        for (trail, body) in self.histories.iter_mut().zip(self.bodies.iter()) {
            trail.push(body.pos_m);
            if trail.len() > MAX_TRAIL {
                trail.drain(..(trail.len() - MAX_TRAIL));
            }
        }
    }
    fn px(&self, pos_m: V2) -> Vec2 {
        let au = pos_m / AU_M;
        vec2(
            (au.x - self.camera_au.x) as f32 * self.pixels_per_au + screen_width() * 0.5,
            (au.y - self.camera_au.y) as f32 * self.pixels_per_au + screen_height() * 0.5,
        )
    }
    fn draw_grid(&self) {
        let span = 18;
        for i in -span..=span {
            let x = (i as f64 - self.camera_au.x) as f32 * self.pixels_per_au + screen_width() * 0.5;
            let y = (i as f64 - self.camera_au.y) as f32 * self.pixels_per_au + screen_height() * 0.5;
            draw_line(x, 0.0, x, screen_height(), 1.0,
                Color::new(0.19, 0.26, 0.37, 0.15));
            draw_line(0.0, y, screen_width(), y, 1.0,
                Color::new(0.19, 0.26, 0.37, 0.15));
        }
        let one_au = self.px(V2::new(AU_M, 0.0));
        let center = self.px(V2::default());
        if self.pixels_per_au > 30.0 {
            draw_circle_lines(center.x, center.y, self.pixels_per_au, 1.0,
                Color::new(0.25, 0.5, 0.65, 0.20));
            draw_text("1 AU", one_au.x, one_au.y - 8.0, 15.0, GRAY);
        }
    }
    fn draw(&self) {
        clear_background(Color::new(0.008, 0.016, 0.040, 1.0));
        self.draw_grid();

        for (i, trail) in self.histories.iter().enumerate() {
            let tone = self.tones[i];
            for segment in trail.windows(2) {
                let a = self.px(segment[0]);
                let b = self.px(segment[1]);
                draw_line(a.x, a.y, b.x, b.y, 1.35,
                    Color::new(tone.r, tone.g, tone.b, 0.6));
            }
        }

        for (i, body) in self.bodies.iter().enumerate() {
            let p = self.px(body.pos_m);
            if p.x < -100.0 || p.y < -100.0 ||
                p.x > screen_width() + 100.0 || p.y > screen_height() + 100.0 { continue; }
            // Visible sizes are visualization markers; physical radii are used
            // in the tide diagnostic and are not remotely to scale on this map.
            let radius = match i {
                0 if self.scenario == Scenario::BlackHoleBinary => 11.0,
                0 => 17.0,
                1 if self.scenario != Scenario::BlackHoleBinary => 5.0,
                1 => 12.0,
                _ => 6.0,
            };
            let tone = self.tones[i];
            for n in 1..=4 {
                draw_circle(p.x, p.y, radius * (1.0 + n as f32 * 0.7),
                    Color::new(tone.r, tone.g, tone.b, 0.035));
            }
            if self.scenario == Scenario::BlackHoleBinary && i == 0 {
                draw_circle(p.x, p.y, radius + 4.0, ORANGE);
                draw_circle(p.x, p.y, radius, BLACK);
                draw_circle_lines(p.x, p.y, radius + 2.0, 1.0, YELLOW);
            } else {
                draw_circle(p.x, p.y, radius, tone);
            }
            let velocity_marker = self.px(body.pos_m + body.vel_m_s * 70_000.0);
            draw_line(p.x, p.y, velocity_marker.x, velocity_marker.y,
                1.2, Color::new(tone.r, tone.g, tone.b, 0.60));
            draw_text(self.names[i], p.x + radius + 10.0, p.y - radius - 8.0, 19.0, tone);
        }

        draw_rectangle(16.0, 16.0, 770.0, 191.0, Color::new(0.01, 0.03, 0.06, 0.9));
        draw_text("BLIGHTED GALAXY / CALIBRATED SI PHYSICS LAB", 29.0, 47.0, 26.0, WHITE);
        draw_text(self.scenario.heading(), 29.0, 80.0, 22.0, SKYBLUE);
        draw_text(self.scenario.notes(), 29.0, 108.0, 17.0, LIGHTGRAY);
        draw_text("1 / 2 / 3 scenarios  |  SPACE pause  |  Q/E time warp  |  ARROWS pan", 29.0, 139.0, 18.0, LIGHTGRAY);
        draw_text("Z/X zoom  |  R reset  |  Esc exit  |  marker radii are schematic", 29.0, 165.0, 18.0, LIGHTGRAY);
        draw_text(&format!("SIM TIME {:.4} yr  |  {:.1} days / real second  |  {}",
            self.elapsed_s / JULIAN_YEAR_S, self.days_per_real_second,
            if self.paused { "PAUSED" } else { "RUNNING" }),
            29.0, 193.0, 17.0, YELLOW);

        let e = newtonian_energy_j(&self.bodies, 0.0);
        let drift = (e - self.reference_energy_j).abs() / self.reference_energy_j.abs().max(1.0);
        let second = &self.bodies[1];
        let first = &self.bodies[0];
        let dist = (second.pos_m - first.pos_m).len();
        let propagation = light_travel_time_s(dist).unwrap_or(0.0);
        let observed = retarded_emission_si(second.pos_m, second.vel_m_s,
            first.pos_m, 0.0);
        let retarded_delay_s = observed.map(|x| x.delay_s).unwrap_or(0.0);
        let tide = near_side_tidal_ratio_si(first, second);
        let black_hole_lensing = if self.scenario == Scenario::BlackHoleBinary {
            weak_lensing_deflection_rad(first.mass_kg, 0.25 * AU_M)
        } else { None };
        let metrics = [
            format!("|delta E|/|E0|   {:10.3e}", drift),
            format!("SEPARATION       {:10.3} AU", dist/AU_M),
            format!("LIGHT DELAY      {:10.2} s", propagation),
            format!("RETARDED DELAY   {:10.2} s", retarded_delay_s),
            format!("TIDAL / SELF G   {:10.3e}", tide.unwrap_or(f64::NAN)),
            format!("MOTION 2         {:10.3} km/s", second.vel_m_s.len()/1000.0),
        ];
        let y_top = screen_height() - 221.0;
        draw_rectangle(16.0, y_top-25.0, 440.0, 210.0,
            Color::new(0.01, 0.03, 0.06, 0.87));
        draw_text("MEASURED DIAGNOSTICS", 29.0, y_top+2.0, 19.0, SKYBLUE);
        for (i, line) in metrics.iter().enumerate() {
            draw_text(line, 29.0, y_top+29.0+i as f32*24.0, 17.0, WHITE);
        }
        if let Some(theta) = black_hole_lensing {
            draw_text(&format!("BH weak-field 0.25 AU: {:.4e} rad", theta),
                29.0, y_top+181.0, 16.0, ORANGE);
        }

        draw_text("Newtonian particle trajectories; SR delays/energies; no Kerr/GR or hydrodynamics",
            24.0, screen_height()-15.0, 18.0, GRAY);
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut lab = Lab::new(Scenario::SunEarth);
    loop {
        if is_key_pressed(KeyCode::Escape) { break; }
        let selected = if is_key_pressed(KeyCode::Key1) { Some(Scenario::SunEarth) }
            else if is_key_pressed(KeyCode::Key2) { Some(Scenario::DwarfFlyby) }
            else if is_key_pressed(KeyCode::Key3) { Some(Scenario::BlackHoleBinary) }
            else { None };
        if let Some(s) = selected { lab = Lab::new(s); }
        if is_key_pressed(KeyCode::R) { lab = Lab::new(lab.scenario); }
        if is_key_pressed(KeyCode::Space) { lab.paused = !lab.paused; }
        if is_key_pressed(KeyCode::Q) {
            lab.days_per_real_second = (lab.days_per_real_second / 2.0).max(0.05);
        }
        if is_key_pressed(KeyCode::E) {
            lab.days_per_real_second = (lab.days_per_real_second * 2.0).min(100.0);
        }
        let dt = get_frame_time().min(0.1) as f64;
        let pan = 1.5 * dt / (lab.pixels_per_au as f64 / 100.0).max(0.1);
        if is_key_down(KeyCode::Left) { lab.camera_au.x -= pan; }
        if is_key_down(KeyCode::Right) { lab.camera_au.x += pan; }
        if is_key_down(KeyCode::Up) { lab.camera_au.y -= pan; }
        if is_key_down(KeyCode::Down) { lab.camera_au.y += pan; }
        let (_, scroll) = mouse_wheel();
        if is_key_down(KeyCode::Z) { lab.pixels_per_au *= 0.98; }
        if is_key_down(KeyCode::X) { lab.pixels_per_au *= 1.02; }
        if scroll.abs() > 0.0 { lab.pixels_per_au *= 1.0 + 0.06*scroll; }
        lab.pixels_per_au = lab.pixels_per_au.clamp(10.0, 2500.0);
        lab.frame_advance(dt);
        lab.draw();
        next_frame().await;
    }
}
