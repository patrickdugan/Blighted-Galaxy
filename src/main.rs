//! Main executable: SI-driven combat, tidal encounter, and M3 contact boxing lab.
//! Legacy normalized prototype remains available via --bin arcade_legacy.
use blighted_galaxy::campaign_save;
use blighted_galaxy::mission::MissionProfile;
use blighted_galaxy::morph_boxing::{
    can_enter_contact_boxing, BoxingInput, MorphBoxingArena, MorphDuelist, MorphForm,
    ARENA_HALF_EXTENT_M,
};
use blighted_galaxy::overworld_nav::OverworldNav;
use blighted_galaxy::physics::{AU_M, C_SI, SOLAR_MASS_KG};
use blighted_galaxy::si_bridge::{
    CameraMeters, Control, Metres, Seconds, SiJoust, SiTidalEncounter, WeaponKind,
};
use blighted_galaxy::strategy::{Faction, TacticalKind, TacticalOutcome};
use blighted_galaxy::V2;
use macroquad::prelude::*;
use std::path::PathBuf;

fn config() -> Conf {
    Conf {
        window_title: "Blighted Galaxy • Galactic Overworld".into(),
        window_width: 1600,
        window_height: 950,
        high_dpi: true,
        ..Default::default()
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Intercept,
    Tidal,
    Boxing,
}
struct Effect {
    pos_m: V2,
    start_s: f64,
    magnitude: f64,
    tone: Color,
}
struct Game {
    mode: Mode,
    intercept: SiJoust,
    tidal: SiTidalEncounter,
    boxing: MorphBoxingArena,
    camera: CameraMeters,
    sim_rate: f64,
    accumulator: f64,
    manual_aim: bool,
    facing_rad: f64,
    effects: Vec<Effect>,
    displayed_latest_time_s: f64,
    mission_name: Option<&'static str>,
    mission_profile: Option<MissionProfile>,
    tidal_tug_budget: u8,
    demo_preview: bool,
    boxing_requested_form: Option<MorphForm>,
    boxing_strike_queued: bool,
    boxing_pulse_queued: bool,
}
const DELTA: f64 = 1.0 / 120.0;
const TIDAL_STEP: f64 = 0.10;
const BOXING_STEP: f64 = 1.0 / 120.0;

impl Game {
    fn new(mode: Mode) -> Self {
        let (center, scale) = match mode {
            Mode::Intercept => (V2::new(AU_M, 0.0), Metres(550_000.0)),
            Mode::Tidal => (V2::default(), Metres(350_000.0)),
            Mode::Boxing => (V2::default(), Metres(0.28)),
        };
        Self {
            mode,
            intercept: SiJoust::new(),
            tidal: SiTidalEncounter::new(),
            boxing: MorphBoxingArena::default(),
            camera: CameraMeters::new(center, scale),
            sim_rate: 1.0,
            accumulator: 0.0,
            manual_aim: false,
            facing_rad: 0.0,
            effects: Vec::new(),
            displayed_latest_time_s: -1.0,
            mission_name: None,
            mission_profile: None,
            tidal_tug_budget: 0,
            demo_preview: false,
            boxing_requested_form: None,
            boxing_strike_queued: false,
            boxing_pulse_queued: false,
        }
    }
    fn for_mission(profile: MissionProfile) -> Self {
        let mode = match profile.kind {
            TacticalKind::Intercept
            | TacticalKind::GufDefense
            | TacticalKind::ConvoyEscort
            | TacticalKind::HaloStrike => Mode::Intercept,
            TacticalKind::Tidal => Mode::Tidal,
        };
        let mut game = Self::new(mode);
        match mode {
            Mode::Intercept => {
                game.intercept = profile.joust();
                if profile.kind == TacticalKind::HaloStrike {
                    game.camera.metres_per_pixel = Metres(1_250_000.0);
                }
            }
            Mode::Tidal => game.tidal = profile.tidal(),
            Mode::Boxing => {}
        }
        game.mission_name = Some(profile.site_name);
        game.tidal_tug_budget = profile.available_tugs();
        game.mission_profile = Some(profile);
        game
    }
    fn for_boxing(site: Option<usize>, site_name: Option<&'static str>) -> Self {
        let mut game = Self::new(Mode::Boxing);
        game.boxing = MorphBoxingArena::new(site, site_name);
        game.mission_name = site_name;
        game
    }
    fn for_boxing_demo() -> Self {
        let mut game = Self::new(Mode::Boxing);
        game.boxing = MorphBoxingArena::demo_showcase();
        game.mission_name = game.boxing.site_name;
        game.demo_preview = true;
        game
    }
    fn screen(&self, pos: V2) -> Vec2 {
        let p = self
            .camera
            .screen_of(pos, screen_width() as f64, screen_height() as f64);
        vec2(p.x as f32, p.y as f32)
    }
    fn current_time(&self) -> f64 {
        match self.mode {
            Mode::Intercept => self.intercept.time_s,
            Mode::Tidal => self.tidal.time_s,
            Mode::Boxing => self.boxing.time_s,
        }
    }
    fn tactical_outcome(&self) -> TacticalOutcome {
        match self.mode {
            Mode::Intercept => {
                if !self.intercept.player.alive() {
                    TacticalOutcome::Defeat
                } else if !self.intercept.opponent.alive()
                    && self.intercept.escorts.iter().all(|e| !e.craft.alive())
                {
                    TacticalOutcome::Victory
                } else {
                    TacticalOutcome::Withdrawn
                }
            }
            Mode::Tidal => {
                if self.tidal.disruption.is_some() {
                    TacticalOutcome::Victory
                } else if self.mission_profile.is_some() && self.tidal.time_s >= 300.0 {
                    TacticalOutcome::Defeat
                } else {
                    TacticalOutcome::Withdrawn
                }
            }
            Mode::Boxing => self.boxing.outcome(),
        }
    }
    fn input(&mut self) {
        if self.mode != Mode::Boxing && self.mission_name.is_none() {
            if is_key_pressed(KeyCode::Key1) {
                *self = Self::new(Mode::Intercept);
                return;
            }
            if is_key_pressed(KeyCode::Key2) {
                *self = Self::new(Mode::Tidal);
                return;
            }
        }
        if is_key_pressed(KeyCode::R) {
            let preview = self.demo_preview;
            *self = match self.mode {
                Mode::Boxing => {
                    if preview {
                        Self::for_boxing_demo()
                    } else {
                        Self::for_boxing(self.boxing.site, self.boxing.site_name)
                    }
                }
                _ => match self.mission_profile {
                    Some(profile) => Self::for_mission(profile),
                    None => Self::new(self.mode),
                },
            };
            self.demo_preview = preview;
            return;
        }
        if self.mode == Mode::Boxing {
            if is_key_pressed(KeyCode::Key1) {
                self.boxing_requested_form = Some(MorphForm::Streamlined);
            } else if is_key_pressed(KeyCode::Key2) {
                self.boxing_requested_form = Some(MorphForm::AegisGuard);
            } else if is_key_pressed(KeyCode::Key3) {
                self.boxing_requested_form = Some(MorphForm::StrikerLimb);
            } else if is_key_pressed(KeyCode::Key4) {
                self.boxing_requested_form = Some(MorphForm::SwarmErode);
            } else if is_key_pressed(KeyCode::Key5) {
                self.boxing_requested_form = Some(MorphForm::SingularityCore);
            }
            if is_key_pressed(KeyCode::Space) {
                self.boxing_strike_queued = true;
            }
            if is_key_pressed(KeyCode::F) || is_key_pressed(KeyCode::E) {
                self.boxing_pulse_queued = true;
            }
            let (_, scroll) = mouse_wheel();
            if scroll != 0.0 {
                self.camera.zoom((1.0 + scroll as f64 * 0.1).clamp(0.4, 2.5));
            }
            if is_key_down(KeyCode::X) {
                self.camera.zoom(1.03);
            }
            if is_key_down(KeyCode::Z) {
                self.camera.zoom(1.0 / 1.03);
            }
            let (mx, my) = mouse_position();
            let cursor = self.camera.world_of(
                V2::new(mx as f64, my as f64),
                screen_width() as f64,
                screen_height() as f64,
            );
            let toward = cursor - self.boxing.player.pos;
            if toward.len() > 1e-3 {
                self.facing_rad = toward.y.atan2(toward.x);
            }
            return;
        }

        if is_key_pressed(KeyCode::Q) {
            self.sim_rate = (self.sim_rate / 2.0).max(0.125);
        }
        if is_key_pressed(KeyCode::E) {
            self.sim_rate = (self.sim_rate * 2.0).min(64.0);
        }
        if is_key_pressed(KeyCode::M) {
            self.manual_aim = !self.manual_aim;
        }
        if self.mode == Mode::Intercept {
            if is_key_pressed(KeyCode::Key7) && can_enter_contact_boxing(&self.intercept) {
                let site = self.mission_profile.map(|p| p.site);
                let name = self.mission_name;
                let preview = self.demo_preview;
                *self = Self::for_boxing(site, name);
                self.demo_preview = preview;
                return;
            }
            let unlocks = self
                .mission_profile
                .map(|p| (p.siege_unlocked, p.fan_unlocked))
                .unwrap_or((true, true));
            if is_key_pressed(KeyCode::Key4) {
                self.intercept.player_weapon = WeaponKind::Standard;
            } else if is_key_pressed(KeyCode::Key5) && unlocks.0 {
                self.intercept.player_weapon = WeaponKind::Siege;
            } else if is_key_pressed(KeyCode::Key6) && unlocks.1 {
                self.intercept.player_weapon = WeaponKind::Fan;
            }
        }
        if self.mode == Mode::Tidal && self.mission_profile.is_some() {
            if is_key_pressed(KeyCode::A) {
                MissionProfile::apply_tidal_tug(
                    &mut self.tidal,
                    true,
                    &mut self.tidal_tug_budget,
                );
            } else if is_key_pressed(KeyCode::D) {
                MissionProfile::apply_tidal_tug(
                    &mut self.tidal,
                    false,
                    &mut self.tidal_tug_budget,
                );
            }
        }
        let (_, scroll) = mouse_wheel();
        if scroll != 0.0 {
            self.camera.zoom((1.0 + scroll as f64 * 0.1).clamp(0.2, 3.0));
        }
        if is_key_down(KeyCode::X) {
            self.camera.zoom(1.03);
        }
        if is_key_down(KeyCode::Z) {
            self.camera.zoom(1.0 / 1.03);
        }
        if self.mode == Mode::Intercept {
            if self.manual_aim {
                if is_key_down(KeyCode::A) {
                    self.facing_rad -= 0.045;
                }
                if is_key_down(KeyCode::D) {
                    self.facing_rad += 0.045;
                }
            } else {
                let (mx, my) = mouse_position();
                let cursor = self.camera.world_of(
                    V2::new(mx as f64, my as f64),
                    screen_width() as f64,
                    screen_height() as f64,
                );
                let toward = cursor - self.intercept.player.pos_m;
                if toward.len() > 0.0 {
                    self.facing_rad = toward.y.atan2(toward.x);
                }
            }
        }
    }
    fn advance(&mut self, real_dt: f64) {
        let step = match self.mode {
            Mode::Intercept => DELTA,
            Mode::Tidal => TIDAL_STEP,
            Mode::Boxing => BOXING_STEP,
        };
        self.accumulator = (self.accumulator + real_dt.min(0.05) * self.sim_rate).min(8.0);
        let mut cycles = 0;
        while self.accumulator >= step && cycles < 400 {
            match self.mode {
                Mode::Intercept => {
                    let thrust = if is_key_down(KeyCode::W) {
                        0.02 * C_SI
                    } else if is_key_down(KeyCode::S) {
                        -0.01 * C_SI
                    } else {
                        0.0
                    };
                    self.intercept.step(
                        Seconds(step),
                        Control {
                            player_thrust_m_s2: thrust,
                            player_facing_rad: self.facing_rad,
                            fire_player: is_key_down(KeyCode::Space),
                        },
                    );
                }
                Mode::Tidal => {
                    // Campaign interdiction has a finite escape deadline.
                    if self.mission_profile.is_none()
                        || self.tidal.time_s < 300.0
                        || self.tidal.disruption.is_some()
                    {
                        self.tidal.step(Seconds(step));
                    }
                }
                Mode::Boxing => {
                    let mut move_dir = V2::default();
                    if is_key_down(KeyCode::W) || is_key_down(KeyCode::Up) {
                        move_dir.y -= 1.0;
                    }
                    if is_key_down(KeyCode::S) || is_key_down(KeyCode::Down) {
                        move_dir.y += 1.0;
                    }
                    if is_key_down(KeyCode::A) || is_key_down(KeyCode::Left) {
                        move_dir.x -= 1.0;
                    }
                    if is_key_down(KeyCode::D) || is_key_down(KeyCode::Right) {
                        move_dir.x += 1.0;
                    }
                    let req = self.boxing_requested_form.take();
                    let strike = self.boxing_strike_queued || is_key_down(KeyCode::Space);
                    let pulse = self.boxing_pulse_queued;
                    self.boxing_strike_queued = false;
                    self.boxing_pulse_queued = false;
                    self.boxing.step(
                        step,
                        BoxingInput {
                            move_dir,
                            facing_rad: self.facing_rad,
                            requested_form: req,
                            strike,
                            memetic_pulse: pulse,
                        },
                    );
                }
            }
            self.accumulator -= step;
            cycles += 1;
        }
        match self.mode {
            Mode::Intercept => {
                for evt in self
                    .intercept
                    .events
                    .iter()
                    .filter(|x| x.time_s > self.displayed_latest_time_s)
                {
                    self.effects.push(Effect {
                        pos_m: evt.pos_m,
                        start_s: evt.time_s,
                        magnitude: evt.energy_j.0 / 1e17,
                        tone: if evt.target_is_player { RED } else { SKYBLUE },
                    });
                }
                if let Some(last) = self.intercept.events.last() {
                    self.displayed_latest_time_s = last.time_s;
                }
                let camera_target = (self.intercept.player.pos_m + self.intercept.opponent.pos_m) * 0.5;
                self.camera.center_m += (camera_target - self.camera.center_m) * 0.08;
            }
            Mode::Tidal => {
                self.camera.center_m += (self.tidal.black_hole.pos_m - self.camera.center_m) * 0.06;
            }
            Mode::Boxing => {
                let target = (self.boxing.player.pos + self.boxing.opponent.pos) * 0.5;
                self.camera.center_m += (target - self.camera.center_m) * 0.12;
            }
        }
        let now = self.current_time();
        self.effects.retain(|e| now - e.start_s < 2.4);
    }
    fn stars(&self) {
        for i in 0..370_u32 {
            let x = ((i as i64 * 2207 + 3181) % 8000) as f32 / 8000.0 * screen_width();
            let y = ((i as i64 * 3911 + 701) % 8000) as f32 / 8000.0 * screen_height();
            let sz = if i % 31 == 0 { 1.7 } else { 0.8 };
            draw_circle(
                x,
                y,
                sz,
                Color::new(0.45, 0.58, 0.79, 0.2 + (i % 5) as f32 * 0.06),
            );
        }
    }
    fn glow(&self, p: Vec2, r: f32, c: Color) {
        for layer in (1..=5).rev() {
            let k = layer as f32;
            draw_circle(
                p.x,
                p.y,
                r * (1.0 + k * 0.48),
                Color::new(c.r, c.g, c.b, 0.032),
            );
        }
    }
    fn actor(&self, pos_m: V2, theta: f64, c: Color, size: f32) {
        let p = self.screen(pos_m);
        let h = vec2(theta.cos() as f32, theta.sin() as f32);
        let w = vec2(-h.y, h.x);
        self.glow(p, size, c);
        draw_triangle(
            p + h * size * 1.4,
            p - h * size + w * size * 0.80,
            p - h * size - w * size * 0.80,
            c,
        );
        draw_circle_lines(p.x, p.y, size * 1.5, 1.0, c);
    }
    fn labeled(&self, title: &str, value: &str, y: f32, c: Color) {
        draw_text(title, 24.0, y, 20.0, LIGHTGRAY);
        draw_text(value, 245.0, y, 20.0, c);
    }
    fn paint_intercept(&self) {
        let player = &self.intercept.player;
        let foe = &self.intercept.opponent;
        let a = self.screen(player.pos_m);
        let b = self.screen(foe.pos_m);
        let shield_player =
            (player.cross_section_radius_m / self.camera.metres_per_pixel.0) as f32;
        let shield_foe = (foe.cross_section_radius_m / self.camera.metres_per_pixel.0) as f32;
        draw_circle_lines(a.x, a.y, shield_player.max(9.0), 1.4, SKYBLUE);
        draw_circle_lines(b.x, b.y, shield_foe.max(9.0), 1.4, ORANGE);
        self.actor(player.pos_m, player.facing_rad, SKYBLUE, 16.0);
        self.actor(foe.pos_m, foe.facing_rad, ORANGE, 19.0);
        for escort in &self.intercept.escorts {
            let craft = &escort.craft;
            let pos = self.screen(craft.pos_m);
            if craft.alive() {
                draw_circle_lines(
                    pos.x,
                    pos.y,
                    (craft.cross_section_radius_m / self.camera.metres_per_pixel.0) as f32,
                    1.2,
                    ORANGE,
                );
                self.actor(craft.pos_m, craft.facing_rad, ORANGE, 13.0);
            } else {
                draw_circle_lines(pos.x, pos.y, 8.0, 1.0, DARKGRAY);
            }
        }
        if let Some((apparent, delay)) = self.intercept.apparent_opponent() {
            let ghost = self.screen(apparent);
            draw_circle_lines(ghost.x, ghost.y, 12.0, 1.8, GOLD);
            draw_line(
                ghost.x,
                ghost.y,
                b.x,
                b.y,
                1.2,
                Color::new(1.0, 0.86, 0.38, 0.54),
            );
            draw_text(
                "LIGHT-CONE OBSERVATION",
                ghost.x + 17.0,
                ghost.y + 5.0,
                15.0,
                GOLD,
            );
            draw_text(
                &format!("T - {:.2} s", delay),
                ghost.x + 17.0,
                ghost.y + 24.0,
                15.0,
                YELLOW,
            );
        }
        for round in &self.intercept.rounds {
            let p = self.screen(round.pos_m);
            let tail = self.screen(round.pos_m - round.vel_m_s * 0.055);
            let tone = if round.is_player { SKYBLUE } else { RED };
            draw_line(tail.x, tail.y, p.x, p.y, 2.5, tone);
            draw_circle(p.x, p.y, 2.4, WHITE);
        }
        let now = self.current_time();
        for effect in &self.effects {
            let p = self.screen(effect.pos_m);
            let age = ((now - effect.start_s) / 2.4).clamp(0.0, 1.0) as f32;
            let radius = (30.0 + effect.magnitude.sqrt() as f32 * 34.0) * (0.3 + age * 2.2);
            draw_circle_lines(
                p.x,
                p.y,
                radius,
                2.3,
                Color::new(effect.tone.r, effect.tone.g, effect.tone.b, 1.0 - age),
            );
        }
        let px = self.camera.metres_per_pixel.0;
        let distance = (foe.pos_m - player.pos_m).len();
        let light_delay = distance / C_SI;
        let last = self.intercept.last_impact();
        self.labeled("RANGE", &format!("{:.0} Mm", distance / 1e6), 315.0, WHITE);
        self.labeled("LIGHT TIME", &format!("{:.3} s", light_delay), 344.0, GOLD);
        self.labeled(
            "SHIP SPEED",
            &format!("{:.5} c", player.vel_m_s.len() / C_SI),
            373.0,
            SKYBLUE,
        );
        self.labeled(
            "LORENTZ",
            &format!("gamma = {:.5}", player.gamma()),
            402.0,
            SKYBLUE,
        );
        self.labeled(
            "PROPER TIME",
            &format!("{:.3} s", player.proper_time_s),
            431.0,
            WHITE,
        );
        self.labeled(
            "PLAYER SHIELD",
            &format!("{:.1}%", 100.0 * player.shield_energy_j / player.shield_capacity_j),
            460.0,
            SKYBLUE,
        );
        self.labeled(
            "ENEMY SHIELD",
            &format!("{:.1}%", 100.0 * foe.shield_energy_j / foe.shield_capacity_j),
            489.0,
            ORANGE,
        );
        self.labeled("RKM SHOTS", &self.intercept.shots_fired.to_string(), 518.0, LIGHTGRAY);
        if let Some(evt) = last {
            self.labeled(
                "IMPACT ENERGY",
                &format!("{:.3e} J", evt.energy_j.0),
                547.0,
                YELLOW,
            );
        }
        self.labeled("SCALE", &format!("{:.3} km/px", px / 1000.0), 576.0, LIGHTGRAY);
        if !self.intercept.escorts.is_empty() {
            let active = self
                .intercept
                .escorts
                .iter()
                .filter(|e| e.craft.alive())
                .count()
                + usize::from(foe.alive());
            self.labeled(
                "HOSTILE FLEET",
                &format!("{} of {} operational", active, 1 + self.intercept.escorts.len()),
                605.0,
                ORANGE,
            );
        }
        if can_enter_contact_boxing(&self.intercept) {
            draw_text(
                "VELOCITY MATCHED — PRESS 7 FOR CLOSE-CONTACT BOXING",
                screen_width() * 0.24,
                screen_height() * 0.74,
                26.0,
                GOLD,
            );
        }
        if !foe.alive() && self.intercept.escorts.iter().all(|e| !e.craft.alive()) {
            draw_text(
                "TARGET SHIELD DISABLED",
                screen_width() * 0.32,
                screen_height() * 0.8,
                39.0,
                SKYBLUE,
            );
        }
        if !player.alive() {
            draw_text(
                "PLAYER SHIELD DISABLED",
                screen_width() * 0.32,
                screen_height() * 0.8,
                39.0,
                RED,
            );
        }
        self.labeled("WEAPON", self.intercept.player_weapon.name(), 634.0, SKYBLUE);
        self.labeled(
            "AMMO",
            &if self.mission_profile.is_some() {
                format!("{} units", self.intercept.ammo_units)
            } else {
                "SANDBOX / UNLIMITED".into()
            },
            663.0,
            GOLD,
        );
        draw_text(
            "SHIELD RADII ARE FICTIONAL DEFENSE ENVELOPES; SHIP MARKERS EXAGGERATED",
            25.0,
            screen_height() - 17.0,
            17.0,
            GRAY,
        );
    }
    fn paint_tidal(&self) {
        let b = &self.tidal.black_hole;
        let star = &self.tidal.white_dwarf;
        let bh = self.screen(b.pos_m);
        self.glow(bh, 19.0, ORANGE);
        draw_circle(bh.x, bh.y, 21.0, ORANGE);
        draw_circle(bh.x, bh.y, 16.0, BLACK);
        draw_circle_lines(bh.x, bh.y, 24.0, 2.0, YELLOW);
        draw_text("10 M_SUN BLACK HOLE", bh.x + 32.0, bh.y + 7.0, 18.0, ORANGE);
        let white = self.screen(star.pos_m);
        if self.tidal.disruption.is_none() {
            self.glow(white, 13.0, SKYBLUE);
            draw_circle(white.x, white.y, 14.0, LIGHTGRAY);
            draw_text("0.65 M_SUN WHITE DWARF", white.x + 19.0, white.y + 7.0, 17.0, WHITE);
        } else {
            let event = self.tidal.disruption.unwrap();
            let locus = self.screen(event.position_m);
            draw_circle_lines(
                locus.x,
                locus.y,
                (self.tidal.time_s - event.time_s).max(0.0) as f32 * 18.0 + 15.0,
                1.5,
                Color::new(1.0, 0.45, 0.1, 0.6),
            );
            for fragment in &self.tidal.fragments {
                let p = self.screen(fragment.pos_m);
                if p.x > 0.0 && p.y > 0.0 && p.x < screen_width() && p.y < screen_height() {
                    draw_circle(p.x, p.y, 2.3, Color::new(0.56, 0.82, 1.0, 0.86));
                }
            }
        }
        let dist = (b.pos_m - star.pos_m).len();
        let q = blighted_galaxy::physics::near_side_tidal_ratio_si(b, star).unwrap_or(f64::NAN);
        self.labeled("SEPARATION", &format!("{:.2} Mm", dist / 1e6), 315.0, WHITE);
        self.labeled(
            "ORBIT SPEED",
            &format!("{:.1} km/s", star.vel_m_s.len() / 1000.0),
            344.0,
            WHITE,
        );
        self.labeled("TIDAL / SELF G", &format!("{:.4}", q), 373.0, YELLOW);
        self.labeled(
            "MAX TIDAL RATIO",
            &format!("{:.4}", self.tidal.peak_tidal_ratio),
            402.0,
            YELLOW,
        );
        self.labeled("ELAPSED", &format!("{:.2} s", self.tidal.time_s), 431.0, LIGHTGRAY);
        self.labeled(
            "MASS TRACED",
            &format!(
                "{:.4} M_SUN",
                self.tidal.total_visible_mass_kg().0 / SOLAR_MASS_KG
            ),
            460.0,
            SKYBLUE,
        );
        if let Some(event) = self.tidal.disruption {
            self.labeled(
                "DISRUPTION",
                &format!("at {:.2} s; q = {:.3}", event.time_s, event.ratio),
                489.0,
                ORANGE,
            );
            self.labeled(
                "DEBRIS",
                &format!("{} tracers", self.tidal.fragments.len()),
                518.0,
                WHITE,
            );
            draw_text(
                "PROXY BREAKUP: DEBRIS IS BALLISTIC, WITHOUT STELLAR HYDRODYNAMICS",
                20.0,
                screen_height() - 17.0,
                16.0,
                GRAY,
            );
        } else {
            draw_text(
                "NEAR-SIDE TIDAL ACCELERATION / SELF-GRAVITY; FAILURE GATE AT 1",
                20.0,
                screen_height() - 17.0,
                17.0,
                GRAY,
            );
        }
    }
    fn paint_morph_duelist(&self, d: &MorphDuelist, base_color: Color, label: &str) {
        let p = self.screen(d.pos);
        let stats = d.stats();
        let px_per_m = (1.0 / self.camera.metres_per_pixel.0) as f32;
        let r = (stats.body_radius_m as f32 * px_per_m).max(10.0);
        let h = vec2(d.facing_rad.cos() as f32, d.facing_rad.sin() as f32);
        let perp = vec2(-h.y, h.x);

        // Form-specific visual geometry & aura.
        match d.form {
            MorphForm::Streamlined => {
                self.glow(p, r * 0.9, base_color);
                draw_triangle(
                    p + h * r * 1.9,
                    p - h * r * 1.1 + perp * r * 0.65,
                    p - h * r * 1.1 - perp * r * 0.65,
                    base_color,
                );
            }
            MorphForm::AegisGuard => {
                self.glow(p, r * 1.1, GOLD);
                draw_poly(p.x, p.y, 6, r * 1.15, d.facing_rad.to_degrees() as f32, Color::new(base_color.r, base_color.g, base_color.b, 0.35));
                draw_poly_lines(p.x, p.y, 6, r * 1.25, d.facing_rad.to_degrees() as f32, 3.0, GOLD);
                draw_circle(p.x, p.y, r * 0.6, base_color);
            }
            MorphForm::StrikerLimb => {
                self.glow(p, r, base_color);
                let reach = (stats.strike_reach_m as f32 * px_per_m).min(110.0);
                draw_line(
                    p.x,
                    p.y,
                    p.x + h.x * reach,
                    p.y + h.y * reach,
                    4.0,
                    Color::new(1.0, 0.85, 0.35, 0.65),
                );
                draw_circle(p.x, p.y, r, base_color);
                draw_circle_lines(p.x + h.x * reach, p.y + h.y * reach, 6.0, 2.0, YELLOW);
            }
            MorphForm::SwarmErode => {
                let cloud_r = stats.erode_radius_m as f32 * px_per_m;
                draw_circle(p.x, p.y, cloud_r, Color::new(0.68, 0.25, 0.95, 0.14));
                draw_circle_lines(p.x, p.y, cloud_r, 1.5, Color::new(0.75, 0.35, 1.0, 0.65));
                draw_circle(p.x, p.y, r, base_color);
            }
            MorphForm::SingularityCore => {
                self.glow(p, r * 1.4, VIOLET);
                draw_circle_lines(p.x, p.y, r * 1.5, 2.2, GOLD);
                draw_circle(p.x, p.y, r * 0.85, BLACK);
                draw_circle_lines(p.x, p.y, r * 0.9, 2.0, VIOLET);
            }
        }

        draw_circle_lines(p.x, p.y, r + 3.0, 1.4, base_color);
        draw_text(label, p.x - 28.0, p.y - r - 18.0, 15.0, base_color);
        draw_text(d.form.label(), p.x - 46.0, p.y - r - 4.0, 13.0, GOLD);
    }
    fn paint_boxing(&self) {
        // Draw arena boundary & spatial partition grid lines.
        let top_left = self.screen(V2::new(-ARENA_HALF_EXTENT_M, -ARENA_HALF_EXTENT_M));
        let bot_right = self.screen(V2::new(ARENA_HALF_EXTENT_M, ARENA_HALF_EXTENT_M));
        draw_rectangle_lines(
            top_left.x,
            top_left.y,
            bot_right.x - top_left.x,
            bot_right.y - top_left.y,
            1.5,
            Color::new(0.35, 0.62, 0.85, 0.55),
        );
        let cell = self.boxing.grid.cell_size;
        let mut g = -ARENA_HALF_EXTENT_M;
        while g <= ARENA_HALF_EXTENT_M {
            let a = self.screen(V2::new(g, -ARENA_HALF_EXTENT_M));
            let b = self.screen(V2::new(g, ARENA_HALF_EXTENT_M));
            draw_line(a.x, a.y, b.x, b.y, 0.7, Color::new(0.25, 0.42, 0.62, 0.16));
            let c = self.screen(V2::new(-ARENA_HALF_EXTENT_M, g));
            let d = self.screen(V2::new(ARENA_HALF_EXTENT_M, g));
            draw_line(c.x, c.y, d.x, d.y, 0.7, Color::new(0.25, 0.42, 0.62, 0.16));
            g += cell;
        }

        // Memetic pulse waves.
        let px_per_m = (1.0 / self.camera.metres_per_pixel.0) as f32;
        for pulse in &self.boxing.pulses {
            let p = self.screen(pulse.origin);
            let alpha = (1.0 - pulse.age_s / 0.45).clamp(0.0, 1.0) as f32;
            let tone = if pulse.from_player {
                Color::new(0.45, 0.95, 1.0, alpha)
            } else {
                Color::new(1.0, 0.35, 0.85, alpha)
            };
            draw_circle_lines(p.x, p.y, pulse.radius_m as f32 * px_per_m, 2.5, tone);
        }

        // Swept continuous strike volumes.
        for strike in &self.boxing.strikes {
            let a = self.screen(strike.start);
            let b = self.screen(strike.end);
            let alpha = (1.0 - strike.age_s / 0.28).clamp(0.0, 1.0) as f32;
            let tone = if strike.hit {
                Color::new(1.0, 0.88, 0.25, alpha)
            } else if strike.from_player {
                Color::new(0.35, 0.85, 1.0, alpha)
            } else {
                Color::new(1.0, 0.42, 0.25, alpha)
            };
            draw_line(a.x, a.y, b.x, b.y, 5.0, tone);
            draw_circle_lines(b.x, b.y, strike.radius_m as f32 * px_per_m, 2.0, tone);
        }

        // Contact sparks.
        for spark in &self.boxing.sparks {
            let p = self.screen(spark.pos);
            let alpha = (1.0 - spark.age_s / 0.35).clamp(0.0, 1.0) as f32;
            let r = (10.0 + spark.age_s as f32 * 70.0) * spark.intensity as f32;
            draw_circle_lines(p.x, p.y, r, 2.2, Color::new(1.0, 0.92, 0.45, alpha));
        }

        self.paint_morph_duelist(&self.boxing.player, SKYBLUE, "ALFONSO / PLAYER");
        self.paint_morph_duelist(&self.boxing.opponent, ORANGE, "BLIGHT AVATAR");

        // Resource & morphogenesis HUD.
        let p = &self.boxing.player;
        let o = &self.boxing.opponent;
        self.labeled("PLAYER FORM", p.form.label(), 315.0, SKYBLUE);
        self.labeled(
            "VIABLE MASS",
            &format!("{:.0} / {:.0} kg", p.mass_kg, p.max_mass_kg),
            344.0,
            SKYBLUE,
        );
        self.labeled(
            "ENERGY / WATT",
            &format!("{:.0} MJ / {:.0} MW", p.energy_j / 1e6, p.wattage_w / 1e6),
            373.0,
            GOLD,
        );
        self.labeled(
            "MATTER STOCK",
            &format!("{:.0} / {:.0} kg", p.matter_kg, p.max_matter_kg),
            402.0,
            WHITE,
        );
        self.labeled(
            "COHERENCE / HEAT",
            &format!("{:.0}% / {:.0}%", p.coherence * 100.0, p.heat * 100.0),
            431.0,
            if p.coherence < 0.35 { RED } else { SKYBLUE },
        );

        self.labeled("OPPONENT FORM", o.form.label(), 474.0, ORANGE);
        self.labeled(
            "OPP MASS / COH",
            &format!("{:.0} kg / {:.0}%", o.mass_kg, o.coherence * 100.0),
            503.0,
            ORANGE,
        );
        self.labeled(
            "CONTACTS / MEMES",
            &format!(
                "{} swept / {} pulses",
                self.boxing.contacts_resolved, self.boxing.memetic_pulses_landed
            ),
            532.0,
            GOLD,
        );

        match self.boxing.outcome() {
            TacticalOutcome::Victory => {
                draw_text(
                    "OPPONENT COHERENCE BROKEN — TAB TO RECORD TRIE CONTACT",
                    screen_width() * 0.22,
                    screen_height() * 0.82,
                    30.0,
                    SKYBLUE,
                );
            }
            TacticalOutcome::Defeat => {
                draw_text(
                    "CHASSIS MASS / COHERENCE LOST — TAB TO WITHDRAW",
                    screen_width() * 0.24,
                    screen_height() * 0.82,
                    30.0,
                    RED,
                );
            }
            TacticalOutcome::Withdrawn => {}
        }

        draw_text(
            "M3 MORPHOGENESIS: UNIFORM-GRID BROAD-PHASE + SWEPT CIRCLE/AABB COLLISION + ENERGY-LIMITED REGEN",
            20.0,
            screen_height() - 17.0,
            16.0,
            GRAY,
        );
    }
    fn draw(&self) {
        clear_background(Color::new(0.010, 0.015, 0.042, 1.0));
        self.stars();
        match self.mode {
            Mode::Intercept => self.paint_intercept(),
            Mode::Tidal => self.paint_tidal(),
            Mode::Boxing => self.paint_boxing(),
        }
        draw_rectangle(13.0, 12.0, 825.0, 207.0, Color::new(0.022, 0.033, 0.058, 0.92));
        draw_rectangle_lines(13.0, 12.0, 825.0, 207.0, 1.0, DARKGRAY);
        draw_text("BLIGHTED GALAXY / SI GAMEPLAY BRIDGE", 23.0, 45.0, 27.0, WHITE);
        let title = match self.mode {
            Mode::Intercept
                if self
                    .mission_profile
                    .is_some_and(|p| p.kind == TacticalKind::GufDefense) =>
            {
                "4  GUF RELAY INTERCEPTION / DEFEND A THREATENED WORLD"
            }
            Mode::Intercept
                if self
                    .mission_profile
                    .is_some_and(|p| p.kind == TacticalKind::ConvoyEscort) =>
            {
                "5  CONVOY ESCORT / REPEL FAST BANDIT RAIDERS"
            }
            Mode::Intercept
                if self
                    .mission_profile
                    .is_some_and(|p| p.kind == TacticalKind::HaloStrike) =>
            {
                "6  HALO STRIKE / 0.76c SIEGE SHIELD BREACH"
            }
            Mode::Intercept => "1  RELATIVISTIC FLEET INTERCEPT • WEAPON LOADOUTS",
            Mode::Tidal => "2  WHITE-DWARF TIDAL INTERDICTION • STELLAR REMNANT",
            Mode::Boxing => "7  M3 TRANSHUMAN CLOSE-CONTACT BOXING • MORPHOGENESIS",
        };
        draw_text(title, 23.0, 80.0, 22.0, SKYBLUE);
        if self.mode == Mode::Boxing {
            draw_text(
                "1-5 MORPH FORM   WASD MOVE   MOUSE AIM   SPACE STRIKE   F MEMETIC PULSE",
                23.0,
                113.0,
                18.0,
                LIGHTGRAY,
            );
            draw_text(
                "1 STREAM  2 AEGIS  3 STRIKER  4 SWARM  5 SINGULARITY  |  TAB DEBRIEF",
                23.0,
                145.0,
                18.0,
                GOLD,
            );
        } else {
            draw_text(
                "1 / 2 SCENARIOS   3 / TAB GALAXY   Q / E TIME   Z / X / WHEEL ZOOM",
                23.0,
                113.0,
                19.0,
                LIGHTGRAY,
            );
            if self.mode == Mode::Intercept {
                draw_text(
                    "MOUSE aim / M keyboard aim / 4 RKM  5 SIEGE  6 FLAK / SPACE fire / 7 BOXING",
                    23.0,
                    145.0,
                    17.5,
                    LIGHTGRAY,
                );
            } else if self.mission_profile.is_some() {
                draw_text(
                    &format!(
                        "A inward tug / D outward tug  |  burns {}  |  300s objective",
                        self.tidal_tug_budget
                    ),
                    23.0,
                    145.0,
                    18.0,
                    LIGHTGRAY,
                );
            } else {
                draw_text(
                    "Observe calibrated freefall and physically scaled tidal ratio",
                    23.0,
                    145.0,
                    18.0,
                    LIGHTGRAY,
                );
            }
        }
        draw_text(
            &format!(
                "SIM RATE {:5.3}x | COORDINATE TIME {:8.3} s | R RESET | ESC QUIT",
                self.sim_rate,
                self.current_time()
            ),
            23.0,
            179.0,
            18.0,
            GOLD,
        );
        draw_text(
            "SI kinematics + Newtonian gravity • exotic shields, morphogenesis and breakup are game rules",
            23.0,
            204.0,
            15.5,
            GRAY,
        );
        if let Some(site) = self.mission_name {
            let status = match self.tactical_outcome() {
                TacticalOutcome::Victory => "OBJECTIVE COMPLETE",
                TacticalOutcome::Defeat => "SHIP / CHASSIS LOST",
                TacticalOutcome::Withdrawn => "OBJECTIVE INCOMPLETE",
            };
            draw_rectangle(13.0, 221.0, 825.0, 75.0, Color::new(0.02, 0.03, 0.06, 0.88));
            draw_text(
                &format!("CAMPAIGN: {site} / {status}"),
                23.0,
                244.0,
                20.0,
                GOLD,
            );
            if let Some(profile) = self.mission_profile {
                let conditions = if self.mode == Mode::Intercept {
                    if profile.kind == TacticalKind::GufDefense {
                        format!(
                            "GUF THREAT {} beams / {:.1} attack units / repel spotter fleet",
                            profile.incoming_guf_count, profile.incoming_guf_power
                        )
                    } else if profile.kind == TacticalKind::ConvoyEscort {
                        format!(
                            "CONVOY #{} ESCORT / {} fast bandit raiders @ {:.3}c",
                            profile.convoy_id.unwrap_or(0),
                            profile.escort_count + 1,
                            profile.enemy_initial_beta
                        )
                    } else if profile.kind == TacticalKind::HaloStrike {
                        format!(
                            "HALO STRIKE APPROACH {:.2}c / SIEGE ARMED / OPP SHIELD {:.2e} J",
                            profile.player_initial_beta, profile.enemy_shield_j
                        )
                    } else {
                        format!(
                            "IND {:.1}  DEF {:.1}  HALO {}   OPP SHIELD {:.2e} J   FIRE {:.2}s",
                            profile.industry,
                            profile.defense,
                            profile.kugelblitz,
                            profile.enemy_shield_j,
                            profile.enemy_fire_interval_s
                        )
                    }
                } else {
                    format!(
                        "IND {:.1} DEF {:.1}  HOLE {:.1} M_sun  DWARF {:.2} M_sun",
                        profile.industry,
                        profile.defense,
                        profile.hole_mass_solar,
                        profile.dwarf_mass_solar
                    )
                };
                draw_text(&conditions, 23.0, 267.0, 16.0, SKYBLUE);
                if self.mode == Mode::Intercept {
                    draw_text(
                        &format!(
                            "FLEET: 1 flagship + {} escorts  |  total {} ships",
                            profile.escort_count,
                            profile.escort_count + 1
                        ),
                        850.0,
                        76.0,
                        17.0,
                        ORANGE,
                    );
                    draw_text(
                        &format!(
                            "ARMORY: 4 RKM / 5 SIEGE {} / 6 FLAK {}",
                            if profile.siege_unlocked { "READY" } else { "LOCKED" },
                            if profile.fan_unlocked { "READY" } else { "LOCKED" }
                        ),
                        850.0,
                        104.0,
                        16.0,
                        SKYBLUE,
                    );
                }
            } else if self.mode == Mode::Boxing {
                draw_text(
                    "CLOSE-CONTACT MORPH DUEL / MEMETIC PULSES FEED PATRICIA TRIE ON DEBRIEF",
                    23.0,
                    267.0,
                    15.5,
                    SKYBLUE,
                );
            }
            draw_text(
                if self.demo_preview {
                    "QA PREVIEW ONLY / TAB RETURNS TO GALAXY WITHOUT REWARDS"
                } else if self.mode == Mode::Boxing {
                    "TAB: RETURN TO GALAXY AND RECORD CONTACT DUEL"
                } else {
                    "3 / TAB: RETURN TO GALAXY AND RECORD MISSION RESULT"
                },
                23.0,
                286.0,
                15.0,
                LIGHTGRAY,
            );
        }
    }
}

#[macroquad::main(config)]
async fn main() {
    let mut game = Game::new(Mode::Intercept);
    let mut overworld = OverworldNav::new();
    let args = std::env::args().collect::<Vec<_>>();
    let save_path = args
        .windows(2)
        .find(|pair| pair[0] == "--save-file")
        .map(|pair| PathBuf::from(&pair[1]))
        .unwrap_or_else(|| PathBuf::from("blighted-galaxy-save.json"));
    // Offline true-game-frame render for CI and art review; the same
    // procedural GPU scene draws during ordinary gameplay.
    let snapshot_path = args
        .windows(2)
        .find(|pair| pair[0] == "--snapshot-overworld")
        .map(|pair| pair[1].clone());
    if args.iter().any(|v| v == "--focus-orion") {
        overworld.focus_frontier();
    }
    // The campaign opens on its galactic navigation layer. Existing
    // physics encounters remain directly reachable with keys 1 and 2.
    // --encounter keeps the previous startup for dedicated combat tests.
    let encounter_start = args.iter().any(|a| a == "--encounter");
    let mut galaxy_mode = !encounter_start;
    let mut active_mission: Option<(usize, TacticalKind)> = None;

    // Deterministic real-engine QA snapshots of named industrial & M3 boxing missions.
    // Preview scenes deliberately have no strategic reward linkage.
    let demo_boxing = args.iter().any(|arg| arg == "--demo-boxing");
    let demo_escort = args.iter().any(|arg| arg == "--demo-convoy-escort");
    let demo_halo = args.iter().any(|arg| arg == "--demo-halo-strike");

    if demo_boxing {
        game = Game::for_boxing_demo();
        galaxy_mode = false;
    } else if let Some(raw_site) = args
        .windows(2)
        .find(|pair| pair[0] == "--demo-mission")
        .map(|pair| &pair[1])
        .cloned()
        .or_else(|| {
            if demo_escort {
                Some("3".into())
            } else if demo_halo {
                Some("6".into())
            } else {
                None
            }
        })
    {
        if let Ok(site) = raw_site.parse::<usize>() {
            let kind = if args.iter().any(|arg| arg == "--demo-tidal") {
                TacticalKind::Tidal
            } else if args.iter().any(|arg| arg == "--demo-guf-defense") {
                for i in [2, 3, 6] {
                    overworld.strategy.set_owner(i, Faction::Weena);
                    overworld.strategy.systems[i].energy = 500.0;
                    overworld.strategy.systems[i].array_node = true;
                }
                overworld.strategy.blight_autonomy = false;
                if let Ok(ready) = overworld.strategy.aim_guf(site) {
                    overworld.strategy.advance_to(ready);
                    overworld.strategy.array_charged = 100.0;
                    let _ = overworld.strategy.fire_guf(site);
                }
                TacticalKind::GufDefense
            } else if demo_escort {
                overworld.strategy.blight_autonomy = false;
                overworld.strategy.systems[3].material = 350.0;
                overworld.strategy.systems[3].energy = 350.0;
                let dest = if site == 3 { 5 } else { 3 };
                let src = if site < overworld.strategy.systems.len()
                    && overworld.strategy.systems[site].owner == Faction::Drummers
                {
                    site
                } else {
                    3
                };
                overworld.strategy.systems[src].material = 350.0;
                overworld.strategy.systems[src].energy = 350.0;
                let _ = overworld.strategy.send_convoy(src, dest);
                TacticalKind::ConvoyEscort
            } else if demo_halo {
                overworld.strategy.blight_autonomy = false;
                overworld.strategy.systems[5].kugelblitz = 2;
                if site < overworld.strategy.systems.len()
                    && !matches!(
                        overworld.strategy.systems[site].owner,
                        Faction::Blight | Faction::Weena
                    )
                {
                    overworld.strategy.set_owner(site, Faction::Weena);
                }
                TacticalKind::HaloStrike
            } else {
                TacticalKind::Intercept
            };
            if let Some(profile) = MissionProfile::from_strategy(&overworld.strategy, site, kind) {
                game = Game::for_mission(profile);
                game.demo_preview = true;
                galaxy_mode = false;
            }
        }
    }

    // Allow the built-in font atlas and GPU textures to settle before
    // capturing an actual frame in headless CI.
    let mut frame_count = 0usize;
    loop {
        if is_key_pressed(KeyCode::Escape) {
            break;
        }
        let dt = get_frame_time() as f64;
        if galaxy_mode {
            if is_key_pressed(KeyCode::F5) {
                let message =
                    match campaign_save::save_file(&save_path, &overworld.model, &overworld.strategy)
                    {
                        Ok(()) => format!("CAMPAIGN SAVED / {}", save_path.display()),
                        Err(e) => format!("SAVE FAILED: {e}"),
                    };
                overworld.report_notice(&message);
            }
            if is_key_pressed(KeyCode::F9) {
                let message = match campaign_save::load_file(&save_path) {
                    Ok((nav, strategy)) => {
                        overworld.import_campaign(nav, strategy);
                        format!("CAMPAIGN RESTORED / {}", save_path.display())
                    }
                    Err(e) => format!("LOAD FAILED: {e}"),
                };
                overworld.report_notice(&message);
                active_mission = None;
            }
            let intercept = is_key_pressed(KeyCode::Key1);
            let tidal = is_key_pressed(KeyCode::Key2);
            let guf = is_key_pressed(KeyCode::Key4);
            let escort = is_key_pressed(KeyCode::Key5);
            let halo = is_key_pressed(KeyCode::Key6);
            let boxing = is_key_pressed(KeyCode::Key7);

            if boxing {
                let site = overworld
                    .selected_tactical_site()
                    .or_else(|| overworld.current_ship_site());
                let name = site.map(|i| overworld.strategy.systems[i].name);
                game = Game::for_boxing(site, name);
                active_mission = None;
                galaxy_mode = false;
            } else if intercept || tidal || guf || escort || halo {
                let kind = if guf {
                    TacticalKind::GufDefense
                } else if escort {
                    TacticalKind::ConvoyEscort
                } else if halo {
                    TacticalKind::HaloStrike
                } else if tidal {
                    TacticalKind::Tidal
                } else {
                    TacticalKind::Intercept
                };
                let site_opt = if escort {
                    overworld.selected_convoy_escort_site()
                } else {
                    overworld.selected_tactical_site()
                };
                let profile = site_opt
                    .and_then(|site| MissionProfile::from_strategy(&overworld.strategy, site, kind));
                if let Some(profile) = profile {
                    game = Game::for_mission(profile);
                    active_mission = Some((profile.site, kind));
                    galaxy_mode = false;
                } else if guf {
                    overworld.report_notice(
                        "NO GUF DEFENSE MISSION: ARRIVE AT A THREATENED DRUMMER SYSTEM",
                    );
                } else if escort {
                    overworld.report_notice(
                        "NO ESCORTABLE CONVOY: DISPATCH A DRUMMER CONVOY (M) NEAR YOUR SHIP",
                    );
                } else if halo {
                    overworld.report_notice(
                        "NO HALO STRIKE: ARRIVE AT HOSTILE SYSTEM WITH STORED DRUMMER KUGELBLITZ",
                    );
                } else {
                    game = Game::new(if tidal { Mode::Tidal } else { Mode::Intercept });
                    active_mission = None;
                    galaxy_mode = false;
                }
            } else {
                overworld.update(dt);
            }
        } else {
            let return_key = if game.mode == Mode::Boxing {
                is_key_pressed(KeyCode::Tab)
            } else {
                is_key_pressed(KeyCode::Key3) || is_key_pressed(KeyCode::Tab)
            };
            if return_key {
                if !game.demo_preview {
                    if game.mode == Mode::Boxing {
                        overworld.boxing_result(game.boxing.debrief());
                        active_mission = None;
                    } else if let Some((site, kind)) = active_mission.take() {
                        overworld.tactical_result(site, kind, game.tactical_outcome());
                    }
                } else {
                    active_mission = None;
                }
                galaxy_mode = true;
            } else {
                game.input();
                game.advance(dt);
            }
        }
        if galaxy_mode {
            overworld.draw();
        } else {
            game.draw();
        }
        frame_count += 1;
        if let Some(path) = snapshot_path.as_ref() {
            if frame_count >= 4 {
                // Real GPU screenshot after font glyphs were loaded.
                macroquad::prelude::get_screen_data().export_png(path);
                break;
            }
        }
        next_frame().await;
    }
}
