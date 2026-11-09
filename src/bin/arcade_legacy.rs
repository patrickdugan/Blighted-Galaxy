use blighted_galaxy::{
    acceleration_at, scatter_swarm, step_bodies, step_swarm, swept_circle_hit, tidal_ratio,
    Body, BodyKind, Duelist, HawkingBomb, SwarmNode, V2,
};
use macroquad::prelude::*;
use std::f64::consts::{PI, TAU};

const STEP: f64 = 1.0 / 120.0;
const MAX_SPARKS: usize = 8_000;

fn window_conf() -> Conf {
    Conf {
        window_title: "Blighted Galaxy: recovered physics sandbox".to_string(),
        window_width: 1600,
        window_height: 900,
        high_dpi: true,
        ..Default::default()
    }
}

#[derive(Clone, Copy)]
struct Bolt {
    pos: V2,
    vel: V2,
    life: f64,
    friendly: bool,
}

#[derive(Clone, Copy)]
struct Spark {
    pos: V2,
    vel: V2,
    life: f32,
    original_life: f32,
    radius: f32,
    tint: Color,
}

struct Game {
    bodies: Vec<Body>,
    sheared: Vec<bool>,
    player: Duelist,
    opponent: Duelist,
    bolts: Vec<Bolt>,
    bombs: Vec<HawkingBomb>,
    swarm: Vec<SwarmNode>,
    sparks: Vec<Spark>,
    time: f64,
    accumulator: f64,
    camera: V2,
    zoom: f32,
    player_cooldown: f64,
    opponent_cooldown: f64,
    punch_cooldown: f64,
    last_event: &'static str,
}

impl Game {
    fn new() -> Self {
        let bodies = vec![
            Body::new(BodyKind::Star, V2::new(-200.0, -170.0), V2::new(0.0, 0.0), 70_000.0, 76.0),
            Body::new(BodyKind::BlackHole, V2::new(650.0, -270.0), V2::new(-2.0, 8.0), 360_000.0, 30.0),
            Body::new(BodyKind::WhiteDwarf, V2::new(870.0, -180.0), V2::new(-12.0, -25.0), 21_000.0, 24.0),
            Body::new(BodyKind::NeutronStar, V2::new(-710.0, 340.0), V2::new(3.0, -6.0), 31_000.0, 13.0),
        ];
        Self {
            sheared: vec![false; bodies.len()],
            bodies,
            player: Duelist::new(V2::new(-240.0, 160.0), -0.3),
            opponent: Duelist::new(V2::new(235.0, 105.0), PI),
            bolts: Vec::new(),
            bombs: Vec::new(),
            swarm: scatter_swarm(V2::new(385.0, 40.0), 1300, 0xB11947, 195.0),
            sparks: Vec::new(),
            time: 0.0,
            accumulator: 0.0,
            camera: V2::new(-180.0, 80.0),
            zoom: 0.85,
            player_cooldown: 0.0,
            opponent_cooldown: 1.0,
            punch_cooldown: 0.0,
            last_event: "Swarm front detected",
        }
    }

    fn screen(&self, world: V2) -> Vec2 {
        vec2(
            (world.x - self.camera.x) as f32 * self.zoom + screen_width() * 0.5,
            (world.y - self.camera.y) as f32 * self.zoom + screen_height() * 0.5,
        )
    }

    fn burst(&mut self, center: V2, count: usize, tint: Color, max_speed: f32) {
        for _ in 0..count {
            let angle: f32 = macroquad::rand::gen_range(0.0, std::f32::consts::TAU);
            let speed: f32 = macroquad::rand::gen_range(0.04, 1.0) * max_speed;
            let life: f32 = macroquad::rand::gen_range(0.6, 3.2);
            let vel = V2::new((angle.cos() * speed) as f64, (angle.sin() * speed) as f64);
            self.sparks.push(Spark {
                pos: center, vel, life, original_life: life,
                radius: macroquad::rand::gen_range(0.6, 3.0), tint,
            });
        }
        if self.sparks.len() > MAX_SPARKS {
            let excess = self.sparks.len() - MAX_SPARKS;
            self.sparks.drain(0..excess);
        }
    }

    fn fire(&mut self, friendly: bool) {
        let pilot = if friendly { self.player } else { self.opponent };
        let heading = pilot.heading();
        self.bolts.push(Bolt {
            pos: pilot.pos + heading * 24.0,
            vel: pilot.vel + heading * 690.0,
            life: 2.4,
            friendly,
        });
        if friendly { self.player_cooldown = 0.14; }
        else { self.opponent_cooldown = 0.86; }
    }

    fn step(&mut self, dt: f64) {
        self.time += dt;
        self.player_cooldown = (self.player_cooldown - dt).max(0.0);
        self.opponent_cooldown = (self.opponent_cooldown - dt).max(0.0);
        self.punch_cooldown = (self.punch_cooldown - dt).max(0.0);

        let player_before = self.player.pos;
        let opponent_before = self.opponent.pos;

        if is_key_down(KeyCode::A) { self.player.angle -= 2.8 * dt; }
        if is_key_down(KeyCode::D) { self.player.angle += 2.8 * dt; }
        let thrust = if is_key_down(KeyCode::W) { 155.0 }
            else if is_key_down(KeyCode::S) { -85.0 } else { 0.0 };

        step_bodies(&mut self.bodies, dt);
        let player_g = acceleration_at(self.player.pos, &self.bodies, None);
        self.player.tick(thrust, player_g, dt);
        self.player.vel = self.player.vel * 0.9999;

        if self.opponent.hp > 0.0 {
            let to_player = self.player.pos - self.opponent.pos;
            let desired = to_player.y.atan2(to_player.x);
            let delta = (desired - self.opponent.angle + PI).rem_euclid(TAU) - PI;
            self.opponent.angle += delta.clamp(-1.7 * dt, 1.7 * dt);
            let enemy_g = acceleration_at(self.opponent.pos, &self.bodies, None);
            let enemy_thrust = if to_player.len() > 170.0 { 88.0 } else { -25.0 };
            self.opponent.tick(enemy_thrust, enemy_g, dt);
            if self.opponent_cooldown <= 0.0 && to_player.len() < 860.0 {
                self.fire(false);
            }
        }

        if is_key_down(KeyCode::Space) && self.player_cooldown <= 0.0 {
            self.fire(true);
        }
        if is_key_pressed(KeyCode::E) && self.punch_cooldown <= 0.0 {
            self.punch_cooldown = 0.55;
            let impact = self.player.pos + self.player.heading() * 34.0;
            if self.opponent.hp > 0.0 &&
                self.player.punch_hitbox().overlaps(self.opponent.box_hitbox()) {
                self.opponent.hp -= 20.0;
                self.opponent.vel += self.player.heading() * 85.0;
                self.burst(impact, 65, SKYBLUE, 130.0);
                self.last_event = "Transhuman contact strike";
            } else {
                self.burst(impact, 12, SKYBLUE, 60.0);
            }
        }
        if is_key_pressed(KeyCode::B) {
            let heading = self.player.heading();
            self.bombs.push(HawkingBomb::new(
                self.player.pos + heading * 46.0,
                self.player.vel + heading * 150.0,
                1.6,
                0.18,
            ));
            self.last_event = "Kugelblitz payload released";
        }
        if is_key_pressed(KeyCode::G) {
            self.swarm = scatter_swarm(self.player.pos + self.player.heading() * 160.0,
                1300, self.time.to_bits(), 175.0);
            self.last_event = "Nanite swarm scattered";
        }

        let mut impacts = Vec::new();
        let mut hit_player = 0.0;
        let mut hit_opponent = 0.0;
        for bolt in &mut self.bolts {
            let previous = bolt.pos;
            bolt.pos += bolt.vel * dt;
            bolt.life -= dt;
            if bolt.friendly && self.opponent.hp > 0.0 &&
                swept_circle_hit(previous, bolt.pos, 3.0,
                    opponent_before, self.opponent.pos, self.opponent.radius).is_some() {
                hit_opponent += 8.0;
                bolt.life = 0.0;
                impacts.push(bolt.pos);
            } else if !bolt.friendly &&
                swept_circle_hit(previous, bolt.pos, 3.0,
                    player_before, self.player.pos, self.player.radius).is_some() {
                hit_player += 6.0;
                bolt.life = 0.0;
                impacts.push(bolt.pos);
            }
        }
        self.bolts.retain(|bolt| bolt.life > 0.0);
        self.player.hp -= hit_player;
        self.opponent.hp -= hit_opponent;
        for pos in impacts { self.burst(pos, 18, PURPLE, 120.0); }

        let mut detonations = Vec::new();
        for bomb in &mut self.bombs {
            if bomb.step(&self.bodies, dt) { detonations.push(bomb.pos); }
        }
        self.bombs.retain(|bomb| !bomb.detonated);
        for location in detonations {
            self.burst(location, 1050, ORANGE, 540.0);
            self.last_event = "Hawking endpoint flash";
            for pilot in [&mut self.player, &mut self.opponent] {
                let offset = pilot.pos - location;
                let range = offset.len();
                if range < 280.0 {
                    pilot.hp -= (280.0 - range) * 0.16;
                    pilot.vel += offset.normalized() * 170.0;
                }
            }
            for node in &mut self.swarm {
                let delta = node.pos - location;
                if delta.len() < 400.0 {
                    node.vel += delta.normalized() * 250.0;
                }
            }
        }

        let mut shearing_events = Vec::new();
        for i in 0..self.bodies.len() {
            if self.sheared[i] || self.bodies[i].kind == BodyKind::BlackHole { continue; }
            let mut ratio: f64 = 0.0;
            for j in 0..self.bodies.len() {
                if i == j { continue; }
                ratio = ratio.max(tidal_ratio(&self.bodies[j], &self.bodies[i]));
            }
            if ratio > 1.0 {
                self.sheared[i] = true;
                self.bodies[i].radius *= 1.3;
                self.bodies[i].mass *= 0.87;
                shearing_events.push(self.bodies[i].pos);
            }
        }
        for location in shearing_events {
            self.burst(location, 1200, YELLOW, 360.0);
            self.last_event = "Stellar tidal shear";
        }

        step_swarm(&mut self.swarm, &self.bodies, dt);
        for spark in &mut self.sparks {
            spark.pos += spark.vel * dt;
            spark.vel = spark.vel * 0.992;
            spark.life -= dt as f32;
        }
        self.sparks.retain(|spark| spark.life > 0.0);
        self.camera += (self.player.pos - self.camera) * (1.0 - (-2.0 * dt).exp());
    }

    fn starfield(&self) {
        for i in 0..450_i32 {
            let x = (((i as i64 * 1801 + 1477) % 6001) - 3000) as f64;
            let y = (((i as i64 * 2539 + 871) % 5003) - 2500) as f64;
            let point = self.screen(V2::new(x, y));
            if point.x < 0.0 || point.y < 0.0 ||
                point.x > screen_width() || point.y > screen_height() { continue; }
            let a = 0.17 + (i % 7) as f32 * 0.055;
            draw_circle(point.x, point.y, if i % 19 == 0 { 1.8 } else { 0.85 },
                Color::new(0.7, 0.85, 1.0, a));
        }
    }

    fn draw_body(&self, index: usize, body: &Body) {
        let p = self.screen(body.pos);
        let r = (body.radius as f32 * self.zoom).max(2.0);
        if p.x < -300.0 || p.y < -300.0 ||
            p.x > screen_width() + 300.0 || p.y > screen_height() + 300.0 { return; }
        match body.kind {
            BodyKind::Star => {
                glow(p, r, Color::new(1.0, 0.64, 0.21, 0.9));
                draw_circle(p.x, p.y, r, Color::new(1.0, 0.87, 0.60, 1.0));
            }
            BodyKind::WhiteDwarf => {
                glow(p, r, Color::new(0.61, 0.86, 1.0, 0.8));
                draw_circle(p.x, p.y, r, WHITE);
            }
            BodyKind::NeutronStar => {
                glow(p, r * 1.2, Color::new(0.45, 0.4, 1.0, 0.8));
                draw_circle(p.x, p.y, r, Color::new(0.75, 0.72, 1.0, 1.0));
                draw_line(p.x - r * 5.0, p.y + r * 2.0,
                    p.x + r * 5.0, p.y - r * 2.0, 1.5, VIOLET);
            }
            BodyKind::BlackHole => {
                glow(p, r * 1.8, Color::new(0.86, 0.42, 0.18, 0.6));
                for k in 0..80 {
                    let angle = k as f32 * std::f32::consts::TAU / 80.0
                        + self.time as f32 * 0.23;
                    let rr = r * (1.5 + (angle * 5.0).sin() * 0.17);
                    draw_circle(p.x + angle.cos() * rr, p.y + angle.sin() * rr * 0.58,
                        2.2 * self.zoom, Color::new(1.0, 0.45, 0.13, 0.67));
                }
                draw_circle(p.x, p.y, r * 0.82, BLACK);
                draw_circle_lines(p.x, p.y, r * 1.03, 2.0, YELLOW);
            }
        }
        if self.sheared[index] {
            draw_circle_lines(p.x, p.y, r * 1.85, 2.0, ORANGE);
        }
    }

    fn draw_duelist(&self, pilot: Duelist, tint: Color) {
        let p = self.screen(pilot.pos);
        let heading = pilot.heading();
        let dir = vec2(heading.x as f32, heading.y as f32);
        let wing = vec2(-dir.y, dir.x);
        let size = 20.0 * self.zoom;
        draw_triangle(p + dir * size * 1.3,
            p - dir * size + wing * size * 0.65,
            p - dir * size - wing * size * 0.65, tint);
        draw_circle_lines(p.x, p.y, size * 1.45, 1.2, tint);
        let hp = (pilot.hp / 100.0).clamp(0.0, 1.0) as f32;
        draw_rectangle(p.x - 24.0, p.y - size - 18.0, 48.0, 4.0, DARKGRAY);
        draw_rectangle(p.x - 24.0, p.y - size - 18.0, 48.0 * hp, 4.0, tint);
    }

    fn draw(&self) {
        clear_background(Color::new(0.019, 0.027, 0.063, 1.0));
        self.starfield();

        for (index, body) in self.bodies.iter().enumerate() {
            self.draw_body(index, body);
        }
        for node in &self.swarm {
            let p = self.screen(node.pos);
            if p.x < -3.0 || p.x > screen_width() + 3.0 ||
                p.y < -3.0 || p.y > screen_height() + 3.0 { continue; }
            draw_circle(p.x, p.y, 1.3 * self.zoom.max(0.4),
                Color::new(0.70, 0.23, 1.0, 0.6));
        }

        for bomb in &self.bombs {
            let p = self.screen(bomb.pos);
            let glow_size = (18.0 + 1.3 / bomb.mass.max(0.2) as f32) * self.zoom;
            glow(p, glow_size, Color::new(1.0, 0.66, 0.18, 0.7));
            draw_circle(p.x, p.y, 2.0, WHITE);
        }
        for bolt in &self.bolts {
            let p = self.screen(bolt.pos);
            let from = self.screen(bolt.pos - bolt.vel * 0.035);
            draw_line(from.x, from.y, p.x, p.y,
                (3.0 * self.zoom).max(1.0), if bolt.friendly { SKYBLUE } else { RED });
        }
        for spark in &self.sparks {
            let p = self.screen(spark.pos);
            if p.x < -5.0 || p.x > screen_width() + 5.0 ||
                p.y < -5.0 || p.y > screen_height() + 5.0 { continue; }
            let alpha = spark.tint.a * (spark.life / spark.original_life);
            draw_circle(p.x, p.y, (spark.radius * self.zoom).max(0.5),
                Color::new(spark.tint.r, spark.tint.g, spark.tint.b, alpha));
        }
        if self.opponent.hp > 0.0 { self.draw_duelist(self.opponent, RED); }
        self.draw_duelist(self.player, SKYBLUE);

        draw_rectangle(16.0, 16.0, 570.0, 170.0, Color::new(0.015, 0.02, 0.042, 0.86));
        draw_rectangle_lines(16.0, 16.0, 570.0, 170.0, 1.0, DARKGRAY);
        draw_text("BLIGHTED GALAXY  /  RUST RECONSTRUCTION", 30.0, 45.0, 23.0, WHITE);
        draw_text("W/S thrust  A/D turn  SPACE plasma  E contact strike", 30.0, 77.0, 19.0, LIGHTGRAY);
        draw_text("B Hawking bomb  G swarm scatter  R restart  wheel zoom", 30.0, 104.0, 19.0, LIGHTGRAY);
        draw_text(&format!("HP {:.0} | rival {:.0} | swarm {} | particles {} | {:.0} fps",
            self.player.hp.max(0.0), self.opponent.hp.max(0.0), self.swarm.len(),
            self.sparks.len(), get_fps()), 30.0, 134.0, 19.0, SKYBLUE);
        draw_text(self.last_event, 30.0, 161.0, 18.0, GOLD);

        draw_text("TOY UNITS  /  NEW CLEAN-ROOM CODE  /  N-BODY + SWEPT COLLISION",
            22.0, screen_height() - 18.0, 18.0, GRAY);
    }
}

fn glow(p: Vec2, radius: f32, tint: Color) {
    for layer in (1..=5).rev() {
        let s = layer as f32;
        draw_circle(p.x, p.y, radius * (1.0 + s * 0.33),
            Color::new(tint.r, tint.g, tint.b, tint.a * 0.042));
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let mut game = Game::new();
    loop {
        if is_key_pressed(KeyCode::R) { game = Game::new(); }
        if is_key_pressed(KeyCode::Escape) { break; }
        let (_, wheel) = mouse_wheel();
        if wheel != 0.0 {
            game.zoom = (game.zoom * (1.0 + wheel * 0.09)).clamp(0.15, 2.2);
        }
        game.accumulator = (game.accumulator + f64::from(get_frame_time())).min(0.16);
        let mut steps = 0;
        while game.accumulator >= STEP && steps < 20 {
            game.step(STEP);
            game.accumulator -= STEP;
            steps += 1;
        }
        game.draw();
        next_frame().await;
    }
}
