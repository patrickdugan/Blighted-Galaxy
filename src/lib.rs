//! Deterministic, dependency-free simulation core for Blighted Galaxy.
//!
//! All positions, times, and masses are presently normalized gameplay units.
//! The limit LIGHT_SPEED is a gameplay constraint, rather than a unit conversion
//! for a physically calibrated general-relativistic simulation.

use std::f64::consts::TAU;
use serde::{Serialize,Deserialize};
use std::ops::{Add, AddAssign, Div, Mul, Sub, SubAssign};

pub const G: f64 = 1.0;
pub const SOFTENING: f64 = 12.0;
pub const LIGHT_SPEED: f64 = 1_800.0;

#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct V2 {
    pub x: f64,
    pub y: f64,
}

impl V2 {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }
    pub fn len_sq(self) -> f64 {
        self.dot(self)
    }
    pub fn len(self) -> f64 {
        self.len_sq().sqrt()
    }
    pub fn normalized(self) -> Self {
        let n = self.len();
        if n > 1e-12 { self / n } else { Self::default() }
    }
    pub fn clamp_len(self, max_length: f64) -> Self {
        let length = self.len();
        if length > max_length && length > 0.0 {
            self * (max_length / length)
        } else {
            self
        }
    }
    pub fn from_angle(angle: f64) -> Self {
        Self::new(angle.cos(), angle.sin())
    }
}
impl Add for V2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}
impl AddAssign for V2 {
    fn add_assign(&mut self, rhs: Self) { *self = *self + rhs; }
}
impl Sub for V2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y)
    }
}
impl SubAssign for V2 {
    fn sub_assign(&mut self, rhs: Self) { *self = *self - rhs; }
}
impl Mul<f64> for V2 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs)
    }
}
impl Div<f64> for V2 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BodyKind {
    Star,
    WhiteDwarf,
    NeutronStar,
    BlackHole,
}

#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub kind: BodyKind,
    pub pos: V2,
    pub vel: V2,
    pub mass: f64,
    pub radius: f64,
}
impl Body {
    pub fn new(kind: BodyKind, pos: V2, vel: V2, mass: f64, radius: f64) -> Self {
        assert!(mass > 0.0 && radius > 0.0);
        Self { kind, pos, vel, mass, radius }
    }
}

/// Plummer-softened inverse-square gravity. The softening radius keeps near
/// encounters finite while reserving full GR for a later module.
pub fn acceleration_at(point: V2, bodies: &[Body], excluded: Option<usize>) -> V2 {
    let mut acc = V2::default();
    for (index, body) in bodies.iter().enumerate() {
        if Some(index) == excluded { continue; }
        let delta = body.pos - point;
        let r2 = delta.len_sq() + SOFTENING * SOFTENING;
        acc += delta * (G * body.mass / (r2 * r2.sqrt()));
    }
    acc
}

/// Velocity Verlet is time-symmetric and has markedly less drift than Euler.
pub fn step_bodies(bodies: &mut [Body], dt: f64) {
    if dt <= 0.0 { return; }
    let initial: Vec<_> = bodies.iter().enumerate()
        .map(|(i, body)| acceleration_at(body.pos, bodies, Some(i)))
        .collect();
    for (body, a) in bodies.iter_mut().zip(initial) {
        body.pos += body.vel * dt + a * (0.5 * dt * dt);
        body.vel += a * (0.5 * dt);
    }
    let final_acc: Vec<_> = bodies.iter().enumerate()
        .map(|(i, body)| acceleration_at(body.pos, bodies, Some(i)))
        .collect();
    for (body, a) in bodies.iter_mut().zip(final_acc) {
        body.vel += a * (0.5 * dt);
    }
}

/// Dimensionless first-order tidal-acceleration / self-gravity ratio.
/// Values above one indicate a likely disruption in this simplified model.
pub fn tidal_ratio(disruptor: &Body, target: &Body) -> f64 {
    let distance = (target.pos - disruptor.pos).len()
        .max(disruptor.radius + target.radius);
    2.0 * disruptor.mass * target.radius.powi(3)
        / (target.mass * distance.powi(3))
}

/// Swept circle collision in the relative reference frame.
/// Returns a normalized time of impact in [0, 1], avoiding tunneling for
/// fast projectiles that never overlap at the discrete frame endpoints.
pub fn swept_circle_hit(
    a_start: V2, a_end: V2, a_radius: f64,
    b_start: V2, b_end: V2, b_radius: f64,
) -> Option<f64> {
    let start = a_start - b_start;
    let displacement = (a_end - a_start) - (b_end - b_start);
    let radius = (a_radius + b_radius).max(0.0);
    let c = start.len_sq() - radius * radius;
    if c <= 0.0 { return Some(0.0); }
    let a = displacement.len_sq();
    if a <= 1e-20 { return None; }
    let b = 2.0 * start.dot(displacement);
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 { return None; }
    let t = (-b - discriminant.sqrt()) / (2.0 * a);
    if (0.0..=1.0).contains(&t) { Some(t) } else { None }
}

#[derive(Clone, Copy, Debug)]
pub struct Aabb {
    pub min: V2,
    pub max: V2,
}
impl Aabb {
    pub fn centered(center: V2, half_size: V2) -> Self {
        Self { min: center - half_size, max: center + half_size }
    }
    pub fn overlaps(self, other: Self) -> bool {
        self.min.x <= other.max.x && self.max.x >= other.min.x
            && self.min.y <= other.max.y && self.max.y >= other.min.y
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Duelist {
    pub pos: V2,
    pub vel: V2,
    pub angle: f64,
    pub hp: f64,
    pub radius: f64,
}
impl Duelist {
    pub fn new(pos: V2, angle: f64) -> Self {
        Self { pos, vel: V2::default(), angle, hp: 100.0, radius: 16.0 }
    }
    pub fn heading(self) -> V2 { V2::from_angle(self.angle) }
    pub fn tick(&mut self, thrust: f64, gravity: V2, dt: f64) {
        self.vel = (self.vel + (self.heading() * thrust + gravity) * dt)
            .clamp_len(LIGHT_SPEED * 0.8);
        self.pos += self.vel * dt;
        self.hp = (self.hp + dt * 0.8).min(100.0);
    }
    pub fn box_hitbox(self) -> Aabb {
        Aabb::centered(self.pos, V2::new(self.radius, self.radius))
    }
    pub fn punch_hitbox(self) -> Aabb {
        Aabb::centered(self.pos + self.heading() * 34.0, V2::new(35.0, 30.0))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct HawkingBomb {
    pub pos: V2,
    pub vel: V2,
    pub mass: f64,
    pub evaporation_k: f64,
    pub detonated: bool,
}
impl HawkingBomb {
    pub fn new(pos: V2, vel: V2, mass: f64, evaporation_k: f64) -> Self {
        assert!(mass > 0.0 && evaporation_k > 0.0);
        Self { pos, vel, mass, evaporation_k, detonated: false }
    }
    pub fn lifetime(self) -> f64 {
        self.mass.powi(3) / (3.0 * self.evaporation_k)
    }
    pub fn luminosity(self) -> f64 {
        if self.detonated { 0.0 }
        else { self.evaporation_k / self.mass.powi(2) }
    }
    /// Exact solution of dm/dt = -k / m^2 for a normalized evaporator.
    /// Returns true only on the step that crosses total evaporation.
    pub fn step(&mut self, bodies: &[Body], dt: f64) -> bool {
        if self.detonated || dt <= 0.0 { return false; }
        let acceleration = acceleration_at(self.pos, bodies, None);
        self.vel = (self.vel + acceleration * dt).clamp_len(LIGHT_SPEED * 0.85);
        self.pos += self.vel * dt;
        let mass_cubed = self.mass.powi(3) - 3.0 * self.evaporation_k * dt;
        if mass_cubed <= 0.0 {
            self.mass = 0.0;
            self.detonated = true;
            true
        } else {
            self.mass = mass_cubed.cbrt();
            false
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SwarmNode {
    pub pos: V2,
    pub vel: V2,
}
fn rand01(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    (*seed >> 33) as f64 / (1_u64 << 31) as f64
}
/// Reproducible initial scatter, suitable for replayable benchmark scenarios.
pub fn scatter_swarm(center: V2, count: usize, seed: u64, radius: f64) -> Vec<SwarmNode> {
    let mut seed = seed;
    (0..count).map(|_| {
        let theta = rand01(&mut seed) * TAU;
        let radial = rand01(&mut seed).sqrt() * radius;
        let direction = V2::from_angle(theta);
        let velocity = direction * (10.0 + 55.0 * rand01(&mut seed));
        SwarmNode { pos: center + direction * radial, vel: velocity }
    }).collect()
}
pub fn step_swarm(nodes: &mut [SwarmNode], bodies: &[Body], dt: f64) {
    if dt <= 0.0 { return; }
    for node in nodes {
        node.vel = (node.vel + acceleration_at(node.pos, bodies, None) * dt)
            .clamp_len(LIGHT_SPEED * 0.8);
        node.pos += node.vel * dt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
    }

    #[test]
    fn zero_vector_normalization_is_finite() {
        assert_eq!(V2::default().normalized(), V2::default());
    }

    #[test]
    fn gravity_is_softened_at_the_source() {
        let body = Body::new(BodyKind::Star, V2::default(), V2::default(), 100.0, 12.0);
        let acc = acceleration_at(V2::default(), &[body], None);
        assert_eq!(acc, V2::default());
    }

    #[test]
    fn symmetric_bodies_accelerate_in_opposite_directions() {
        let bodies = vec![
            Body::new(BodyKind::Star, V2::new(-100.0, 0.0), V2::default(), 100.0, 12.0),
            Body::new(BodyKind::Star, V2::new(100.0, 0.0), V2::default(), 100.0, 12.0),
        ];
        let a = acceleration_at(bodies[0].pos, &bodies, Some(0));
        let b = acceleration_at(bodies[1].pos, &bodies, Some(1));
        near(a.x, -b.x, 1e-12);
        near(a.y, -b.y, 1e-12);
        assert!(a.x > 0.0);
    }

    #[test]
    fn zero_step_does_not_change_n_body_state() {
        let mut bodies = vec![Body::new(BodyKind::Star, V2::new(2.0, 4.0), V2::new(7.0, 11.0), 30.0, 4.0)];
        step_bodies(&mut bodies, 0.0);
        assert_eq!(bodies[0].pos, V2::new(2.0, 4.0));
        assert_eq!(bodies[0].vel, V2::new(7.0, 11.0));
    }

    #[test]
    fn projectile_sweep_catches_fast_crossing() {
        let hit = swept_circle_hit(V2::new(-100.0, 0.0), V2::new(100.0, 0.0), 2.0,
            V2::default(), V2::default(), 9.0);
        near(hit.expect("swept collision"), 0.445, 1e-10);
    }

    #[test]
    fn projectile_sweep_rejects_misses() {
        assert!(swept_circle_hit(V2::new(-100.0, 100.0), V2::new(100.0, 100.0), 2.0,
            V2::default(), V2::default(), 9.0).is_none());
    }

    #[test]
    fn projectile_sweep_reports_initial_overlap() {
        assert_eq!(swept_circle_hit(V2::default(), V2::new(100.0, 0.0), 2.0,
            V2::default(), V2::default(), 9.0), Some(0.0));
    }

    #[test]
    fn boxing_hitboxes_overlap() {
        let a = Aabb::centered(V2::default(), V2::new(10.0, 10.0));
        let b = Aabb::centered(V2::new(20.0, 0.0), V2::new(10.0, 10.0));
        assert!(a.overlaps(b));
    }

    #[test]
    fn tides_increase_rapidly_with_proximity() {
        let target = Body::new(BodyKind::WhiteDwarf, V2::default(), V2::default(), 50.0, 10.0);
        let far = Body::new(BodyKind::BlackHole, V2::new(100.0, 0.0), V2::default(), 100.0, 5.0);
        let near_body = Body { pos: V2::new(30.0, 0.0), ..far };
        assert!(tidal_ratio(&near_body, &target) > tidal_ratio(&far, &target));
    }

    #[test]
    fn hawking_evaporation_is_step_partition_independent() {
        let mut whole = HawkingBomb::new(V2::default(), V2::default(), 3.0, 0.5);
        let mut pieces = whole;
        assert!(!whole.step(&[], 5.0));
        for _ in 0..5 { assert!(!pieces.step(&[], 1.0)); }
        near(whole.mass, pieces.mass, 1e-12);
    }

    #[test]
    fn hawking_bomb_detonates_exactly_once() {
        let mut bomb = HawkingBomb::new(V2::default(), V2::default(), 1.0, 1.0);
        assert!(bomb.step(&[], 1.0));
        assert!(!bomb.step(&[], 1.0));
        assert!(bomb.detonated);
    }

    #[test]
    fn swarm_scatter_is_reproducible() {
        let a = scatter_swarm(V2::new(3.0, 4.0), 100, 42, 20.0);
        let b = scatter_swarm(V2::new(3.0, 4.0), 100, 42, 20.0);
        assert_eq!(a[12].pos, b[12].pos);
        assert_eq!(a[12].vel, b[12].vel);
    }

    #[test]
    fn duelist_speed_respects_gameplay_light_cone() {
        let mut player = Duelist::new(V2::default(), 0.0);
        player.tick(1e9, V2::default(), 1.0);
        assert!(player.vel.len() <= LIGHT_SPEED * 0.8 + 1e-9);
    }
}

// Calibrated SI reference physics is implemented independently of gameplay units.
pub mod physics;

// SI-to-gameplay adapter, relativistic combat and tidal-breakup demo.
pub mod si_bridge;

// Schwarzschild geodesics, bounded causal histories and symmetric 2D SPH debris.
pub mod schwarzschild;
pub mod causal;
pub mod debris;

// BG-PHY-05: Kerr equatorial light rays, 3D SPH fluid and GPU-backed particle meshes.
pub mod kerr;
pub mod sph3;
pub mod particle_batch;

// BG-PHY-06: cold electron-degenerate stellar structure reference.
pub mod white_dwarf;

// BG-PHY-06: non-equatorial Kerr null-ray camera.
pub mod kerr_camera;

// BG-PHY-07: curved Kerr optical acquisition of recorded orbital combatants.
pub mod kerr_optics;

// BG-PHY-07: damped isolated SPH white-dwarf relaxation diagnostics.
pub mod stellar_relax;

// BG-PHY-08: variable-h grad-h-corrected SPH force diagnostics.
pub mod sph_adaptive;

// BG-PHY-08: equatorial moving Kerr observers and aberrated causal optics.
pub mod kerr_moving_observer;

// BG-PHY-08: local emitter rest-frame Kerr target envelopes.
pub mod kerr_worldtube;

// BG-PHY-09: radial-shell-balanced degenerate stellar particle sampling.
pub mod sph_sampling;

// BG-PHY-10: Poisson-consistent M4 cubic spline softened gravitation.
pub mod sph_gravity_kernel;

// BG-PHY-11: deterministic compact-support SPH cell neighbors.
pub mod sph_neighbors;

// BG-PHY-11: mass-radius percentiles and hydrostatic profile diagnostics.
pub mod sph_massprofile;

// BG-PHY-12: hydrostatic radial force budgets and local virial diagnostics.
pub mod stellar_balance;

// BG-OVR-01: causal galactic campaign and self-rendered Macroquad navigation.
pub mod overworld;
pub mod overworld_nav;

// BG-STR-01: persistent industrial strategy and causal attack network.
pub mod strategy;

// BG-STR-03: reproducible encounter profiles from industrial star systems.
pub mod mission;

// BG-STR-05: versioned, validated, atomic campaign persistence.
pub mod campaign_save;
