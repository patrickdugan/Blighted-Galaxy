//! M3: Transhuman close-contact boxing and morphogenesis (GDD §4.2, §10).
//!
//! Models deformable bio-nanoborg duelists with distinct combat morphologies,
//! deterministic 2D uniform-grid broad-phase spatial partitioning, continuous
//! swept circle/AABB contact collision (zero tunneling at high closing speeds),
//! persistent viable mass / energy / wattage / matter / coherence / heat pools,
//! energy-limited regeneration, adaptive opponent AI, and memetic disruption
//! pulses that feed signed records back into the Patricia Trie on debrief.
use crate::physics::C_SI;
use crate::si_bridge::SiJoust;
use crate::strategy::TacticalOutcome;
use crate::{swept_circle_hit, Aabb, V2};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const CONTACT_ENTRY_MAX_RANGE_M: f64 = 80.0e6;
pub const CONTACT_ENTRY_MAX_REL_BETA: f64 = 0.015;
pub const ARENA_HALF_EXTENT_M: f64 = 120.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum MorphForm {
    /// High thrust and low profile; small collision cross-section.
    Streamlined,
    /// Wide damage-absorbing nanite barrier at reduced maneuver speed.
    AegisGuard,
    /// Extended kinetic striker limb with long continuous-swept reach.
    StrikerLimb,
    /// Short-range corrosive nanite cloud that erodes matter and shields.
    SwarmErode,
    /// High-wattage micro-singularity surge at sustained heat/coherence cost.
    SingularityCore,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MorphStats {
    pub body_radius_m: f64,
    pub thrust_accel_m_s2: f64,
    pub max_speed_m_s: f64,
    pub strike_reach_m: f64,
    pub strike_radius_m: f64,
    pub strike_damage_kg: f64,
    pub strike_energy_cost_j: f64,
    pub damage_mitigation: f64,
    pub wattage_w: f64,
    pub heat_per_s: f64,
    pub coherence_drain_per_s: f64,
    pub morph_matter_cost_kg: f64,
    pub erode_radius_m: f64,
    pub erode_dps_kg: f64,
}

impl MorphForm {
    pub const ALL: [MorphForm; 5] = [
        Self::Streamlined,
        Self::AegisGuard,
        Self::StrikerLimb,
        Self::SwarmErode,
        Self::SingularityCore,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::Streamlined => "1 STREAMLINED",
            Self::AegisGuard => "2 AEGIS GUARD",
            Self::StrikerLimb => "3 STRIKER LIMB",
            Self::SwarmErode => "4 SWARM ERODE",
            Self::SingularityCore => "5 SINGULARITY CORE",
        }
    }

    pub fn stats(self) -> MorphStats {
        match self {
            Self::Streamlined => MorphStats {
                body_radius_m: 2.8,
                thrust_accel_m_s2: 145.0,
                max_speed_m_s: 115.0,
                strike_reach_m: 8.5,
                strike_radius_m: 2.6,
                strike_damage_kg: 68.0,
                strike_energy_cost_j: 45.0e6,
                damage_mitigation: 0.12,
                wattage_w: 145.0e6,
                heat_per_s: 0.03,
                coherence_drain_per_s: 0.0,
                morph_matter_cost_kg: 8.0,
                erode_radius_m: 0.0,
                erode_dps_kg: 0.0,
            },
            Self::AegisGuard => MorphStats {
                body_radius_m: 5.4,
                thrust_accel_m_s2: 52.0,
                max_speed_m_s: 42.0,
                strike_reach_m: 7.0,
                strike_radius_m: 4.0,
                strike_damage_kg: 52.0,
                strike_energy_cost_j: 40.0e6,
                damage_mitigation: 0.65,
                wattage_w: 165.0e6,
                heat_per_s: 0.02,
                coherence_drain_per_s: 0.0,
                morph_matter_cost_kg: 14.0,
                erode_radius_m: 0.0,
                erode_dps_kg: 0.0,
            },
            Self::StrikerLimb => MorphStats {
                body_radius_m: 3.6,
                thrust_accel_m_s2: 95.0,
                max_speed_m_s: 78.0,
                strike_reach_m: 19.5,
                strike_radius_m: 3.8,
                strike_damage_kg: 142.0,
                strike_energy_cost_j: 85.0e6,
                damage_mitigation: 0.22,
                wattage_w: 175.0e6,
                heat_per_s: 0.05,
                coherence_drain_per_s: 0.0,
                morph_matter_cost_kg: 12.0,
                erode_radius_m: 0.0,
                erode_dps_kg: 0.0,
            },
            Self::SwarmErode => MorphStats {
                body_radius_m: 4.5,
                thrust_accel_m_s2: 80.0,
                max_speed_m_s: 64.0,
                strike_reach_m: 11.0,
                strike_radius_m: 4.5,
                strike_damage_kg: 76.0,
                strike_energy_cost_j: 55.0e6,
                damage_mitigation: 0.28,
                wattage_w: 160.0e6,
                heat_per_s: 0.06,
                coherence_drain_per_s: 0.008,
                morph_matter_cost_kg: 16.0,
                erode_radius_m: 22.0,
                erode_dps_kg: 78.0,
            },
            Self::SingularityCore => MorphStats {
                body_radius_m: 3.2,
                thrust_accel_m_s2: 125.0,
                max_speed_m_s: 98.0,
                strike_reach_m: 13.5,
                strike_radius_m: 3.4,
                strike_damage_kg: 178.0,
                strike_energy_cost_j: 110.0e6,
                damage_mitigation: 0.18,
                wattage_w: 440.0e6,
                heat_per_s: 0.22,
                coherence_drain_per_s: 0.045,
                morph_matter_cost_kg: 18.0,
                erode_radius_m: 0.0,
                erode_dps_kg: 0.0,
            },
        }
    }
}

/// Deterministic 2D uniform spatial hash grid for broad-phase contact filtering.
/// Uses `BTreeMap` so cell iteration and candidate pairs are strictly reproducible.
#[derive(Clone, Debug)]
pub struct UniformGrid2D {
    pub cell_size: f64,
    cells: BTreeMap<(i64, i64), Vec<usize>>,
}

impl UniformGrid2D {
    pub fn new(cell_size: f64) -> Self {
        assert!(cell_size.is_finite() && cell_size > 0.0);
        Self {
            cell_size,
            cells: BTreeMap::new(),
        }
    }

    pub fn clear(&mut self) {
        self.cells.clear();
    }

    fn cell_coord(&self, v: f64) -> i64 {
        (v / self.cell_size).floor() as i64
    }

    pub fn insert_aabb(&mut self, id: usize, aabb: Aabb) {
        let min_x = self.cell_coord(aabb.min.x);
        let max_x = self.cell_coord(aabb.max.x);
        let min_y = self.cell_coord(aabb.min.y);
        let max_y = self.cell_coord(aabb.max.y);
        for gx in min_x..=max_x {
            for gy in min_y..=max_y {
                self.cells.entry((gx, gy)).or_default().push(id);
            }
        }
    }

    pub fn query_aabb(&self, aabb: Aabb) -> Vec<usize> {
        let min_x = self.cell_coord(aabb.min.x);
        let max_x = self.cell_coord(aabb.max.x);
        let min_y = self.cell_coord(aabb.min.y);
        let max_y = self.cell_coord(aabb.max.y);
        let mut found = BTreeSet::new();
        for gx in min_x..=max_x {
            for gy in min_y..=max_y {
                if let Some(bucket) = self.cells.get(&(gx, gy)) {
                    for &id in bucket {
                        found.insert(id);
                    }
                }
            }
        }
        found.into_iter().collect()
    }

    pub fn candidate_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = BTreeSet::new();
        for bucket in self.cells.values() {
            for i in 0..bucket.len() {
                for j in (i + 1)..bucket.len() {
                    let a = bucket[i];
                    let b = bucket[j];
                    if a != b {
                        pairs.insert(if a < b { (a, b) } else { (b, a) });
                    }
                }
            }
        }
        pairs.into_iter().collect()
    }

    pub fn occupied_cells(&self) -> usize {
        self.cells.len()
    }
}

/// Swept bounding box for inserting moving circles/limbs into `UniformGrid2D`.
pub fn swept_aabb(start: V2, end: V2, radius: f64) -> Aabb {
    let r = radius.max(0.0);
    Aabb {
        min: V2::new(start.x.min(end.x) - r, start.y.min(end.y) - r),
        max: V2::new(start.x.max(end.x) + r, start.y.max(end.y) + r),
    }
}

/// Continuous swept AABB collision solver in the relative frame using Minkowski
/// sum slab intersection. Returns normalized time of impact `t in [0, 1]`.
pub fn swept_aabb_hit(
    a_start: V2,
    a_end: V2,
    a_half: V2,
    b_start: V2,
    b_end: V2,
    b_half: V2,
) -> Option<f64> {
    let half = V2::new(
        (a_half.x + b_half.x).max(0.0),
        (a_half.y + b_half.y).max(0.0),
    );
    let p0 = a_start - b_start;
    let rel_v = (a_end - a_start) - (b_end - b_start);

    if p0.x.abs() <= half.x && p0.y.abs() <= half.y {
        return Some(0.0);
    }

    let mut t_enter: f64 = 0.0;
    let mut t_exit: f64 = 1.0;

    for (pos, vel, ext) in [(p0.x, rel_v.x, half.x), (p0.y, rel_v.y, half.y)] {
        if vel.abs() <= 1e-12 {
            if pos.abs() > ext {
                return None;
            }
        } else {
            let t1 = (-ext - pos) / vel;
            let t2 = (ext - pos) / vel;
            let axis_enter = t1.min(t2);
            let axis_exit = t1.max(t2);
            t_enter = t_enter.max(axis_enter);
            t_exit = t_exit.min(axis_exit);
            if t_enter > t_exit {
                return None;
            }
        }
    }

    if (0.0..=1.0).contains(&t_enter) {
        Some(t_enter)
    } else {
        None
    }
}

/// Velocity-matching gate for transitioning from long-range relativistic joust
/// into close-contact transhuman boxing (GDD §2, §3.1).
pub fn can_enter_contact_boxing(joust: &SiJoust) -> bool {
    if !joust.player.alive() || !joust.opponent.alive() {
        return false;
    }
    let separation_m = (joust.opponent.pos_m - joust.player.pos_m).len();
    let rel_speed_m_s = (joust.opponent.vel_m_s - joust.player.vel_m_s).len();
    separation_m <= CONTACT_ENTRY_MAX_RANGE_M && rel_speed_m_s <= CONTACT_ENTRY_MAX_REL_BETA * C_SI
}

#[derive(Clone, Debug)]
pub struct MorphDuelist {
    pub pos: V2,
    pub vel: V2,
    pub facing_rad: f64,
    pub form: MorphForm,
    pub mass_kg: f64,
    pub max_mass_kg: f64,
    pub energy_j: f64,
    pub max_energy_j: f64,
    pub wattage_w: f64,
    pub matter_kg: f64,
    pub max_matter_kg: f64,
    pub coherence: f64,
    pub heat: f64,
    pub morph_cooldown_s: f64,
    pub strike_cooldown_s: f64,
    pub pulse_cooldown_s: f64,
    pub stagger_s: f64,
}

impl MorphDuelist {
    pub fn new(pos: V2, facing_rad: f64, form: MorphForm) -> Self {
        let stats = form.stats();
        Self {
            pos,
            vel: V2::default(),
            facing_rad,
            form,
            mass_kg: 1_000.0,
            max_mass_kg: 1_000.0,
            energy_j: 850.0e6,
            max_energy_j: 1_000.0e6,
            wattage_w: stats.wattage_w,
            matter_kg: 400.0,
            max_matter_kg: 400.0,
            coherence: 1.0,
            heat: 0.0,
            morph_cooldown_s: 0.0,
            strike_cooldown_s: 0.0,
            pulse_cooldown_s: 0.0,
            stagger_s: 0.0,
        }
    }

    pub fn stats(&self) -> MorphStats {
        self.form.stats()
    }

    pub fn alive(&self) -> bool {
        self.mass_kg > 0.0 && self.coherence > 0.05
    }

    pub fn heading(&self) -> V2 {
        V2::from_angle(self.facing_rad)
    }

    pub fn body_aabb(&self) -> Aabb {
        let r = self.stats().body_radius_m;
        Aabb::centered(self.pos, V2::new(r, r))
    }

    /// Refactor body morphology if cooldown, stagger, matter and energy permit.
    pub fn morph_to(&mut self, target: MorphForm) -> bool {
        if self.form == target || self.morph_cooldown_s > 0.0 || self.stagger_s > 0.0 {
            return false;
        }
        let next = target.stats();
        if self.matter_kg < next.morph_matter_cost_kg || self.energy_j < 25.0e6 {
            return false;
        }
        self.matter_kg -= next.morph_matter_cost_kg;
        self.energy_j -= 25.0e6;
        self.form = target;
        self.wattage_w = next.wattage_w;
        self.morph_cooldown_s = 0.55;
        true
    }

    /// Apply thermal, wattage, and energy-limited structural regeneration updates.
    pub fn tick_resources(&mut self, dt: f64) {
        if dt <= 0.0 || !self.alive() {
            return;
        }
        let s = self.stats();
        self.wattage_w = s.wattage_w;
        self.morph_cooldown_s = (self.morph_cooldown_s - dt).max(0.0);
        self.strike_cooldown_s = (self.strike_cooldown_s - dt).max(0.0);
        self.pulse_cooldown_s = (self.pulse_cooldown_s - dt).max(0.0);
        self.stagger_s = (self.stagger_s - dt).max(0.0);

        // Energy recharge bounded by wattage and max capacitor.
        self.energy_j = (self.energy_j + s.wattage_w * dt).min(self.max_energy_j);

        // Thermal balance: radiative rejection vs morph heat load.
        let cooling_per_s = 0.10;
        self.heat = (self.heat + (s.heat_per_s - cooling_per_s) * dt).clamp(0.0, 1.0);

        // SingularityCore and overheating erode cognitive coherence.
        let overheat_drain = if self.heat > 0.78 {
            (self.heat - 0.78) * 0.25
        } else {
            0.0
        };
        let total_drain = s.coherence_drain_per_s + overheat_drain;
        if total_drain > 0.0 {
            self.coherence = (self.coherence - total_drain * dt).clamp(0.0, 1.0);
        } else if self.heat < 0.45 && self.stagger_s <= 0.0 {
            self.coherence = (self.coherence + 0.015 * dt).min(1.0);
        }

        // Energy- and matter-limited structural regeneration.
        let energy_reserve = 40.0e6;
        let energy_per_kg = 6.0e6;
        if self.mass_kg < self.max_mass_kg
            && self.matter_kg > 0.0
            && self.energy_j > energy_reserve
            && self.heat < 0.85
        {
            let max_rate_kg = 24.0 * dt;
            let energy_limited_kg = (self.energy_j - energy_reserve) / energy_per_kg;
            let regen_kg = (self.max_mass_kg - self.mass_kg)
                .min(self.matter_kg)
                .min(max_rate_kg)
                .min(energy_limited_kg)
                .max(0.0);
            if regen_kg > 0.0 {
                self.mass_kg += regen_kg;
                self.matter_kg -= regen_kg;
                self.energy_j -= regen_kg * energy_per_kg;
                self.heat = (self.heat + regen_kg * 0.002).min(1.0);
            }
        }
    }

    pub fn apply_damage(&mut self, raw_kg: f64, bypass_mitigation: f64) -> f64 {
        let s = self.stats();
        let eff_mitigation = (s.damage_mitigation * (1.0 - bypass_mitigation.clamp(0.0, 1.0)))
            .clamp(0.0, 0.85);
        // Low coherence increases structural vulnerability.
        let coherence_penalty = 1.0 + (1.0 - self.coherence) * 0.45;
        let actual = (raw_kg * (1.0 - eff_mitigation) * coherence_penalty).max(0.0);
        self.mass_kg = (self.mass_kg - actual).max(0.0);
        actual
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct BoxingInput {
    pub move_dir: V2,
    pub facing_rad: f64,
    pub requested_form: Option<MorphForm>,
    pub strike: bool,
    pub memetic_pulse: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct StrikeSweep {
    pub from_player: bool,
    pub start: V2,
    pub end: V2,
    pub radius_m: f64,
    pub age_s: f64,
    pub hit: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct MemeticPulse {
    pub from_player: bool,
    pub origin: V2,
    pub radius_m: f64,
    pub max_radius_m: f64,
    pub age_s: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ContactSpark {
    pub pos: V2,
    pub age_s: f64,
    pub intensity: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct BoxingDebrief {
    pub site: Option<usize>,
    pub outcome: TacticalOutcome,
    pub memetic_pulses_landed: u32,
    pub player_coherence: f64,
    pub opponent_coherence: f64,
    pub contacts_resolved: u64,
}

pub struct MorphBoxingArena {
    pub time_s: f64,
    pub player: MorphDuelist,
    pub opponent: MorphDuelist,
    pub grid: UniformGrid2D,
    pub strikes: Vec<StrikeSweep>,
    pub pulses: Vec<MemeticPulse>,
    pub sparks: Vec<ContactSpark>,
    pub contacts_resolved: u64,
    pub memetic_pulses_landed: u32,
    pub site: Option<usize>,
    pub site_name: Option<&'static str>,
}

impl Default for MorphBoxingArena {
    fn default() -> Self {
        Self::new(None, None)
    }
}

impl MorphBoxingArena {
    pub fn new(site: Option<usize>, site_name: Option<&'static str>) -> Self {
        Self {
            time_s: 0.0,
            player: MorphDuelist::new(V2::new(-38.0, 0.0), 0.0, MorphForm::Streamlined),
            opponent: MorphDuelist::new(
                V2::new(38.0, 0.0),
                std::f64::consts::PI,
                MorphForm::StrikerLimb,
            ),
            grid: UniformGrid2D::new(16.0),
            strikes: Vec::new(),
            pulses: Vec::new(),
            sparks: Vec::new(),
            contacts_resolved: 0,
            memetic_pulses_landed: 0,
            site,
            site_name,
        }
    }

    /// Deterministic mid-exchange showcase fixture for `--demo-boxing` QA frames.
    pub fn demo_showcase() -> Self {
        let mut arena = Self::new(Some(5), Some("BEL-AIR"));
        arena.player.form = MorphForm::StrikerLimb;
        arena.player.pos = V2::new(-16.0, -4.0);
        arena.player.facing_rad = 0.14;
        arena.player.mass_kg = 860.0;
        arena.player.energy_j = 710.0e6;
        arena.player.matter_kg = 325.0;
        arena.player.heat = 0.34;

        arena.opponent.form = MorphForm::AegisGuard;
        arena.opponent.pos = V2::new(14.0, 2.0);
        arena.opponent.facing_rad = std::f64::consts::PI + 0.12;
        arena.opponent.mass_kg = 640.0;
        arena.opponent.energy_j = 580.0e6;
        arena.opponent.matter_kg = 270.0;
        arena.opponent.coherence = 0.68;
        arena.opponent.heat = 0.42;

        arena.strikes.push(StrikeSweep {
            from_player: true,
            start: arena.player.pos,
            end: arena.player.pos + arena.player.heading() * 19.5,
            radius_m: 3.8,
            age_s: 0.06,
            hit: true,
        });
        arena.pulses.push(MemeticPulse {
            from_player: true,
            origin: arena.player.pos,
            radius_m: 28.0,
            max_radius_m: 62.0,
            age_s: 0.18,
        });
        arena.sparks.push(ContactSpark {
            pos: V2::new(2.0, -1.0),
            age_s: 0.05,
            intensity: 1.0,
        });
        arena.memetic_pulses_landed = 1;
        arena.contacts_resolved = 3;
        arena.time_s = 12.4;
        arena
    }

    pub fn outcome(&self) -> TacticalOutcome {
        if !self.player.alive() {
            TacticalOutcome::Defeat
        } else if !self.opponent.alive() {
            TacticalOutcome::Victory
        } else {
            TacticalOutcome::Withdrawn
        }
    }

    pub fn debrief(&self) -> BoxingDebrief {
        BoxingDebrief {
            site: self.site,
            outcome: self.outcome(),
            memetic_pulses_landed: self.memetic_pulses_landed,
            player_coherence: self.player.coherence,
            opponent_coherence: self.opponent.coherence,
            contacts_resolved: self.contacts_resolved,
        }
    }

    fn choose_opponent_input(&self) -> BoxingInput {
        let to_player = self.player.pos - self.opponent.pos;
        let dist = to_player.len();
        let dir = to_player.normalized();
        let facing_rad = dir.y.atan2(dir.x);

        let desired_form = if dist > 48.0 {
            MorphForm::Streamlined
        } else if self.player.form == MorphForm::AegisGuard && dist < 24.0 {
            MorphForm::SwarmErode
        } else if self.player.form == MorphForm::StrikerLimb
            && self.opponent.mass_kg < self.opponent.max_mass_kg * 0.65
            && dist < 25.0
        {
            MorphForm::AegisGuard
        } else if self.opponent.mass_kg < 420.0
            && self.opponent.coherence > 0.45
            && self.opponent.heat < 0.65
        {
            MorphForm::SingularityCore
        } else {
            MorphForm::StrikerLimb
        };

        let move_dir = if dist > self.opponent.stats().strike_reach_m * 0.85 {
            dir
        } else if dist < 9.0 {
            V2::new(-dir.y, dir.x) * 0.6 - dir * 0.4
        } else {
            V2::new(-dir.y, dir.x) * 0.8 + dir * 0.2
        };

        let strike = dist <= self.opponent.stats().strike_reach_m + self.player.stats().body_radius_m + 2.0
            && self.opponent.strike_cooldown_s <= 0.0;
        let memetic_pulse = dist <= 45.0
            && self.opponent.pulse_cooldown_s <= 0.0
            && self.opponent.energy_j >= 140.0e6
            && self.player.coherence > 0.55;

        BoxingInput {
            move_dir,
            facing_rad,
            requested_form: (desired_form != self.opponent.form).then_some(desired_form),
            strike,
            memetic_pulse,
        }
    }

    fn execute_strike(
        attacker: &mut MorphDuelist,
        defender: &mut MorphDuelist,
        attacker_prev: V2,
        defender_prev: V2,
        from_player: bool,
        grid: &UniformGrid2D,
        strikes: &mut Vec<StrikeSweep>,
        sparks: &mut Vec<ContactSpark>,
        contacts: &mut u64,
    ) {
        let s = attacker.stats();
        if attacker.strike_cooldown_s > 0.0
            || attacker.stagger_s > 0.0
            || attacker.energy_j < s.strike_energy_cost_j
        {
            return;
        }
        attacker.energy_j -= s.strike_energy_cost_j;
        attacker.strike_cooldown_s = 0.42;
        attacker.heat = (attacker.heat + 0.06).min(1.0);

        let limb_start = attacker_prev;
        let limb_end = attacker.pos + attacker.heading() * s.strike_reach_m;
        let sweep_box = swept_aabb(limb_start, limb_end, s.strike_radius_m);
        let defender_id = if from_player { 1usize } else { 0usize };
        let broad_candidates = grid.query_aabb(sweep_box);

        let mut landed = false;
        if broad_candidates.contains(&defender_id) {
            let def_r = defender.stats().body_radius_m;
            let circle_toi = swept_circle_hit(
                limb_start,
                limb_end,
                s.strike_radius_m,
                defender_prev,
                defender.pos,
                def_r,
            );
            let aabb_toi = swept_aabb_hit(
                limb_start,
                limb_end,
                V2::new(s.strike_radius_m, s.strike_radius_m),
                defender_prev,
                defender.pos,
                V2::new(def_r, def_r),
            );
            if let Some(toi) = circle_toi.or(aabb_toi) {
                landed = true;
                *contacts += 1;
                let impact_pos = limb_start + (limb_end - limb_start) * toi;
                let rel_speed = (attacker.vel - defender.vel).len();
                let kinetic_bonus = 1.0 + (rel_speed / 120.0).min(0.8);
                let bypass = if attacker.form == MorphForm::SingularityCore {
                    0.35
                } else {
                    0.0
                };
                defender.apply_damage(s.strike_damage_kg * kinetic_bonus, bypass);
                let knock_dir = (defender.pos - attacker.pos).normalized();
                defender.vel += knock_dir * 26.0;
                sparks.push(ContactSpark {
                    pos: impact_pos,
                    age_s: 0.0,
                    intensity: 1.0,
                });
            }
        }
        strikes.push(StrikeSweep {
            from_player,
            start: limb_start,
            end: limb_end,
            radius_m: s.strike_radius_m,
            age_s: 0.0,
            hit: landed,
        });
    }

    fn execute_memetic_pulse(
        caster: &mut MorphDuelist,
        target: &mut MorphDuelist,
        from_player: bool,
        pulses: &mut Vec<MemeticPulse>,
        sparks: &mut Vec<ContactSpark>,
        landed_counter: &mut u32,
    ) {
        let energy_cost = 120.0e6;
        if caster.pulse_cooldown_s > 0.0 || caster.stagger_s > 0.0 || caster.energy_j < energy_cost
        {
            return;
        }
        caster.energy_j -= energy_cost;
        caster.pulse_cooldown_s = 2.6;
        let max_range = 62.0;
        pulses.push(MemeticPulse {
            from_player,
            origin: caster.pos,
            radius_m: 4.0,
            max_radius_m: max_range,
            age_s: 0.0,
        });
        let dist = (target.pos - caster.pos).len();
        if dist <= max_range {
            let falloff = (1.0 - dist / max_range).clamp(0.25, 1.0);
            let impact = 0.24 * falloff * (0.65 + 0.35 * caster.coherence);
            target.coherence = (target.coherence - impact).max(0.0);
            target.stagger_s = target.stagger_s.max(0.55 * falloff);
            if from_player {
                *landed_counter += 1;
            }
            sparks.push(ContactSpark {
                pos: target.pos,
                age_s: 0.0,
                intensity: 0.85,
            });
        }
    }

    fn integrate_duelist(duelist: &mut MorphDuelist, input: BoxingInput, dt: f64) {
        if let Some(req) = input.requested_form {
            duelist.morph_to(req);
        }
        duelist.facing_rad = input.facing_rad;
        let s = duelist.stats();
        let thrust_dir = input.move_dir.normalized();
        if thrust_dir.len_sq() > 1e-9 && duelist.energy_j > 5.0e6 {
            duelist.vel += thrust_dir * (s.thrust_accel_m_s2 * dt);
            duelist.energy_j = (duelist.energy_j - 18.0e6 * dt).max(0.0);
        } else {
            // Mild inertial drag / stabilization jets.
            duelist.vel = duelist.vel * (1.0 - (1.8 * dt).min(0.5));
        }
        duelist.vel = duelist.vel.clamp_len(s.max_speed_m_s);
        duelist.pos += duelist.vel * dt;

        // Soft containment within close-contact bounding arena.
        for (coord, vel) in [
            (&mut duelist.pos.x, &mut duelist.vel.x),
            (&mut duelist.pos.y, &mut duelist.vel.y),
        ] {
            if *coord < -ARENA_HALF_EXTENT_M {
                *coord = -ARENA_HALF_EXTENT_M;
                *vel = vel.abs() * 0.5;
            } else if *coord > ARENA_HALF_EXTENT_M {
                *coord = ARENA_HALF_EXTENT_M;
                *vel = -vel.abs() * 0.5;
            }
        }
    }

    pub fn step(&mut self, dt: f64, player_input: BoxingInput) {
        if dt <= 0.0 || !dt.is_finite() {
            return;
        }
        let dt = dt.min(0.05);
        if !self.player.alive() || !self.opponent.alive() {
            self.time_s += dt;
            return;
        }

        let ai_input = self.choose_opponent_input();
        let p_prev = self.player.pos;
        let o_prev = self.opponent.pos;

        self.player.tick_resources(dt);
        self.opponent.tick_resources(dt);
        Self::integrate_duelist(&mut self.player, player_input, dt);
        Self::integrate_duelist(&mut self.opponent, ai_input, dt);

        // Rebuild UniformGrid2D broad-phase with swept body volumes.
        self.grid.clear();
        self.grid.insert_aabb(
            0,
            swept_aabb(p_prev, self.player.pos, self.player.stats().body_radius_m),
        );
        self.grid.insert_aabb(
            1,
            swept_aabb(
                o_prev,
                self.opponent.pos,
                self.opponent.stats().body_radius_m,
            ),
        );

        // Continuous body-to-body contact resolution via broad-phase + swept circle/AABB.
        if self.grid.candidate_pairs().contains(&(0, 1)) {
            let pr = self.player.stats().body_radius_m;
            let or = self.opponent.stats().body_radius_m;
            if let Some(toi) = swept_circle_hit(
                p_prev,
                self.player.pos,
                pr,
                o_prev,
                self.opponent.pos,
                or,
            ) {
                self.contacts_resolved += 1;
                let pa = p_prev + (self.player.pos - p_prev) * toi;
                let pb = o_prev + (self.opponent.pos - o_prev) * toi;
                let mut normal = (pb - pa).normalized();
                if normal.len_sq() < 1e-9 {
                    normal = V2::new(1.0, 0.0);
                }
                let rel_v = self.player.vel - self.opponent.vel;
                let closing = rel_v.dot(normal);
                if closing > 0.0 {
                    let impact_dmg = (closing * 0.45).min(60.0);
                    self.player.apply_damage(impact_dmg, 0.0);
                    self.opponent.apply_damage(impact_dmg, 0.0);
                    let impulse = normal * (closing * 0.75);
                    self.player.vel -= impulse;
                    self.opponent.vel += impulse;
                    self.sparks.push(ContactSpark {
                        pos: (pa + pb) * 0.5,
                        age_s: 0.0,
                        intensity: (closing / 80.0).clamp(0.3, 1.0),
                    });
                }
                // Separate overlapping cores cleanly.
                let sep = (self.opponent.pos - self.player.pos).len();
                let min_sep = pr + or;
                if sep < min_sep {
                    let push = normal * ((min_sep - sep) * 0.55);
                    self.player.pos -= push;
                    self.opponent.pos += push;
                }
            }
        }

        // SwarmErode nanite cloud continuous proximity erosion.
        let dist_now = (self.opponent.pos - self.player.pos).len();
        let ps = self.player.stats();
        if ps.erode_radius_m > 0.0 && dist_now <= ps.erode_radius_m && self.player.matter_kg > 0.0 {
            let eroded = ps.erode_dps_kg * dt;
            self.opponent.apply_damage(eroded, 0.55);
            self.opponent.matter_kg = (self.opponent.matter_kg - 12.0 * dt).max(0.0);
            self.player.matter_kg = (self.player.matter_kg - 5.0 * dt).max(0.0);
        }
        let os = self.opponent.stats();
        if os.erode_radius_m > 0.0 && dist_now <= os.erode_radius_m && self.opponent.matter_kg > 0.0
        {
            let eroded = os.erode_dps_kg * dt;
            self.player.apply_damage(eroded, 0.55);
            self.player.matter_kg = (self.player.matter_kg - 12.0 * dt).max(0.0);
            self.opponent.matter_kg = (self.opponent.matter_kg - 5.0 * dt).max(0.0);
        }

        // Continuous-swept limb strikes.
        if player_input.strike {
            Self::execute_strike(
                &mut self.player,
                &mut self.opponent,
                p_prev,
                o_prev,
                true,
                &self.grid,
                &mut self.strikes,
                &mut self.sparks,
                &mut self.contacts_resolved,
            );
        }
        if ai_input.strike && self.opponent.alive() {
            Self::execute_strike(
                &mut self.opponent,
                &mut self.player,
                o_prev,
                p_prev,
                false,
                &self.grid,
                &mut self.strikes,
                &mut self.sparks,
                &mut self.contacts_resolved,
            );
        }

        // Memetic disruption pulses.
        if player_input.memetic_pulse {
            Self::execute_memetic_pulse(
                &mut self.player,
                &mut self.opponent,
                true,
                &mut self.pulses,
                &mut self.sparks,
                &mut self.memetic_pulses_landed,
            );
        }
        if ai_input.memetic_pulse && self.opponent.alive() {
            let mut ai_counter = 0u32;
            Self::execute_memetic_pulse(
                &mut self.opponent,
                &mut self.player,
                false,
                &mut self.pulses,
                &mut self.sparks,
                &mut ai_counter,
            );
        }

        // Age transient visual/telemetry effects.
        for s in &mut self.strikes {
            s.age_s += dt;
        }
        self.strikes.retain(|s| s.age_s < 0.28);

        for p in &mut self.pulses {
            p.age_s += dt;
            p.radius_m = (p.radius_m + 145.0 * dt).min(p.max_radius_m);
        }
        self.pulses.retain(|p| p.age_s < 0.45);

        for sp in &mut self.sparks {
            sp.age_s += dt;
        }
        self.sparks.retain(|sp| sp.age_s < 0.35);

        self.time_s += dt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_grid_2d_broad_phase_is_deterministic_and_finds_overlaps() {
        let mut g = UniformGrid2D::new(10.0);
        g.insert_aabb(0, Aabb::centered(V2::new(0.0, 0.0), V2::new(6.0, 6.0)));
        g.insert_aabb(1, Aabb::centered(V2::new(8.0, 2.0), V2::new(4.0, 4.0)));
        g.insert_aabb(2, Aabb::centered(V2::new(100.0, 100.0), V2::new(3.0, 3.0)));
        assert_eq!(g.candidate_pairs(), vec![(0, 1)]);
        assert_eq!(
            g.query_aabb(Aabb::centered(V2::new(5.0, 0.0), V2::new(2.0, 2.0))),
            vec![0, 1]
        );
    }

    #[test]
    fn swept_aabb_catches_high_speed_crossing_with_zero_tunneling() {
        // Object moves 2,000 m in one step past a 4 m stationary target.
        let toi = swept_aabb_hit(
            V2::new(-1_000.0, 0.0),
            V2::new(1_000.0, 0.0),
            V2::new(2.0, 2.0),
            V2::default(),
            V2::default(),
            V2::new(4.0, 4.0),
        );
        let t = toi.expect("must detect swept AABB impact without tunneling");
        assert!((t - 0.497).abs() < 1e-6);

        // Parallel miss must return None.
        assert!(swept_aabb_hit(
            V2::new(-1_000.0, 20.0),
            V2::new(1_000.0, 20.0),
            V2::new(2.0, 2.0),
            V2::default(),
            V2::default(),
            V2::new(4.0, 4.0),
        )
        .is_none());
    }

    #[test]
    fn morph_forms_enforce_matter_cost_and_distinct_mitigation() {
        let mut a = MorphDuelist::new(V2::default(), 0.0, MorphForm::Streamlined);
        let mut b = MorphDuelist::new(V2::default(), 0.0, MorphForm::Streamlined);
        let start_matter = b.matter_kg;
        assert!(b.morph_to(MorphForm::AegisGuard));
        assert!(b.matter_kg < start_matter);

        let dmg_streamlined = a.apply_damage(100.0, 0.0);
        let dmg_aegis = b.apply_damage(100.0, 0.0);
        assert!(dmg_aegis < dmg_streamlined * 0.5);
    }

    #[test]
    fn energy_and_matter_limit_structural_regeneration() {
        let mut d = MorphDuelist::new(V2::default(), 0.0, MorphForm::Streamlined);
        d.mass_kg = 800.0;
        d.energy_j = 10.0e6; // below reserve -> cannot regenerate until wattage charges
        d.matter_kg = 100.0;
        d.wattage_w = 0.0;
        let m0 = d.mass_kg;
        // Manually test with depleted matter too:
        d.matter_kg = 0.0;
        d.energy_j = 900.0e6;
        d.tick_resources(0.1);
        assert_eq!(d.mass_kg, m0);

        // With both energy and matter, regeneration consumes matter & energy.
        d.matter_kg = 50.0;
        let e_before = d.energy_j;
        d.tick_resources(0.2);
        assert!(d.mass_kg > m0);
        assert!(d.matter_kg < 50.0);
        assert!(d.energy_j <= d.max_energy_j);
        let _ = e_before;
    }

    #[test]
    fn velocity_matching_gate_requires_proximity_and_matched_speed() {
        let mut joust = SiJoust::new();
        // Default joust is 500 Mm apart -> too far for contact boxing.
        assert!(!can_enter_contact_boxing(&joust));

        // Close range but huge relative speed -> still blocked.
        joust.opponent.pos_m = joust.player.pos_m + V2::new(20.0e6, 0.0);
        joust.opponent.vel_m_s = V2::new(-0.35 * C_SI, 0.0);
        assert!(!can_enter_contact_boxing(&joust));

        // Matched velocity within 0.015c and <= 80 Mm -> allowed!
        joust.opponent.vel_m_s = joust.player.vel_m_s + V2::new(0.005 * C_SI, 0.0);
        assert!(can_enter_contact_boxing(&joust));
    }

    #[test]
    fn memetic_pulse_and_strikes_resolve_and_produce_debrief() {
        let mut arena = MorphBoxingArena::new(Some(5), Some("BEL-AIR"));
        arena.player.pos = V2::new(-10.0, 0.0);
        arena.opponent.pos = V2::new(6.0, 0.0);
        arena.player.form = MorphForm::StrikerLimb;
        let opp_mass_0 = arena.opponent.mass_kg;
        let opp_coh_0 = arena.opponent.coherence;

        arena.step(
            1.0 / 60.0,
            BoxingInput {
                move_dir: V2::new(1.0, 0.0),
                facing_rad: 0.0,
                requested_form: None,
                strike: true,
                memetic_pulse: true,
            },
        );
        assert!(arena.opponent.mass_kg < opp_mass_0);
        assert!(arena.opponent.coherence < opp_coh_0);
        assert_eq!(arena.memetic_pulses_landed, 1);
        assert!(arena.contacts_resolved >= 1);
        let debrief = arena.debrief();
        assert_eq!(debrief.site, Some(5));
        assert_eq!(debrief.memetic_pulses_landed, 1);
    }
}
