//! BG-STR-06: Episodes XVII-XXII campaign spine and Great Aunt Weena's
//! autonomous industrial / Chamber of Guf planner (GDD §8).
//!
//! Narrative provenance: episode titles and `[CANON]` summaries derive from
//! the recovered `macro-history.txt` and surviving manuscripts. Objective
//! thresholds and star-system triggers are explicitly marked `[STAGING]` as
//! playable campaign mechanics rather than recovered canon dates.
use crate::patricia_trie::PatriciaTrieWorld;
use crate::strategy::{Faction, Project, SignalKind, StrategicState, System};
use serde::{Deserialize, Serialize};

pub const EPISODE_COUNT: usize = 6;
pub const MAX_DISPATCHES: usize = 64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EpisodeSpec {
    pub index: usize,
    pub numeral: &'static str,
    pub title: &'static str,
    pub canon_basis: &'static str,
    pub staging_objective: &'static str,
}

pub const EPISODES: [EpisodeSpec; EPISODE_COUNT] = [
    EpisodeSpec {
        index: 0,
        numeral: "XVII",
        title: "You Never Go Full Blight",
        canon_basis: "[CANON] Remnant hides beyond the Outer Arm, fabricating kugelblitz black holes to reclaim Halo Drive capacity.",
        staging_objective: "[STAGING] Expand any Drummer system to IND >= 2.5 and contain >= 1 kugelblitz (or launch a Halo strike).",
    },
    EpisodeSpec {
        index: 1,
        numeral: "XVIII",
        title: "Patricia Tries",
        canon_basis: "[CANON] Alfonso and Tatyana establish first contact with local Patricia Tries nodes and intercept a Cluster Bomber.",
        staging_objective: "[STAGING] Win a tactical/boxing engagement and log >= 1 contact or memetic record in the Patricia Trie.",
    },
    EpisodeSpec {
        index: 2,
        numeral: "XIX",
        title: "Blowing Up The Galaxy",
        canon_basis: "[CANON] Great Aunt Weena constructs the Chamber of Guf in the Nebular Wastes from a 12 M_sun remnant.",
        staging_objective: "[STAGING] Defend a Drummer system with a Guf screen after an array fire, or capture Nebular Wastes.",
    },
    EpisodeSpec {
        index: 3,
        numeral: "XX",
        title: "Death Star Hop-Scotch",
        canon_basis: "[CANON] Cascading light-speed counter-attacks force coordinated anticipation and logistics across the Outer Arm.",
        staging_objective: "[STAGING] Hold Guf defense/Nebular Wastes and deliver or escort a Drummer convoy (or control >= 6 systems).",
    },
    EpisodeSpec {
        index: 4,
        numeral: "XXI",
        title: "Double-barreled White Dwarf Shotguns",
        canon_basis: "[CANON] Blight avatars corral white dwarfs for stellar-remnant sieges across the protected Outer Arm sectors.",
        staging_objective: "[STAGING] Win a White-Dwarf Tidal Interdiction or a tactical Halo Strike.",
    },
    EpisodeSpec {
        index: 5,
        numeral: "XXII",
        title: "Amicus Brief",
        canon_basis: "[CANON] Two Prime Branches of the Patricia Tries wage light-speed consensus war: coexistence vs assimilation.",
        staging_objective: "[STAGING] Reconcile the Patricia Trie to 1 unified root (Shift+M) with avg coexist >= 45%.",
    },
];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EpisodeTracker {
    pub active_index: usize,
    pub completed: [bool; EPISODE_COUNT],
    pub completion_years: [Option<f64>; EPISODE_COUNT],
    pub tactical_victories: u32,
    pub tidal_victories: u32,
    pub guf_defenses_won: u32,
    pub convoy_escorts_won: u32,
    pub halo_strikes_won: u32,
    pub boxing_duels_won: u32,
    pub convoys_delivered: u32,
    pub halo_strikes_launched: u32,
    pub dispatches: Vec<String>,
}

impl Default for EpisodeTracker {
    fn default() -> Self {
        Self {
            active_index: 0,
            completed: [false; EPISODE_COUNT],
            completion_years: [None; EPISODE_COUNT],
            tactical_victories: 0,
            tidal_victories: 0,
            guf_defenses_won: 0,
            convoy_escorts_won: 0,
            halo_strikes_won: 0,
            boxing_duels_won: 0,
            convoys_delivered: 0,
            halo_strikes_launched: 0,
            dispatches: vec![format!(
                "EP {}: {} — {}",
                EPISODES[0].numeral, EPISODES[0].title, EPISODES[0].staging_objective
            )],
        }
    }
}

impl EpisodeTracker {
    pub fn current_spec(&self) -> Option<&'static EpisodeSpec> {
        EPISODES.get(self.active_index)
    }

    pub fn completed_count(&self) -> usize {
        self.completed.iter().filter(|&&c| c).count()
    }

    pub fn all_complete(&self) -> bool {
        self.active_index >= EPISODE_COUNT
    }

    pub fn predicate_met(
        &self,
        index: usize,
        systems: &[System],
        pulses: u32,
        trie: &PatriciaTrieWorld,
    ) -> bool {
        match index {
            0 => {
                let has_kugelblitz = self.halo_strikes_launched >= 1
                    || systems
                        .iter()
                        .any(|s| s.owner == Faction::Drummers && s.kugelblitz >= 1);
                let has_industry = systems
                    .iter()
                    .any(|s| s.owner == Faction::Drummers && s.industry >= 2.5);
                has_kugelblitz && has_industry
            }
            1 => {
                (self.tactical_victories >= 1 || self.boxing_duels_won >= 1)
                    && (trie.contact_records >= 1 || trie.memetic_broadcasts >= 1)
            }
            2 => {
                let warded_after_fire = pulses >= 1
                    && (self.guf_defenses_won >= 1
                        || systems
                            .iter()
                            .any(|s| s.owner == Faction::Drummers && s.guf_ward > 0.0));
                let wastes_taken = systems
                    .get(6)
                    .is_some_and(|s| s.owner == Faction::Drummers);
                warded_after_fire || wastes_taken
            }
            3 => {
                let guf_countered = self.guf_defenses_won >= 1
                    || systems
                        .get(6)
                        .is_some_and(|s| s.owner == Faction::Drummers);
                let logistics_or_reach = self.convoys_delivered >= 1
                    || self.convoy_escorts_won >= 1
                    || systems
                        .iter()
                        .filter(|s| s.owner == Faction::Drummers)
                        .count()
                        >= 6;
                guf_countered && logistics_or_reach
            }
            4 => self.tidal_victories >= 1 || self.halo_strikes_won >= 1,
            5 => trie.branch_divergence_count() <= 1 && trie.average_coexist() >= 0.45,
            _ => false,
        }
    }

    /// Evaluate sequential episode progression against live strategic + trie state.
    /// Returns any newly issued narrative dispatches.
    pub fn evaluate(
        &mut self,
        systems: &[System],
        pulses: u32,
        year: f64,
        trie: &PatriciaTrieWorld,
    ) -> Vec<String> {
        let mut issued = Vec::new();
        while self.active_index < EPISODE_COUNT {
            let idx = self.active_index;
            if !self.predicate_met(idx, systems, pulses, trie) {
                break;
            }
            self.completed[idx] = true;
            self.completion_years[idx] = Some(year);
            let spec = &EPISODES[idx];
            let done_msg = format!(
                "EP {} COMPLETE (T+{year:.0}): {} — {}",
                spec.numeral, spec.title, spec.canon_basis
            );
            self.push_dispatch(done_msg.clone());
            issued.push(done_msg);
            self.active_index += 1;
            if let Some(next) = EPISODES.get(self.active_index) {
                let next_msg = format!(
                    "EP {} ACTIVE: {} — {}",
                    next.numeral, next.title, next.staging_objective
                );
                self.push_dispatch(next_msg.clone());
                issued.push(next_msg);
            }
        }
        issued
    }

    fn push_dispatch(&mut self, msg: String) {
        self.dispatches.push(msg);
        if self.dispatches.len() > MAX_DISPATCHES {
            let extra = self.dispatches.len() - MAX_DISPATCHES;
            self.dispatches.drain(..extra);
        }
    }

    pub fn valid(&self, year: f64) -> bool {
        if self.active_index > EPISODE_COUNT
            || self.dispatches.len() > MAX_DISPATCHES
            || self.dispatches.iter().any(|d| d.len() > 384)
        {
            return false;
        }
        for i in 0..EPISODE_COUNT {
            if i < self.active_index && !self.completed[i] {
                return false;
            }
            if let Some(t) = self.completion_years[i] {
                if !self.completed[i] || !t.is_finite() || t < 0.0 || t > year {
                    return false;
                }
            } else if self.completed[i] {
                return false;
            }
        }
        true
    }
}

/// Great Aunt Weena's autonomous industrial and Chamber of Guf planner.
///
/// Strictly causal:
/// 1. Expands industry at Weena systems with `industry < 2.0` when stores permit.
/// 2. Builds synchronized Guf relay nodes at qualified Weena systems.
/// 3. Causally recruits reachable Independent/Drummer systems (`c + 2,000 yr`).
/// 4. Issues finite-slew `aim_guf` commands and fires `fire_guf` ONLY when
///    at least 3 mounts are aligned and `array_charged >= 100.0`.
pub fn plan_weena_expansion(state: &mut StrategicState) -> Option<String> {
    // 1. Ensure Weena systems have enough industry to construct Guf nodes.
    let ind_target = state
        .systems
        .iter()
        .enumerate()
        .find(|(_, s)| {
            s.owner == Faction::Weena
                && s.industry < 2.0
                && s.material >= 45.0
                && s.energy >= 20.0
        })
        .map(|(i, _)| i);
    if let Some(site) = ind_target {
        if let Ok(msg) = state.project(site, Project::Industry) {
            return Some(format!("WEENA AI: {msg}"));
        }
    }

    // 2. Construct Guf relay nodes where qualified.
    let node_target = state
        .systems
        .iter()
        .enumerate()
        .find(|(_, s)| {
            s.owner == Faction::Weena
                && !s.array_node
                && s.industry >= 2.0
                && s.material >= 65.0
                && s.energy >= 75.0
        })
        .map(|(i, _)| i);
    if let Some(site) = node_target {
        if let Ok(msg) = state.project(site, Project::GufNode) {
            return Some(format!("WEENA AI: {msg}"));
        }
    }

    // 3. Recruit nearby Independent or Drummer systems until Weena holds >= 3 worlds.
    let weena_count = state
        .systems
        .iter()
        .filter(|s| s.owner == Faction::Weena)
        .count();
    let claims_in_flight = state
        .signals
        .iter()
        .filter(|p| p.kind == SignalKind::WeenaClaim)
        .count();
    if weena_count + claims_in_flight < 4 && claims_in_flight < 1 {
        let controller_pos = state
            .systems
            .iter()
            .find(|s| s.owner == Faction::Weena)
            .map(|s| s.pos);
        if let Some(origin) = controller_pos {
            let mut candidates: Vec<(u8, f64, usize)> = state
                .systems
                .iter()
                .enumerate()
                .filter(|(idx, s)| {
                    matches!(s.owner, Faction::Independent | Faction::Drummers)
                        && !state
                            .signals
                            .iter()
                            .any(|p| p.kind == SignalKind::WeenaClaim && p.target == *idx)
                })
                .map(|(idx, s)| {
                    // Prefer Independent worlds first, then nearest Drummer systems.
                    let priority = if s.owner == Faction::Independent { 0 } else { 1 };
                    (priority, (s.pos - origin).len(), idx)
                })
                .collect();
            candidates.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
            for (_, _, target) in candidates {
                if let Ok(eta) = state.claim_for_weena(target) {
                    return Some(format!(
                        "WEENA AI: recruitment relay to {} (ETA T+{eta:.0})",
                        state.systems[target].name
                    ));
                }
            }
        }
    }

    // 4. Coordinate Guf array aiming and firing once >= 3 powered relays exist.
    let powered_nodes = state
        .systems
        .iter()
        .filter(|s| s.owner == Faction::Weena && s.array_node && s.energy >= 50.0)
        .count();
    if powered_nodes >= 3 {
        // Choose priority hostile/non-Weena target without an active inbound Guf beam.
        let weena_center = state
            .systems
            .iter()
            .find(|s| s.owner == Faction::Weena && s.array_node)
            .map(|s| s.pos)
            .unwrap_or_default();
        let mut targets: Vec<(u8, f64, usize)> = state
            .systems
            .iter()
            .enumerate()
            .filter(|(idx, s)| {
                s.owner != Faction::Weena && state.incoming_guf(*idx).is_none()
            })
            .map(|(idx, s)| {
                let prio = match s.owner {
                    Faction::Blight => 0,
                    Faction::Drummers => 1,
                    Faction::Independent => 2,
                    Faction::Weena => 3,
                };
                (prio, (s.pos - weena_center).len(), idx)
            })
            .collect();
        targets.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));

        if let Some((_, _, target)) = targets.first().copied() {
            if let Some((assigned, aligned, _)) = state.guf_alignment_status(target) {
                if aligned >= 3 && state.array_charged >= 100.0 {
                    if let Ok(eta) = state.fire_guf(target) {
                        return Some(format!(
                            "WEENA AI: Chamber of Guf fired at {} (ETA T+{eta:.0})",
                            state.systems[target].name
                        ));
                    }
                } else if assigned < 3 {
                    if let Ok(ready) = state.aim_guf(target) {
                        return Some(format!(
                            "WEENA AI: slewing Guf array toward {} (READY T+{ready:.0})",
                            state.systems[target].name
                        ));
                    }
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strategy::{TacticalKind, TacticalOutcome};

    #[test]
    fn episodes_xvii_through_xxii_progress_deterministically() {
        let mut s = StrategicState::new(100_000.0);
        s.blight_autonomy = false;
        assert_eq!(s.episodes.active_index, 0);

        // Ep XVII: expand industry >= 2.5 and fabricate kugelblitz.
        s.systems[5].industry = 2.5;
        s.systems[5].kugelblitz = 1;
        s.evaluate_episodes();
        assert!(s.episodes.completed[0]);
        assert_eq!(s.episodes.active_index, 1);

        // Ep XVIII: win a tactical encounter and broadcast a contact/memetic record.
        s.resolve_tactical(5, TacticalKind::Intercept, TacticalOutcome::Victory)
            .unwrap();
        s.broadcast_memetic_record(5, "ALFONSO CONTACT").unwrap();
        s.evaluate_episodes();
        assert!(s.episodes.completed[1]);
        assert_eq!(s.episodes.active_index, 2);

        // Ep XIX: capture Nebular Wastes (site 6) or defend against Guf fire.
        s.set_owner(6, Faction::Drummers);
        s.evaluate_episodes();
        assert!(s.episodes.completed[2]);
        assert_eq!(s.episodes.active_index, 3);

        // Ep XX: deliver a convoy / hold >= 5 Drummer systems.
        s.set_owner(0, Faction::Drummers);
        s.evaluate_episodes();
        assert!(s.episodes.completed[3]);
        assert_eq!(s.episodes.active_index, 4);

        // Ep XXI: tidal victory.
        s.advance_to(100_300.0);
        s.resolve_tactical(5, TacticalKind::Tidal, TacticalOutcome::Victory)
            .unwrap();
        assert!(s.episodes.completed[4]);
        assert_eq!(s.episodes.active_index, 5);

        // Ep XXII: advance until light-speed memetic reconciliation converges all nodes.
        while let Some(next_t) = s.patricia_trie.next_arrival_after(s.year) {
            s.advance_to(next_t);
        }
        assert_eq!(s.patricia_trie.branch_divergence_count(), 1);
        assert!(s.patricia_trie.average_coexist() >= 0.45);
        assert!(s.episodes.completed[5]);
        assert!(s.episodes.all_complete());
        assert!(s.episodes.valid(s.year));
    }

    #[test]
    fn weena_autonomous_planner_recruits_builds_aims_and_fires_causally() {
        let mut s = StrategicState::new(100_000.0);
        s.blight_autonomy = false;
        s.weena_autonomy = true;
        // Seed 3 Weena systems with plenty of resources so the planner builds nodes,
        // slews mounts causally, charges, and fires only once >= 3 mounts align.
        for i in [2, 3, 6] {
            s.set_owner(i, Faction::Weena);
            s.systems[i].industry = 2.5;
            s.systems[i].material = 300.0;
            s.systems[i].energy = 400.0;
        }
        // First 3 planning ticks build the 3 Guf relay nodes.
        for _ in 0..3 {
            let action = s.plan_weena_expansion();
            assert!(action.is_some());
        }
        assert_eq!(
            s.systems
                .iter()
                .filter(|sys| sys.owner == Faction::Weena && sys.array_node)
                .count(),
            3
        );
        // Next planning tick aims (cannot fire yet because mounts are not aligned!).
        s.array_charged = 100.0;
        let aim_msg = s.plan_weena_expansion().unwrap();
        assert!(aim_msg.contains("slewing") || aim_msg.contains("recruitment"));
        // Even if we call plan again immediately, it must NOT fire while slewing!
        assert_eq!(s.pulses, 0);

        // Advance time until Weena's autonomous 2,500-yr ticks finish recruiting/slewing & firing.
        s.advance_to(135_000.0);
        assert!(s.pulses >= 1);
    }
}
