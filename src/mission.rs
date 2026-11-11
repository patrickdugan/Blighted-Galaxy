//! BG-STR-03: snapshot a star system into a reproducible SI encounter.
//! Galactic logistics are abstract game values; the scene integrators retain
//! their SI units, causal observation histories, and existing collision rules.
use crate::causal::{WorldlineHistory, WorldlineSample};
use crate::physics::{AU_M, C_SI, G_SI, SOLAR_MASS_KG, SiBody};
use crate::si_bridge::{SiJoust, SiTidalEncounter};
use crate::strategy::{Faction, StrategicState, TacticalKind};
use crate::V2;

#[derive(Clone, Copy, Debug)]
pub struct MissionProfile {
    pub site: usize,
    pub site_name: &'static str,
    pub kind: TacticalKind,
    pub owner: Faction,
    pub year: f64,
    pub industry: f64,
    pub defense: f64,
    pub kugelblitz: u32,
    pub enemy_shield_j: f64,
    pub own_shield_j: f64,
    pub enemy_fire_interval_s: f64,
    pub enemy_initial_beta: f64,
    pub intercept_range_m: f64,
    pub hole_mass_solar: f64,
    pub dwarf_mass_solar: f64,
    pub tidal_separation_m: f64,
}

impl MissionProfile {
    /// A single immutable mission seed, derived from the selected system.
    /// All caps are game balance limits, independent of Galactic scale.
    pub fn from_strategy(state: &StrategicState, site: usize,
        kind: TacticalKind) -> Option<Self> {
        let s=state.systems.get(site)?;
        let industry=s.industry.clamp(0.0, 10.0);
        let defense=s.defense.clamp(0.0, 25.0);
        let stock=s.kugelblitz.min(10);
        let hostile=matches!(s.owner,Faction::Blight|Faction::Weena);
        let shield_factor=(0.65+defense*0.12+industry*0.10+
            if hostile {stock as f64*0.07}else{0.0}).clamp(0.5, 4.0);
        let support=if s.owner==Faction::Drummers{
            (stock as f64*0.14+industry*0.035).min(0.8)
        }else{0.0};
        Some(Self{
            site,site_name:s.name,kind,owner:s.owner,year:state.year,
            industry,defense,kugelblitz:stock,
            enemy_shield_j:7.0e17*shield_factor,
            own_shield_j:7.0e17*(1.0+support),
            enemy_fire_interval_s:(3.2-industry*0.15-stock as f64*0.055)
                .clamp(1.15,3.2),
            enemy_initial_beta:(0.008+industry*0.0017).clamp(0.008,0.032),
            intercept_range_m:(410.0e6+defense*12.0e6)
                .clamp(410.0e6,710.0e6),
            // Weena's 12-solar-mass remnant corresponds to Episode XIX.
            // Other encounters use plausible stellar-remnant masses.
            hole_mass_solar:if s.owner==Faction::Weena {12.0}
                else{(8.0+industry*0.55).clamp(8.0,13.5)},
            dwarf_mass_solar:(0.59+defense*0.012).clamp(0.59,0.89),
            tidal_separation_m:(1.30e8+industry*6.0e6+defense*2.0e6)
                .clamp(1.30e8,2.00e8),
        })
    }

    /// Build a joust with initial causal histories matching the new fleet
    /// positions and velocities. This avoids impossible teleport samples.
    pub fn joust(self) -> SiJoust {
        let mut g=SiJoust::new();
        let center=V2::new(AU_M,0.0);
        let offset=V2::new(self.intercept_range_m*0.5,20.0e6);
        g.player.pos_m=center-V2::new(self.intercept_range_m*0.5,20.0e6);
        g.opponent.pos_m=center+offset;
        g.opponent.vel_m_s=V2::new(-self.enemy_initial_beta*C_SI,0.0);
        g.opponent.shield_capacity_j=self.enemy_shield_j;
        g.opponent.shield_energy_j=self.enemy_shield_j;
        g.opponent.cross_section_radius_m=4.0e6*
            (1.0+0.045*self.industry).clamp(1.0,1.5);
        g.player.shield_capacity_j=self.own_shield_j;
        g.player.shield_energy_j=self.own_shield_j;
        g.enemy_fire_interval_s=self.enemy_fire_interval_s;
        g.player_history=WorldlineHistory::new(8192);
        g.opponent_history=WorldlineHistory::new(8192);
        g.player_history.push(WorldlineSample{
            time_s:0.0,pos_m:g.player.pos_m,vel_m_s:g.player.vel_m_s
        }).expect("valid player mission initial history");
        g.opponent_history.push(WorldlineSample{
            time_s:0.0,pos_m:g.opponent.pos_m,vel_m_s:g.opponent.vel_m_s
        }).expect("valid enemy mission initial history");
        g
    }

    /// A black-hole/white-dwarf interdiction lab. The black hole here is a
    /// stellar remnant, distinct from the manufactured micro-kugelblitz
    /// stored in the strategic inventory.
    pub fn tidal(self) -> SiTidalEncounter {
        let mut g=SiTidalEncounter::new();
        let hole_mass=self.hole_mass_solar*SOLAR_MASS_KG;
        let dwarf_mass=self.dwarf_mass_solar*SOLAR_MASS_KG;
        // Non-relativistic radius scaling is an approximate scene choice.
        let dwarf_radius=9.0e6*(0.65/self.dwarf_mass_solar).powf(1.0/3.0);
        g.black_hole=SiBody::new(hole_mass,2.0*G_SI*hole_mass/C_SI.powi(2),
            V2::default(),V2::default());
        g.white_dwarf=SiBody::new(dwarf_mass,dwarf_radius,
            V2::new(-self.tidal_separation_m,2.0e7),
            V2::new(6.0e5,-8.0e4));
        g
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn disparate_systems_seed_distinct_shields_and_hole_masses() {
        let mut state=StrategicState::default();
        let low=MissionProfile::from_strategy(&state,0,TacticalKind::Intercept).unwrap();
        state.systems[0].industry=9.0;
        state.systems[0].defense=14.0;
        let high=MissionProfile::from_strategy(&state,0,TacticalKind::Intercept).unwrap();
        assert!(high.enemy_shield_j>low.enemy_shield_j);
        assert!(high.enemy_fire_interval_s<low.enemy_fire_interval_s);
        assert!(high.hole_mass_solar>low.hole_mass_solar);
        assert!(high.intercept_range_m>low.intercept_range_m);
    }
    #[test]
    fn stored_kugelblitz_supports_homeworld_defenses_without_consuming_stock() {
        let mut state=StrategicState::default();
        state.systems[5].kugelblitz=4;
        let profile=MissionProfile::from_strategy(&state,5,TacticalKind::Intercept).unwrap();
        let battle=profile.joust();
        assert!(battle.player.shield_capacity_j>7.0e17);
        assert_eq!(state.systems[5].kugelblitz,4);
    }
    #[test]
    fn mission_joust_history_starts_at_actual_scene_positions() {
        let state=StrategicState::default();
        let profile=MissionProfile::from_strategy(&state,1,TacticalKind::Intercept).unwrap();
        let mut g=profile.joust();
        assert!(g.player_history.validate_segments());
        assert!(g.opponent_history.validate_segments());
        for _ in 0..400 {
            g.step(crate::si_bridge::Seconds(1.0/120.0),
                crate::si_bridge::Control::default());
        }
        assert!(g.player_history.validate_segments());
        assert!(g.opponent_history.validate_segments());
        assert!(g.player.vel_m_s.len()<C_SI);
    }
    #[test]
    fn stellar_remnant_remains_distinct_from_strategic_micro_hole() {
        let mut state=StrategicState::default();
        state.systems[6].kugelblitz=3;
        let profile=MissionProfile::from_strategy(&state,6,TacticalKind::Tidal).unwrap();
        let scene=profile.tidal();
        assert!((scene.black_hole.mass_kg/SOLAR_MASS_KG-12.0).abs()<1e-10);
        assert!(scene.white_dwarf.mass_kg>0.0);
        assert_eq!(state.systems[6].kugelblitz,3);
    }
    #[test]
    fn scene_snapshot_is_stable_after_system_development() {
        let mut state=StrategicState::default();
        let original=MissionProfile::from_strategy(&state,5,TacticalKind::Intercept).unwrap();
        state.systems[5].industry=10.0;
        let next=MissionProfile::from_strategy(&state,5,TacticalKind::Intercept).unwrap();
        assert!(next.enemy_shield_j>original.enemy_shield_j);
        assert_eq!(original.joust().opponent.shield_energy_j,original.enemy_shield_j);
    }
}
