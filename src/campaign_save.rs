//! BG-STR-05 / BG-STR-06: versioned JSON campaign persistence (schema v3 + v2 migration).
//! Save the authoritative model state, not GPU textures or derived map pixels.
//! All untrusted save content is checked before it is admitted into gameplay.
use crate::episodes::EpisodeTracker;
use crate::guf_geometry::SlewMount;
use crate::overworld::{GalacticNavModel, TravelLeg, TIMELINE_MAX_YR};
use crate::patricia_trie::PatriciaTrieWorld;
use crate::strategy::{CargoConvoy, Faction, Signal, SignalKind, StrategicState, System};
use crate::V2;
use serde::{Deserialize, Serialize};
use std::{error::Error, fmt, fs, path::Path};

pub const SAVE_FORMAT: &str = "BLIGHTED-GALAXY-CAMPAIGN";
pub const SAVE_VERSION: u32 = 3;
const LEGACY_V2_VERSION: u32 = 2;
const MAX_PACKETS: usize = 4096;
const MAX_EVENTS: usize = 256;

#[derive(Debug)]
pub enum SaveError {
    Io(std::io::Error),
    Json(serde_json::Error),
    Invalid(String),
}
impl fmt::Display for SaveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "save I/O: {e}"),
            Self::Json(e) => write!(f, "save JSON: {e}"),
            Self::Invalid(e) => write!(f, "save invalid: {e}"),
        }
    }
}
impl Error for SaveError {}
impl From<std::io::Error> for SaveError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
impl From<serde_json::Error> for SaveError {
    fn from(e: serde_json::Error) -> Self {
        Self::Json(e)
    }
}
fn bad(msg: &str) -> SaveError {
    SaveError::Invalid(msg.into())
}
fn finite(v: f64) -> bool {
    v.is_finite()
}
fn world(p: V2) -> bool {
    finite(p.x) && finite(p.y) && p.len() < 1.0e8
}
fn interval(v: f64, lo: f64, hi: f64) -> bool {
    finite(v) && v >= lo && v <= hi
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SystemSave {
    owner: Faction,
    industry: f64,
    energy: f64,
    material: f64,
    kugelblitz: u32,
    array_node: bool,
    defense: f64,
    guf_ward: f64,
}
impl From<&System> for SystemSave {
    fn from(s: &System) -> Self {
        Self {
            owner: s.owner,
            industry: s.industry,
            energy: s.energy,
            material: s.material,
            kugelblitz: s.kugelblitz,
            array_node: s.array_node,
            defense: s.defense,
            guf_ward: s.guf_ward,
        }
    }
}
impl SystemSave {
    fn apply(self, site: &mut System) {
        site.owner = self.owner;
        site.industry = self.industry;
        site.energy = self.energy;
        site.material = self.material;
        site.kugelblitz = self.kugelblitz;
        site.array_node = self.array_node;
        site.defense = self.defense;
        site.guf_ward = self.guf_ward;
    }
    fn valid(&self) -> bool {
        interval(self.industry, 0.0, 1000.0)
            && interval(self.energy, 0.0, 1000.0)
            && interval(self.material, 0.0, 600.0)
            && self.kugelblitz <= 1_000_000
            && interval(self.defense, -100_000.0, 100_000.0)
            && interval(self.guf_ward, 0.0, 1_000_000.0)
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CargoConvoyV2 {
    id: u64,
    source: usize,
    target: usize,
    owner: Faction,
    material: f64,
    energy: f64,
    launch_year: f64,
    arrival_year: f64,
    beta: f64,
}
impl From<CargoConvoyV2> for CargoConvoy {
    fn from(c: CargoConvoyV2) -> Self {
        Self {
            id: c.id,
            source: c.source,
            target: c.target,
            owner: c.owner,
            material: c.material,
            energy: c.energy,
            launch_year: c.launch_year,
            arrival_year: c.arrival_year,
            beta: c.beta,
            escorted: false,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StrategySaveV2 {
    year: f64,
    systems: Vec<SystemSave>,
    signals: Vec<Signal>,
    convoys: Vec<CargoConvoyV2>,
    next_convoy_id: u64,
    events: Vec<String>,
    array_charged: f64,
    pulses: u32,
    historical: bool,
    blight_autonomy: bool,
    guf_mounts: Vec<SlewMount>,
    territory_revision: u64,
    last_tactical_year: Vec<Option<f64>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CampaignFileV2 {
    format: String,
    version: u32,
    site_fingerprint: u64,
    nav: GalacticNavModel,
    strategy: StrategySaveV2,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StrategySaveV3 {
    year: f64,
    systems: Vec<SystemSave>,
    signals: Vec<Signal>,
    convoys: Vec<CargoConvoy>,
    next_convoy_id: u64,
    events: Vec<String>,
    array_charged: f64,
    pulses: u32,
    historical: bool,
    blight_autonomy: bool,
    weena_autonomy: bool,
    episodes: EpisodeTracker,
    patricia_trie: PatriciaTrieWorld,
    guf_mounts: Vec<SlewMount>,
    territory_revision: u64,
    last_tactical_year: Vec<Option<f64>>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CampaignFileV3 {
    format: String,
    version: u32,
    site_fingerprint: u64,
    nav: GalacticNavModel,
    strategy: StrategySaveV3,
}

/// A stable fingerprint of canonical waypoint identifiers and geometry.
/// Versioned save slots deliberately reject incompatible maps.
fn fingerprint(sites: &[System]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for site in sites {
        for byte in site
            .name
            .as_bytes()
            .iter()
            .copied()
            .chain(site.pos.x.to_bits().to_le_bytes())
            .chain(site.pos.y.to_bits().to_le_bytes())
        {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}

fn capture(nav: &GalacticNavModel, strategy: &StrategicState) -> CampaignFileV3 {
    CampaignFileV3 {
        format: SAVE_FORMAT.into(),
        version: SAVE_VERSION,
        site_fingerprint: fingerprint(&strategy.systems),
        nav: nav.clone(),
        strategy: StrategySaveV3 {
            year: strategy.year,
            systems: strategy.systems.iter().map(SystemSave::from).collect(),
            signals: strategy.signals.clone(),
            convoys: strategy.convoys.clone(),
            next_convoy_id: strategy.next_convoy_id,
            events: strategy.events.clone(),
            array_charged: strategy.array_charged,
            pulses: strategy.pulses,
            historical: strategy.historical,
            blight_autonomy: strategy.blight_autonomy,
            weena_autonomy: strategy.weena_autonomy,
            episodes: strategy.episodes.clone(),
            patricia_trie: strategy.patricia_trie.clone(),
            guf_mounts: strategy.guf_mounts.clone(),
            territory_revision: strategy.territory_revision,
            last_tactical_year: strategy
                .last_tactical_year
                .iter()
                .map(|&t| if t == f64::NEG_INFINITY { None } else { Some(t) })
                .collect(),
        },
    }
}

fn valid_leg(t: TravelLeg) -> bool {
    world(t.start)
        && world(t.destination)
        && interval(t.cruise_beta, 0.001, 0.999999)
        && interval(t.total_distance_ly, 1.0, 200_000.0)
        && interval(t.distance_travelled_ly, 0.0, t.total_distance_ly)
        && interval(t.proper_elapsed_yr, 0.0, TIMELINE_MAX_YR)
        && (t.total_distance_ly - (t.destination - t.start).len()).abs() < 1e-6
}

fn validate_nav(nav: &GalacticNavModel) -> Result<(), SaveError> {
    if !interval(nav.year, 0.0, TIMELINE_MAX_YR)
        || !interval(nav.rate_yr_per_second, 0.01, 100_000.0)
        || !interval(nav.proper_time_total_yr, 0.0, 1e9)
        || !world(nav.position)
        || nav.destination.is_some_and(|p| !world(p))
        || nav.last_arrival.is_some_and(|p| !world(p))
    {
        return Err(bad("invalid navigation clock, position or rate"));
    }
    if let Some(t) = nav.travel {
        if !valid_leg(t) || (nav.position - t.position()).len() > 1e-4 {
            return Err(bad("invalid in-flight relativistic leg"));
        }
    }
    Ok(())
}

fn validate_packet(p: &Signal, systems: &[System], year: f64) -> bool {
    if p.source >= systems.len()
        || p.target >= systems.len()
        || !interval(p.launch_year, 0.0, TIMELINE_MAX_YR)
        || !interval(p.arrival_year, p.launch_year, TIMELINE_MAX_YR)
        || p.arrival_year < year
        || !interval(p.power, 0.0, 1e7)
    {
        return false;
    }
    let length = (systems[p.source].pos - systems[p.target].pos).len();
    let duration = p.arrival_year - p.launch_year;
    let limit = match p.kind {
        SignalKind::ArrayPulse | SignalKind::WeenaClaim => length,
        SignalKind::BlightSeed => length / 0.10,
        SignalKind::HaloStrike => length / 0.76,
    };
    let aim_valid = match p.kind {
        SignalKind::ArrayPulse => p.beam_azimuth_rad.is_some_and(|az| finite(az) && az.abs() < 100.0),
        _ => p.beam_azimuth_rad.is_none(),
    };
    aim_valid && duration + 1e-6 >= limit
}

fn validate_convoy(c: &CargoConvoy, systems: &[System], year: f64) -> bool {
    if c.source >= systems.len()
        || c.target >= systems.len()
        || c.source == c.target
        || c.id == 0
        || !interval(c.launch_year, 0.0, year)
        || !interval(c.arrival_year, year, TIMELINE_MAX_YR)
        || !interval(c.beta, 0.01, 0.999999)
        || !interval(c.material, 0.0, 600.0)
        || !interval(c.energy, 0.0, 1000.0)
    {
        return false;
    }
    let duration = (systems[c.source].pos - systems[c.target].pos).len() / c.beta;
    (c.arrival_year - c.launch_year - duration).abs() <= 1e-6
}

fn validate_and_apply_common(
    site_fingerprint: u64,
    nav: GalacticNavModel,
    year: f64,
    systems_save: Vec<SystemSave>,
    signals: Vec<Signal>,
    convoys: Vec<CargoConvoy>,
    next_convoy_id: u64,
    events: Vec<String>,
    array_charged: f64,
    pulses: u32,
    historical: bool,
    blight_autonomy: bool,
    guf_mounts: Vec<SlewMount>,
    territory_revision: u64,
    last_tactical_year: Vec<Option<f64>>,
) -> Result<(GalacticNavModel, StrategicState), SaveError> {
    validate_nav(&nav)?;
    if !interval(year, 0.0, TIMELINE_MAX_YR) || (nav.year - year).abs() > 1e-7 {
        return Err(bad("mismatched galactic clocks"));
    }
    let mut state = StrategicState::new(year);
    if site_fingerprint != fingerprint(&state.systems) {
        return Err(bad("incompatible canonical star map"));
    }
    if systems_save.len() != state.systems.len()
        || last_tactical_year.len() != state.systems.len()
        || guf_mounts.len() != state.systems.len()
    {
        return Err(bad("incorrect number of industrial systems"));
    }
    if signals.len() > MAX_PACKETS
        || convoys.len() > MAX_PACKETS
        || events.len() > MAX_EVENTS
        || events.iter().any(|e| e.len() > 512)
    {
        return Err(bad("oversized packets, convoy queue or event log"));
    }
    if !interval(array_charged, 0.0, 100.0) || next_convoy_id == 0 {
        return Err(bad("invalid Guf charge or convoy sequence"));
    }
    for (saved, site) in systems_save.into_iter().zip(state.systems.iter_mut()) {
        if !saved.valid() {
            return Err(bad("invalid industrial resource balance"));
        }
        saved.apply(site);
    }
    if !signals
        .iter()
        .all(|p| validate_packet(p, &state.systems, state.year))
    {
        return Err(bad("invalid attack route or light-cone time"));
    }
    if !convoys
        .iter()
        .all(|c| validate_convoy(c, &state.systems, state.year))
    {
        return Err(bad("invalid convoy route or sublight time"));
    }
    let newest = convoys.iter().map(|c| c.id).max().unwrap_or(0);
    if next_convoy_id <= newest {
        return Err(bad("convoy sequence may duplicate existing cargo"));
    }
    if !guf_mounts
        .iter()
        .all(|m| m.valid(state.year, state.systems.len()))
    {
        return Err(bad("Guf mount targeting or slew clock invalid"));
    }
    if !last_tactical_year
        .iter()
        .all(|t| t.is_none_or(|v| interval(v, 0.0, state.year)))
    {
        return Err(bad("tactical cooldown timestamp is invalid"));
    }
    state.signals = signals;
    state.convoys = convoys;
    state.next_convoy_id = next_convoy_id;
    state.events = events;
    state.array_charged = array_charged;
    state.pulses = pulses;
    state.historical = historical;
    state.blight_autonomy = blight_autonomy;
    state.territory_revision = territory_revision;
    state.guf_mounts = guf_mounts;
    state.last_tactical_year = last_tactical_year
        .into_iter()
        .map(|x| x.unwrap_or(f64::NEG_INFINITY))
        .collect();
    Ok((nav, state))
}

pub fn encode(nav: &GalacticNavModel, strategy: &StrategicState) -> Result<String, SaveError> {
    validate_nav(nav)?;
    if (nav.year - strategy.year).abs() > 1e-7 {
        return Err(bad("navigation and strategic clocks differ"));
    }
    let raw = serde_json::to_string_pretty(&capture(nav, strategy))?;
    // Validate the emitted save too, including its finite fields.
    let _ = decode(&raw)?;
    Ok(raw)
}

pub fn decode(raw: &str) -> Result<(GalacticNavModel, StrategicState), SaveError> {
    if raw.len() > 8 * 1024 * 1024 {
        return Err(bad("oversized campaign file"));
    }
    let value: serde_json::Value = serde_json::from_str(raw)?;
    let format = value
        .get("format")
        .and_then(|v| v.as_str())
        .ok_or_else(|| bad("missing save format header"))?;
    let version = value
        .get("version")
        .and_then(|v| v.as_u64())
        .and_then(|v| u32::try_from(v).ok())
        .ok_or_else(|| bad("missing or invalid save version"))?;
    if format != SAVE_FORMAT {
        return Err(bad("unsupported save format/version"));
    }

    match version {
        SAVE_VERSION => {
            let file: CampaignFileV3 = serde_json::from_value(value)?;
            if file.format != SAVE_FORMAT || file.version != SAVE_VERSION {
                return Err(bad("unsupported save format/version"));
            }
            let d = file.strategy;
            let (nav, mut state) = validate_and_apply_common(
                file.site_fingerprint,
                file.nav,
                d.year,
                d.systems,
                d.signals,
                d.convoys,
                d.next_convoy_id,
                d.events,
                d.array_charged,
                d.pulses,
                d.historical,
                d.blight_autonomy,
                d.guf_mounts,
                d.territory_revision,
                d.last_tactical_year,
            )?;
            if !d.episodes.valid(state.year) {
                return Err(bad("invalid episode progression state"));
            }
            if !d.patricia_trie.valid(&state.systems, state.year) {
                return Err(bad("invalid Patricia Trie root, policy or causal queue"));
            }
            state.weena_autonomy = d.weena_autonomy;
            state.episodes = d.episodes;
            state.patricia_trie = d.patricia_trie;
            Ok((nav, state))
        }
        LEGACY_V2_VERSION => {
            let file: CampaignFileV2 = serde_json::from_value(value)?;
            if file.format != SAVE_FORMAT || file.version != LEGACY_V2_VERSION {
                return Err(bad("unsupported save format/version"));
            }
            let d = file.strategy;
            let convoys: Vec<CargoConvoy> = d.convoys.into_iter().map(CargoConvoy::from).collect();
            let (nav, mut state) = validate_and_apply_common(
                file.site_fingerprint,
                file.nav,
                d.year,
                d.systems,
                d.signals,
                convoys,
                d.next_convoy_id,
                d.events,
                d.array_charged,
                d.pulses,
                d.historical,
                d.blight_autonomy,
                d.guf_mounts,
                d.territory_revision,
                d.last_tactical_year,
            )?;
            // Clean v2 -> v3 migration defaults.
            state.weena_autonomy = false;
            state.episodes = EpisodeTracker::default();
            state.patricia_trie = PatriciaTrieWorld::new_for_systems(&state.systems, state.year);
            Ok((nav, state))
        }
        _ => Err(bad("unsupported save format/version")),
    }
}

pub fn save_file(
    path: &Path,
    nav: &GalacticNavModel,
    strategy: &StrategicState,
) -> Result<(), SaveError> {
    let text = encode(nav, strategy)?;
    let tmp = path.with_extension("json.tmp");
    let backup = path.with_extension("json.bak");
    fs::write(&tmp, text.as_bytes())?;
    // Keep the previous valid slot as a recovery point on Windows too.
    if path.exists() {
        if backup.exists() {
            fs::remove_file(&backup)?;
        }
        fs::rename(path, &backup)?;
    }
    if let Err(err) = fs::rename(&tmp, path) {
        if backup.exists() {
            let _ = fs::rename(&backup, path);
        }
        return Err(SaveError::Io(err));
    }
    Ok(())
}

pub fn load_file(path: &Path) -> Result<(GalacticNavModel, StrategicState), SaveError> {
    let data = fs::read_to_string(path)?;
    decode(&data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strategy::{TacticalKind, TacticalOutcome};

    #[test]
    fn round_trip_preserves_resources_in_flight_convoys_guf_and_cooldowns() {
        let mut nav = GalacticNavModel::default();
        let mut s = StrategicState::default();
        s.blight_autonomy = false;
        s.weena_autonomy = true;
        s.systems[3].material = 350.0;
        s.systems[3].energy = 500.0;
        s.send_convoy(3, 5).unwrap();
        s.resolve_tactical(3, TacticalKind::ConvoyEscort, TacticalOutcome::Victory)
            .unwrap();
        s.broadcast_memetic_record(5, "TEST BROADCAST").unwrap();
        for i in [2, 3, 6] {
            s.set_owner(i, Faction::Weena);
            s.systems[i].array_node = true;
            s.systems[i].energy = 500.0;
        }
        s.blight_autonomy = false;
        s.weena_autonomy = false;
        let ready = s.aim_guf(5).unwrap();
        s.advance_to(ready);
        s.array_charged = 100.0;
        s.fire_guf(5).unwrap();
        s.systems[5].energy = 500.0;
        s.resolve_tactical(5, TacticalKind::GufDefense, TacticalOutcome::Victory)
            .unwrap();
        nav.seek(s.year);
        nav.choose(crate::overworld::ORION_OUTBREAK);
        assert!(nav.depart(0.76));
        nav.tick(0.05);
        s.advance_to(nav.year);
        let original = encode(&nav, &s).unwrap();
        let (n, t) = decode(&original).unwrap();
        let restored = encode(&n, &t).unwrap();
        assert_eq!(original, restored);
        assert_eq!(t.convoys.len(), 1);
        assert!(t.convoys[0].escorted);
        assert_eq!(t.incoming_guf(5).unwrap().1, 3);
        assert_eq!(t.last_tactical_year[5], s.last_tactical_year[5]);
        assert_eq!(t.patricia_trie, s.patricia_trie);
        assert_eq!(t.episodes, s.episodes);
        assert!(n.travel.is_some());
    }

    #[test]
    fn save_rejects_forged_resource_and_bad_version() {
        let nav = GalacticNavModel::default();
        let s = StrategicState::default();
        let data = encode(&nav, &s).unwrap();
        assert!(decode(&data.replace("\"version\": 3", "\"version\": 1")).is_err());
        assert!(decode(&data.replace("\"version\": 3", "\"version\": 99")).is_err());
        assert!(decode(&data.replace("\"energy\": 100.0", "\"energy\": -100.0")).is_err());
    }

    #[test]
    fn v2_save_migrates_cleanly_to_v3_defaults() {
        let nav = GalacticNavModel::default();
        let mut s = StrategicState::default();
        s.blight_autonomy = false;
        s.systems[3].material = 300.0;
        s.systems[3].energy = 300.0;
        s.send_convoy(3, 5).unwrap();

        // Construct an authentic schema-v2 JSON payload (no weena_autonomy, episodes,
        // patricia_trie, or convoy.escorted).
        let v3_val: serde_json::Value = serde_json::from_str(&encode(&nav, &s).unwrap()).unwrap();
        let mut v2 = v3_val.clone();
        v2["version"] = serde_json::json!(2);
        let strat = v2["strategy"].as_object_mut().unwrap();
        strat.remove("weena_autonomy");
        strat.remove("episodes");
        strat.remove("patricia_trie");
        for convoy in strat["convoys"].as_array_mut().unwrap() {
            convoy.as_object_mut().unwrap().remove("escorted");
        }
        let v2_json = serde_json::to_string_pretty(&v2).unwrap();

        let (migrated_nav, migrated_state) = decode(&v2_json).unwrap();
        assert_eq!(migrated_nav.year, nav.year);
        assert!(!migrated_state.weena_autonomy);
        assert_eq!(migrated_state.episodes.active_index, 0);
        assert_eq!(migrated_state.convoys.len(), 1);
        assert!(!migrated_state.convoys[0].escorted);
        assert!(migrated_state
            .patricia_trie
            .valid(&migrated_state.systems, migrated_state.year));
    }

    #[test]
    fn v3_rejects_tampered_trie_root_policy_and_episode_bounds() {
        let nav = GalacticNavModel::default();
        let s = StrategicState::default();
        let base_json = encode(&nav, &s).unwrap();

        // Tamper trie root_hash.
        let mut bad_root: serde_json::Value = serde_json::from_str(&base_json).unwrap();
        bad_root["strategy"]["patricia_trie"]["nodes"][0]["root_hash"] = serde_json::json!(12345);
        assert!(decode(&serde_json::to_string(&bad_root).unwrap()).is_err());

        // Tamper policy weights out of sum=1.0 bounds.
        let mut bad_policy: serde_json::Value = serde_json::from_str(&base_json).unwrap();
        bad_policy["strategy"]["patricia_trie"]["nodes"][0]["policy"]["coexist"] =
            serde_json::json!(0.99);
        assert!(decode(&serde_json::to_string(&bad_policy).unwrap()).is_err());

        // Tamper episode active_index past EPISODE_COUNT.
        let mut bad_ep: serde_json::Value = serde_json::from_str(&base_json).unwrap();
        bad_ep["strategy"]["episodes"]["active_index"] = serde_json::json!(99);
        assert!(decode(&serde_json::to_string(&bad_ep).unwrap()).is_err());
    }

    #[test]
    fn saving_before_any_tactical_engagement_handles_negative_infinity() {
        let nav = GalacticNavModel::default();
        let s = StrategicState::default();
        let raw = encode(&nav, &s).unwrap();
        assert!(raw.contains("\"last_tactical_year\""));
        let (_, r) = decode(&raw).unwrap();
        assert_eq!(r.last_tactical_year[0], f64::NEG_INFINITY);
    }

    #[test]
    fn in_progress_guf_slew_replays_after_loading() {
        let mut nav = GalacticNavModel::default();
        let mut s = StrategicState::default();
        s.blight_autonomy = false;
        for i in [2, 3, 6] {
            s.set_owner(i, Faction::Weena);
            s.systems[i].array_node = true;
            s.systems[i].energy = 500.0;
        }
        let target = 5;
        let ready = s.aim_guf(target).unwrap();
        s.advance_to(100_250.0);
        nav.seek(100_250.0);
        let (mut n, mut restored) = decode(&encode(&nav, &s).unwrap()).unwrap();
        assert!(restored
            .guf_mounts
            .iter()
            .any(|m| m.assigned_target == Some(target) && m.slew_end_year > restored.year));
        s.advance_to(ready);
        restored.advance_to(ready);
        nav.seek(ready);
        n.seek(ready);
        s.array_charged = 100.0;
        restored.array_charged = 100.0;
        assert_eq!(
            s.fire_guf(target).unwrap(),
            restored.fire_guf(target).unwrap()
        );
        assert_eq!(encode(&nav, &s).unwrap(), encode(&n, &restored).unwrap());
    }

    #[test]
    fn file_save_keeps_previous_version_in_backup() {
        let path = std::env::temp_dir().join(format!(
            "bg-str-06-{}-{}.json",
            std::process::id(),
            SAVE_VERSION
        ));
        let backup = path.with_extension("json.bak");
        let temp = path.with_extension("json.tmp");
        for p in [&path, &backup, &temp] {
            let _ = fs::remove_file(p);
        }
        let mut nav = GalacticNavModel::default();
        let mut s = StrategicState::default();
        s.blight_autonomy = false;
        save_file(&path, &nav, &s).unwrap();
        nav.seek(100_500.0);
        s.advance_to(nav.year);
        save_file(&path, &nav, &s).unwrap();
        let (current, _) = load_file(&path).unwrap();
        let (previous, _) = load_file(&backup).unwrap();
        assert_eq!(current.year, 100_500.0);
        assert_eq!(previous.year, 100_000.0);
        assert!(!temp.exists());
        for p in [&path, &backup, &temp] {
            let _ = fs::remove_file(p);
        }
    }

    #[test]
    fn convoy_id_cannot_be_replayed_by_file_edit() {
        let nav = GalacticNavModel::default();
        let mut s = StrategicState::default();
        s.systems[3].material = 400.0;
        s.systems[3].energy = 400.0;
        s.send_convoy(3, 5).unwrap();
        s.next_convoy_id = 1;
        assert!(encode(&nav, &s).is_err());
    }
}
