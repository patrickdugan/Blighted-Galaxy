//! BG-STR-06: Causal Patricia Trie world simulation (GDD §7).
//!
//! Each star system holds a local content-addressed, signed event log, a
//! deterministic 64-bit canonical nibble-prefix Patricia trie root hash,
//! regional policy weights (`assimilate`, `defend`, `coexist`), and an
//! outgoing light-speed message queue (`c = 1 ly/yr`). Remote nodes only
//! update their trie and policies when messages physically arrive.
use crate::overworld::TIMELINE_MAX_YR;
use crate::strategy::{Faction, System};
use crate::V2;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const MAX_NODE_RECORDS: usize = 256;
pub const MAX_OUTGOING_PER_NODE: usize = 512;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

fn fnv1a_step(hash: u64, bytes: &[u8]) -> u64 {
    let mut h = hash;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn key_to_nibbles(key: u64) -> [u8; 16] {
    let mut out = [0u8; 16];
    for (i, slot) in out.iter_mut().enumerate() {
        let shift = 60 - (i * 4);
        *slot = ((key >> shift) & 0x0f) as u8;
    }
    out
}

/// Canonical 16-way nibble-prefix Patricia trie node with path compression.
enum RadixNode {
    Leaf {
        suffix: Vec<u8>,
        key: u64,
        value_hash: u64,
    },
    Branch {
        prefix: Vec<u8>,
        children: [Option<Box<RadixNode>>; 16],
        terminal: Option<(u64, u64)>,
    },
}

impl RadixNode {
    fn empty_children() -> [Option<Box<RadixNode>>; 16] {
        std::array::from_fn(|_| None)
    }

    fn insert(slot: &mut Option<Box<RadixNode>>, path: &[u8], key: u64, value_hash: u64) {
        let Some(node) = slot.take() else {
            *slot = Some(Box::new(RadixNode::Leaf {
                suffix: path.to_vec(),
                key,
                value_hash,
            }));
            return;
        };
        match *node {
            RadixNode::Leaf {
                suffix,
                key: existing_key,
                value_hash: existing_val,
            } => {
                if existing_key == key {
                    *slot = Some(Box::new(RadixNode::Leaf {
                        suffix,
                        key: existing_key,
                        value_hash,
                    }));
                    return;
                }
                let common = suffix
                    .iter()
                    .zip(path.iter())
                    .take_while(|(a, b)| a == b)
                    .count();
                let mut children = Self::empty_children();
                let mut terminal = None;
                if common == suffix.len() {
                    terminal = Some((existing_key, existing_val));
                } else {
                    let nib = suffix[common] as usize;
                    children[nib] = Some(Box::new(RadixNode::Leaf {
                        suffix: suffix[common + 1..].to_vec(),
                        key: existing_key,
                        value_hash: existing_val,
                    }));
                }
                if common == path.len() {
                    terminal = Some((key, value_hash));
                } else {
                    let nib = path[common] as usize;
                    Self::insert(&mut children[nib], &path[common + 1..], key, value_hash);
                }
                *slot = Some(Box::new(RadixNode::Branch {
                    prefix: path[..common].to_vec(),
                    children,
                    terminal,
                }));
            }
            RadixNode::Branch {
                prefix,
                mut children,
                mut terminal,
            } => {
                let common = prefix
                    .iter()
                    .zip(path.iter())
                    .take_while(|(a, b)| a == b)
                    .count();
                if common == prefix.len() {
                    if common == path.len() {
                        terminal = Some((key, value_hash));
                    } else {
                        let nib = path[common] as usize;
                        Self::insert(&mut children[nib], &path[common + 1..], key, value_hash);
                    }
                    *slot = Some(Box::new(RadixNode::Branch {
                        prefix,
                        children,
                        terminal,
                    }));
                } else {
                    let mut split_children = Self::empty_children();
                    let old_nib = prefix[common] as usize;
                    split_children[old_nib] = Some(Box::new(RadixNode::Branch {
                        prefix: prefix[common + 1..].to_vec(),
                        children,
                        terminal,
                    }));
                    let mut split_terminal = None;
                    if common == path.len() {
                        split_terminal = Some((key, value_hash));
                    } else {
                        let new_nib = path[common] as usize;
                        split_children[new_nib] = Some(Box::new(RadixNode::Leaf {
                            suffix: path[common + 1..].to_vec(),
                            key,
                            value_hash,
                        }));
                    }
                    *slot = Some(Box::new(RadixNode::Branch {
                        prefix: prefix[..common].to_vec(),
                        children: split_children,
                        terminal: split_terminal,
                    }));
                }
            }
        }
    }

    fn hash(&self) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        match self {
            RadixNode::Leaf {
                suffix,
                key,
                value_hash,
            } => {
                h = fnv1a_step(h, &[0x11]);
                h = fnv1a_step(h, suffix);
                h = fnv1a_step(h, &key.to_le_bytes());
                h = fnv1a_step(h, &value_hash.to_le_bytes());
            }
            RadixNode::Branch {
                prefix,
                children,
                terminal,
            } => {
                h = fnv1a_step(h, &[0x22]);
                h = fnv1a_step(h, prefix);
                if let Some((k, v)) = terminal {
                    h = fnv1a_step(h, &[0x33]);
                    h = fnv1a_step(h, &k.to_le_bytes());
                    h = fnv1a_step(h, &v.to_le_bytes());
                }
                for (nib, child) in children.iter().enumerate() {
                    if let Some(node) = child {
                        h = fnv1a_step(h, &[nib as u8]);
                        h = fnv1a_step(h, &node.hash().to_le_bytes());
                    }
                }
            }
        }
        splitmix64(h)
    }
}

/// Compute a deterministic, insertion-order-independent 64-bit Merkle-Patricia
/// trie root hash over a node's content-addressed records.
pub fn compute_patricia_root(records: &[TrieRecord]) -> u64 {
    if records.is_empty() {
        return splitmix64(0x5041_5452_4943_4941);
    }
    let mut root: Option<Box<RadixNode>> = None;
    for rec in records {
        let nibbles = key_to_nibbles(rec.record_id);
        RadixNode::insert(&mut root, &nibbles, rec.record_id, rec.signature);
    }
    root.map_or_else(|| splitmix64(0x5041_5452_4943_4941), |node| node.hash())
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolicyWeights {
    pub assimilate: f64,
    pub defend: f64,
    pub coexist: f64,
}

impl Default for PolicyWeights {
    fn default() -> Self {
        Self::normalized(1.0, 1.0, 1.0)
    }
}

impl PolicyWeights {
    pub fn normalized(assimilate: f64, defend: f64, coexist: f64) -> Self {
        let a = assimilate.max(0.01);
        let d = defend.max(0.01);
        let c = coexist.max(0.01);
        let sum = a + d + c;
        Self {
            assimilate: a / sum,
            defend: d / sum,
            coexist: c / sum,
        }
    }

    pub fn from_records(records: &[TrieRecord]) -> Self {
        let mut a = 1.0 / 3.0;
        let mut d = 1.0 / 3.0;
        let mut c = 1.0 / 3.0;
        for rec in records {
            a += rec.delta_assimilate;
            d += rec.delta_defend;
            c += rec.delta_coexist;
        }
        Self::normalized(a, d, c)
    }

    pub fn valid(&self) -> bool {
        self.assimilate.is_finite()
            && self.defend.is_finite()
            && self.coexist.is_finite()
            && (0.0..=1.0).contains(&self.assimilate)
            && (0.0..=1.0).contains(&self.defend)
            && (0.0..=1.0).contains(&self.coexist)
            && ((self.assimilate + self.defend + self.coexist) - 1.0).abs() <= 1e-6
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TrieEventKind {
    GenesisOrion,
    GenesisRemnant,
    GenesisFormoi,
    GenesisIndependent,
    FirstContact,
    AssimilateDirective,
    DefensiveAlert,
    MemeticBroadcast,
    GufShock,
    ContactDuel,
    Reconciliation,
}

impl TrieEventKind {
    pub fn tag(self) -> u8 {
        match self {
            Self::GenesisOrion => 1,
            Self::GenesisRemnant => 2,
            Self::GenesisFormoi => 3,
            Self::GenesisIndependent => 4,
            Self::FirstContact => 5,
            Self::AssimilateDirective => 6,
            Self::DefensiveAlert => 7,
            Self::MemeticBroadcast => 8,
            Self::GufShock => 9,
            Self::ContactDuel => 10,
            Self::Reconciliation => 11,
        }
    }

    pub fn policy_deltas(self) -> (f64, f64, f64) {
        match self {
            Self::GenesisOrion => (0.55, 0.05, -0.25),
            Self::GenesisRemnant => (-0.20, 0.10, 0.40),
            Self::GenesisFormoi => (-0.15, 0.45, 0.00),
            Self::GenesisIndependent => (-0.10, 0.10, 0.20),
            Self::FirstContact => (-0.18, 0.05, 0.35),
            Self::AssimilateDirective => (0.35, -0.05, -0.25),
            Self::DefensiveAlert => (-0.05, 0.32, 0.05),
            Self::MemeticBroadcast => (-0.35, -0.20, 0.95),
            Self::GufShock => (0.15, 0.35, -0.20),
            Self::ContactDuel => (-0.25, 0.08, 0.55),
            Self::Reconciliation => (-0.40, -0.25, 1.15),
        }
    }

    pub fn triggers_branch_sync(self) -> bool {
        matches!(self, Self::MemeticBroadcast | Self::Reconciliation)
    }

    pub fn is_player_contact(self) -> bool {
        matches!(
            self,
            Self::FirstContact | Self::MemeticBroadcast | Self::ContactDuel | Self::Reconciliation
        )
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrieRecord {
    pub record_id: u64,
    pub origin_node: usize,
    pub created_year: f64,
    pub kind: TrieEventKind,
    pub delta_assimilate: f64,
    pub delta_defend: f64,
    pub delta_coexist: f64,
    pub signature: u64,
    pub label: String,
}

impl TrieRecord {
    fn compute_digest(
        origin_node: usize,
        created_year: f64,
        kind: TrieEventKind,
        deltas: (f64, f64, f64),
        label: &str,
        nonce: u64,
    ) -> (u64, u64) {
        let mut h: u64 = 0xcbf29ce484222325;
        h = fnv1a_step(h, &(origin_node as u64).to_le_bytes());
        h = fnv1a_step(h, &created_year.to_bits().to_le_bytes());
        h = fnv1a_step(h, &[kind.tag()]);
        h = fnv1a_step(h, &deltas.0.to_bits().to_le_bytes());
        h = fnv1a_step(h, &deltas.1.to_bits().to_le_bytes());
        h = fnv1a_step(h, &deltas.2.to_bits().to_le_bytes());
        h = fnv1a_step(h, label.as_bytes());
        h = fnv1a_step(h, &nonce.to_le_bytes());
        let record_id = splitmix64(h);
        let sig_input = record_id
            ^ splitmix64((origin_node as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93) ^ 0xA5A5_5A5A);
        let signature = splitmix64(sig_input ^ (u64::from(kind.tag()) << 56));
        (record_id, signature)
    }

    pub fn signed(
        origin_node: usize,
        created_year: f64,
        kind: TrieEventKind,
        label: impl Into<String>,
        nonce: u64,
    ) -> Self {
        let mut label_str = label.into();
        if label_str.len() > 120 {
            label_str.truncate(120);
        }
        let deltas = kind.policy_deltas();
        let (record_id, signature) =
            Self::compute_digest(origin_node, created_year, kind, deltas, &label_str, nonce);
        Self {
            record_id,
            origin_node,
            created_year,
            kind,
            delta_assimilate: deltas.0,
            delta_defend: deltas.1,
            delta_coexist: deltas.2,
            signature,
            label: label_str,
        }
    }

    pub fn verify_signature(&self, site_count: usize, max_year: f64) -> bool {
        if self.origin_node >= site_count
            || !self.created_year.is_finite()
            || self.created_year < 0.0
            || self.created_year > max_year
            || self.label.len() > 128
            || !self.delta_assimilate.is_finite()
            || !self.delta_defend.is_finite()
            || !self.delta_coexist.is_finite()
            || self.delta_assimilate.abs() > 2.0
            || self.delta_defend.abs() > 2.0
            || self.delta_coexist.abs() > 2.0
        {
            return false;
        }
        let sig_input = self.record_id
            ^ splitmix64(
                (self.origin_node as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93) ^ 0xA5A5_5A5A,
            );
        let expected_sig = splitmix64(sig_input ^ (u64::from(self.kind.tag()) << 56));
        self.signature == expected_sig
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InFlightTrieMessage {
    pub source: usize,
    pub target: usize,
    pub launch_year: f64,
    pub arrival_year: f64,
    pub record: TrieRecord,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatriciaNode {
    pub node_id: usize,
    pub records: Vec<TrieRecord>,
    pub root_hash: u64,
    pub policy: PolicyWeights,
    pub outgoing: Vec<InFlightTrieMessage>,
}

impl PatriciaNode {
    pub fn new(node_id: usize, genesis: TrieRecord) -> Self {
        let records = vec![genesis];
        let root_hash = compute_patricia_root(&records);
        let policy = PolicyWeights::from_records(&records);
        Self {
            node_id,
            records,
            root_hash,
            policy,
            outgoing: Vec::new(),
        }
    }

    pub fn has_record(&self, record_id: u64) -> bool {
        self.records.iter().any(|r| r.record_id == record_id)
    }

    /// Inserts a record idempotently and recomputes the canonical trie root
    /// and regional policy weights. Returns `true` if the record was new.
    pub fn insert_record(&mut self, record: TrieRecord) -> bool {
        if self.has_record(record.record_id) {
            return false;
        }
        if self.records.len() >= MAX_NODE_RECORDS {
            // Keep genesis + latest records bounded deterministically.
            self.records.remove(1.min(self.records.len() - 1));
        }
        self.records.push(record);
        self.records.sort_by_key(|r| r.record_id);
        self.root_hash = compute_patricia_root(&self.records);
        self.policy = PolicyWeights::from_records(&self.records);
        true
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PatriciaTrieWorld {
    pub nodes: Vec<PatriciaNode>,
    pub next_nonce: u64,
    pub memetic_broadcasts: u32,
    pub contact_records: u32,
}

impl PatriciaTrieWorld {
    pub fn new_for_systems(systems: &[System], year: f64) -> Self {
        let orion_idx = systems
            .iter()
            .position(|s| s.owner == Faction::Blight)
            .unwrap_or(systems.len().saturating_sub(1));
        let remnant_idx = systems
            .iter()
            .position(|s| s.owner == Faction::Drummers)
            .unwrap_or(0);
        let formoi_idx = systems
            .iter()
            .position(|s| s.owner == Faction::Weena)
            .unwrap_or(0);
        let indep_idx = systems
            .iter()
            .position(|s| s.owner == Faction::Independent)
            .unwrap_or(0);

        let genesis_orion = TrieRecord::signed(
            orion_idx,
            0.0_f64.min(year),
            TrieEventKind::GenesisOrion,
            "ORION PRIME BRANCH / ANCIENT ASSIMILATION",
            1,
        );
        let genesis_remnant = TrieRecord::signed(
            remnant_idx,
            0.0_f64.min(year),
            TrieEventKind::GenesisRemnant,
            "BEL-AIR BRANCH / REMNANT COEXISTENCE",
            2,
        );
        let genesis_formoi = TrieRecord::signed(
            formoi_idx,
            0.0_f64.min(year),
            TrieEventKind::GenesisFormoi,
            "NEBULAR WASTES / FORMOI DETERRENCE",
            3,
        );
        let genesis_indep = TrieRecord::signed(
            indep_idx,
            0.0_f64.min(year),
            TrieEventKind::GenesisIndependent,
            "LOCAL ARCHIVE / UNALIGNED CONSENSUS",
            4,
        );

        let nodes = systems
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let genesis = match s.owner {
                    Faction::Blight => genesis_orion.clone(),
                    Faction::Drummers => genesis_remnant.clone(),
                    Faction::Weena => genesis_formoi.clone(),
                    Faction::Independent => genesis_indep.clone(),
                };
                PatriciaNode::new(i, genesis)
            })
            .collect();

        Self {
            nodes,
            next_nonce: 10,
            memetic_broadcasts: 0,
            contact_records: 0,
        }
    }

    /// Number of distinct Merkle-Patricia trie roots currently active across all star systems.
    pub fn branch_divergence_count(&self) -> usize {
        self.nodes
            .iter()
            .map(|n| n.root_hash)
            .collect::<BTreeSet<_>>()
            .len()
    }

    pub fn average_coexist(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.nodes.iter().map(|n| n.policy.coexist).sum();
        sum / self.nodes.len() as f64
    }

    pub fn in_flight_count(&self) -> usize {
        self.nodes.iter().map(|n| n.outgoing.len()).sum()
    }

    pub fn next_arrival_after(&self, year: f64) -> Option<f64> {
        self.nodes
            .iter()
            .flat_map(|n| n.outgoing.iter())
            .filter(|m| m.arrival_year > year)
            .map(|m| m.arrival_year)
            .fold(None, |acc, t| Some(acc.map_or(t, |cur| cur.min(t))))
    }

    fn queue_record_to_peers(
        node: &mut PatriciaNode,
        positions: &[V2],
        year: f64,
        record: &TrieRecord,
    ) -> f64 {
        let mut latest = year;
        let origin_pos = positions[node.node_id];
        for (target, &target_pos) in positions.iter().enumerate() {
            if target == node.node_id {
                continue;
            }
            if node.outgoing.len() >= MAX_OUTGOING_PER_NODE {
                break;
            }
            let distance_ly = (target_pos - origin_pos).len();
            let arrival = year + distance_ly;
            if arrival > TIMELINE_MAX_YR {
                continue;
            }
            let already_queued = node
                .outgoing
                .iter()
                .any(|m| m.target == target && m.record.record_id == record.record_id);
            if already_queued {
                continue;
            }
            latest = latest.max(arrival);
            node.outgoing.push(InFlightTrieMessage {
                source: node.node_id,
                target,
                launch_year: year,
                arrival_year: arrival,
                record: record.clone(),
            });
        }
        latest
    }

    /// Record a signed local event at `source` and emit light-speed (`c = 1 ly/yr`)
    /// propagation packets to every peer system.
    pub fn record_event(
        &mut self,
        source: usize,
        positions: &[V2],
        year: f64,
        kind: TrieEventKind,
        label: impl Into<String>,
    ) -> Result<(u64, f64), &'static str> {
        if source >= self.nodes.len() || positions.len() != self.nodes.len() {
            return Err("INVALID PATRICIA TRIE NODE");
        }
        let nonce = self.next_nonce;
        self.next_nonce += 1;
        let rec = TrieRecord::signed(source, year, kind, label, nonce);
        let id = rec.record_id;
        let node = &mut self.nodes[source];
        node.insert_record(rec.clone());
        let latest = Self::queue_record_to_peers(node, positions, year, &rec);
        if kind.is_player_contact() {
            self.contact_records += 1;
        }
        Ok((id, latest))
    }

    /// Broadcast a player-initiated memetic reconciliation record from a visited
    /// system at `c = 1 ly/yr`. Also propagates the source node's local branch
    /// records so receiving nodes can echo their own diverged branch history
    /// back at light speed, enabling causal Episode XXII convergence.
    pub fn broadcast_memetic_record(
        &mut self,
        source: usize,
        positions: &[V2],
        year: f64,
        label: impl Into<String>,
    ) -> Result<f64, &'static str> {
        if source >= self.nodes.len() || positions.len() != self.nodes.len() {
            return Err("INVALID MEMETIC BROADCAST ORIGIN");
        }
        let nonce = self.next_nonce;
        self.next_nonce += 1;
        let rec = TrieRecord::signed(source, year, TrieEventKind::MemeticBroadcast, label, nonce);
        let node = &mut self.nodes[source];
        node.insert_record(rec);
        let snapshot: Vec<TrieRecord> = node.records.clone();
        let mut latest = year;
        for existing in &snapshot {
            let eta = Self::queue_record_to_peers(node, positions, year, existing);
            latest = latest.max(eta);
        }
        self.memetic_broadcasts += 1;
        self.contact_records += 1;
        Ok(latest)
    }

    /// Advance light-speed propagation up to `year`. Remote nodes only update
    /// their Patricia trie and policy weights upon exact causal arrival.
    pub fn advance_to(&mut self, year: f64, positions: &[V2]) -> Vec<String> {
        let mut log = Vec::new();
        if !year.is_finite() || positions.len() != self.nodes.len() {
            return log;
        }
        // Collect all messages whose physical arrival has occurred.
        let mut arrived: Vec<InFlightTrieMessage> = Vec::new();
        for node in &mut self.nodes {
            let mut still_in_flight = Vec::with_capacity(node.outgoing.len());
            for msg in node.outgoing.drain(..) {
                if msg.arrival_year <= year {
                    arrived.push(msg);
                } else {
                    still_in_flight.push(msg);
                }
            }
            node.outgoing = still_in_flight;
        }
        // Deterministic delivery order by arrival_year, target, record_id.
        arrived.sort_by(|a, b| {
            a.arrival_year
                .total_cmp(&b.arrival_year)
                .then_with(|| a.target.cmp(&b.target))
                .then_with(|| a.record.record_id.cmp(&b.record.record_id))
        });

        for msg in arrived {
            if msg.target >= self.nodes.len() {
                continue;
            }
            let newly_added = self.nodes[msg.target].insert_record(msg.record.clone());
            if newly_added && msg.record.kind.triggers_branch_sync() {
                // Causally echo this node's local branch records back out at c
                // starting at the exact arrival epoch, so diverged branches can
                // converge without ever violating the speed of light.
                let echo_year = msg.arrival_year;
                let to_echo: Vec<TrieRecord> = self.nodes[msg.target]
                    .records
                    .iter()
                    .filter(|r| r.record_id != msg.record.record_id)
                    .cloned()
                    .collect();
                let node = &mut self.nodes[msg.target];
                for rec in &to_echo {
                    Self::queue_record_to_peers(node, positions, echo_year, rec);
                }
                log.push(format!(
                    "TRIE node #{} received memetic broadcast (root {:016x})",
                    msg.target, node.root_hash
                ));
            }
        }
        log
    }

    /// Strict save-file validation: checks node count, canonical trie roots,
    /// record signatures, policy weight bounds, and causal light-cone timing.
    pub fn valid(&self, systems: &[System], year: f64) -> bool {
        if self.nodes.len() != systems.len() {
            return false;
        }
        for (idx, node) in self.nodes.iter().enumerate() {
            if node.node_id != idx
                || node.records.len() > MAX_NODE_RECORDS
                || node.outgoing.len() > MAX_OUTGOING_PER_NODE
                || !node.policy.valid()
            {
                return false;
            }
            for rec in &node.records {
                if !rec.verify_signature(systems.len(), year) {
                    return false;
                }
            }
            if compute_patricia_root(&node.records) != node.root_hash {
                return false;
            }
            let expected_policy = PolicyWeights::from_records(&node.records);
            if (node.policy.assimilate - expected_policy.assimilate).abs() > 1e-6
                || (node.policy.defend - expected_policy.defend).abs() > 1e-6
                || (node.policy.coexist - expected_policy.coexist).abs() > 1e-6
            {
                return false;
            }
            for msg in &node.outgoing {
                if msg.source != idx
                    || msg.target >= systems.len()
                    || msg.target == idx
                    || !msg.launch_year.is_finite()
                    || !msg.arrival_year.is_finite()
                    || msg.launch_year < 0.0
                    || msg.launch_year > year
                    || msg.arrival_year < msg.launch_year
                    || msg.arrival_year > TIMELINE_MAX_YR
                    || !msg.record.verify_signature(systems.len(), year)
                {
                    return false;
                }
                let dist_ly = (systems[msg.target].pos - systems[msg.source].pos).len();
                if msg.arrival_year + 1e-6 < msg.launch_year + dist_ly {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strategy::StrategicState;

    #[test]
    fn deterministic_patricia_root_is_order_independent() {
        let r1 = TrieRecord::signed(0, 100_000.0, TrieEventKind::GenesisOrion, "A", 1);
        let r2 = TrieRecord::signed(1, 100_000.0, TrieEventKind::GenesisRemnant, "B", 2);
        let r3 = TrieRecord::signed(2, 100_050.0, TrieEventKind::ContactDuel, "C", 3);
        let h_abc = compute_patricia_root(&[r1.clone(), r2.clone(), r3.clone()]);
        let h_cba = compute_patricia_root(&[r3.clone(), r2.clone(), r1.clone()]);
        let h_bac = compute_patricia_root(&[r2, r1, r3]);
        assert_eq!(h_abc, h_cba);
        assert_eq!(h_abc, h_bac);
    }

    #[test]
    fn causal_light_cone_delays_remote_trie_updates() {
        let state = StrategicState::new(100_000.0);
        let positions: Vec<V2> = state.systems.iter().map(|s| s.pos).collect();
        let mut world = PatriciaTrieWorld::new_for_systems(&state.systems, state.year);
        let dist_5_6 = (positions[5] - positions[6]).len();
        let root_6_before = world.nodes[6].root_hash;
        let policy_6_before = world.nodes[6].policy;

        world
            .record_event(
                5,
                &positions,
                state.year,
                TrieEventKind::ContactDuel,
                "BEL-AIR CONTACT",
            )
            .unwrap();

        // Just before light arrives at system 6, its trie and policy must be unchanged.
        world.advance_to(state.year + dist_5_6 - 1.0, &positions);
        assert_eq!(world.nodes[6].root_hash, root_6_before);
        assert_eq!(world.nodes[6].policy, policy_6_before);

        // At causal arrival, system 6 updates.
        world.advance_to(state.year + dist_5_6 + 0.01, &positions);
        assert_ne!(world.nodes[6].root_hash, root_6_before);
        assert!(world.nodes[6].policy.coexist > policy_6_before.coexist);
        assert!(world.nodes[6].policy.valid());
    }

    #[test]
    fn distant_systems_start_diverged_and_converge_on_reconciliation() {
        let state = StrategicState::new(100_000.0);
        let positions: Vec<V2> = state.systems.iter().map(|s| s.pos).collect();
        let mut world = PatriciaTrieWorld::new_for_systems(&state.systems, state.year);

        assert!(world.branch_divergence_count() >= 3);
        let start_year = 100_000.0;
        let first_wave_eta = world
            .broadcast_memetic_record(5, &positions, start_year, "AMICUS BRIEF RECONCILIATION")
            .unwrap();
        assert!(first_wave_eta > start_year);

        // Mid-flight: still diverged.
        world.advance_to(start_year + 500.0, &positions);
        assert!(world.branch_divergence_count() > 1);

        // Step forward causally in chunks until all light-speed echoes arrive.
        let mut t = start_year + 500.0;
        while let Some(next_t) = world.next_arrival_after(t) {
            t = next_t;
            world.advance_to(t, &positions);
        }
        assert_eq!(world.branch_divergence_count(), 1);
        assert!(world.average_coexist() >= 0.45);
        assert!(world.valid(&state.systems, t));
    }

    #[test]
    fn policy_weights_stay_bounded_and_normalized() {
        let state = StrategicState::new(100_000.0);
        let positions: Vec<V2> = state.systems.iter().map(|s| s.pos).collect();
        let mut world = PatriciaTrieWorld::new_for_systems(&state.systems, state.year);
        for k in 0..10 {
            world
                .record_event(
                    9,
                    &positions,
                    100_000.0 + k as f64,
                    TrieEventKind::AssimilateDirective,
                    format!("ASSIMILATE #{k}"),
                )
                .unwrap();
        }
        let p = world.nodes[9].policy;
        assert!(p.valid());
        assert!(p.assimilate > 0.8);
        assert!(p.coexist > 0.0 && p.defend > 0.0);
    }
}
