//! Deterministic 3D cell-list neighbor search for compact-support SPH kernels.
 //!
 //! Cell size is chosen as the *maximum support radius* (2*h_max).
 //! Any pair that interacts through W_i or W_j lies in the same cell or
 //! one of 26 adjacent cells. Query results are sorted into particle-index
 //! order so reductions remain reproducible between runs and tests.
 //!
 //! This accelerates density/grad-h and pressure neighbor discovery. The
 //! Newtonian long-range gravity remains all-pairs until a validated
 //! multipole tree or FMM backend is available.
use crate::sph3::{Fluid3,V3};
use std::collections::HashMap;

#[derive(Debug)]
pub struct CellNeighbors {
    edge_m:f64,
    cells:HashMap<(i64,i64,i64),Vec<usize>>,
    cell_of:Vec<(i64,i64,i64)>,
}
impl CellNeighbors {
    pub fn new(particles:&[Fluid3],maximum_h_m:f64)->Self {
        assert!(maximum_h_m.is_finite()&&maximum_h_m>0.0);
        let edge=2.0*maximum_h_m;
        let mut cells=HashMap::<(i64,i64,i64),Vec<usize>>::new();
        let mut cell_of=Vec::with_capacity(particles.len());
        for (i,p) in particles.iter().enumerate() {
            assert!(p.pos_m.finite());
            let cell=Self::key(p.pos_m,edge);
            cells.entry(cell).or_default().push(i);
            cell_of.push(cell);
        }
        Self{edge_m:edge,cells,cell_of}
    }
    fn key(p:V3,edge:f64)->(i64,i64,i64) {
        ((p.x/edge).floor() as i64,
         (p.y/edge).floor() as i64,
         (p.z/edge).floor() as i64)
    }
    pub fn candidates(&self,index:usize)->Vec<usize> {
        let (x,y,z)=self.cell_of[index];
        let mut list=Vec::new();
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if let Some(candidates)=self.cells.get(&(
                        x.saturating_add(dx),y.saturating_add(dy),
                        z.saturating_add(dz),
                    )) {
                        list.extend_from_slice(candidates);
                    }
                }
            }
        }
        list.sort_unstable();
        list
    }
    /// Prepare neighbor candidates once per force evaluation. Positions
    /// remain fixed during the density-h iterative closure.
    pub fn all_candidates(&self)->Vec<Vec<usize>> {
        (0..self.cell_of.len()).map(|i|self.candidates(i)).collect()
    }
    pub fn mean_candidate_count(&self)->f64 {
        if self.cell_of.is_empty(){return 0.0;}
        (0..self.cell_of.len()).map(|i|self.candidates(i).len())
            .sum::<usize>() as f64/self.cell_of.len() as f64
    }
    pub fn cell_edge_m(&self)->f64 {self.edge_m}
    pub fn occupied_cells(&self)->usize {self.cells.len()}
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    use crate::sph3;
    use crate::white_dwarf::{DegenerateElectrons,profile_for_mass};
    fn fixture(n:usize)->Vec<Fluid3> {
        let eos=DegenerateElectrons::new(2.0).unwrap();
        let star=profile_for_mass(0.65*SOLAR_MASS_KG,eos,10_000.0).unwrap();
        sph3::seed_from_hydrostatic(&star,V3::new(1e8,0.0,-2e8),
            V3::default(),n)
    }
    #[test]
    fn all_kernel_pairs_are_covered_including_cell_boundaries() {
        for n in [32,128,512] {
            let particles=fixture(n);
            let max_h=2.4e6;
            let grid=CellNeighbors::new(&particles,max_h);
            let all=grid.all_candidates();
            for i in 0..n {
                assert!(all[i].contains(&i));
                assert!(all[i].windows(2).all(|x|x[0]<x[1]));
                for j in 0..n {
                    let distance=(particles[i].pos_m-particles[j].pos_m).len();
                    if distance<2.0*max_h {
                        assert!(all[i].contains(&j),
                            "missed pair {i} {j} at r {distance}");
                    }
                }
            }
        }
    }
    #[test]
    fn negative_and_large_coordinates_are_indexed_consistently() {
        let points=[
            V3::new(-1e12,-1e12,-1e12),
            V3::new(-1e12+100.0,-1e12+90.0,-1e12+80.0),
            V3::new(1e12,1e12,1e12),
        ].map(|p|Fluid3{pos_m:p,vel_m_s:V3::default(),mass_kg:1.0});
        let cells=CellNeighbors::new(&points,1e6);
        assert_eq!(cells.candidates(0),vec![0,1]);
        assert_eq!(cells.candidates(2),vec![2]);
    }
    #[test]
    fn deterministic_candidate_order_across_hash_table_instances() {
        let p=fixture(192);
        let a=CellNeighbors::new(&p,3e6);
        let b=CellNeighbors::new(&p,3e6);
        assert_eq!(a.all_candidates(),b.all_candidates());
        assert!(a.occupied_cells()>1);
        assert!(a.mean_candidate_count()<p.len() as f64);
    }
}
