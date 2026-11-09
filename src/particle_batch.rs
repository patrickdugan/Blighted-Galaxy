//! GPU-backed batched 2D particle billboards for macroquad.
 //!
 //! Particle positions are computed on the CPU; this builder emits indexed
 //! triangles for GPU rasterization with bounded draw calls. This is not a
 //! GPU compute-particle implementation or an instanced renderer.
use macroquad::models::{draw_mesh, Mesh, Vertex};
use macroquad::prelude::{Color,Vec2};

/// 4 vertices/particle. Restrict per-batch indices to u16 range.
pub const MAX_PARTICLES_PER_MESH:usize=4096;
#[derive(Debug,Clone,Copy)]
pub struct ParticleBillboard {
    pub center_px:Vec2,
    pub radius_px:f32,
    pub tint:Color,
    pub depth:f32,
}
impl ParticleBillboard {
    fn valid(&self)->bool {
        self.center_px.x.is_finite()&&self.center_px.y.is_finite()&&
        self.radius_px.is_finite()&&self.radius_px>0.0&&
        self.depth.is_finite()
    }
}
pub fn build_meshes(particles:&[ParticleBillboard])->Vec<Mesh> {
    let mut result=Vec::with_capacity(
        particles.len().div_ceil(MAX_PARTICLES_PER_MESH));
    for batch in particles.chunks(MAX_PARTICLES_PER_MESH) {
        let mut vertices=Vec::with_capacity(batch.len()*4);
        let mut indices=Vec::with_capacity(batch.len()*6);
        for p in batch {
            if !p.valid() {continue;}
            let cx=p.center_px.x;let cy=p.center_px.y;
            let r=p.radius_px;
            let b=vertices.len() as u16;
            vertices.push(Vertex::new(cx-r,cy-r,p.depth,0.0,0.0,p.tint));
            vertices.push(Vertex::new(cx+r,cy-r,p.depth,1.0,0.0,p.tint));
            vertices.push(Vertex::new(cx+r,cy+r,p.depth,1.0,1.0,p.tint));
            vertices.push(Vertex::new(cx-r,cy+r,p.depth,0.0,1.0,p.tint));
            indices.extend_from_slice(&[b,b+1,b+2,b,b+2,b+3]);
        }
        if !vertices.is_empty() {
            result.push(Mesh{vertices,indices,texture:None});
        }
    }
    result
}
/// Submit indexed particle triangles, one GPU draw per mesh segment.
pub fn draw_billboards(meshes:&[Mesh]) {
    for mesh in meshes {draw_mesh(mesh);}
}
pub fn count_indexed_quads(meshes:&[Mesh])->usize {
    meshes.iter().map(|m|m.indices.len()/6).sum()
}
#[cfg(test)]
mod tests {
    use super::*;
    fn sprite(i:usize)->ParticleBillboard {
        ParticleBillboard {
            center_px:Vec2::new(i as f32*0.1,300.0),
            radius_px:2.0,
            tint:Color::new(0.2,0.8,1.0,1.0),
            depth:0.0,
        }
    }
    #[test]
    fn triangles_per_billboard_and_valid_indices() {
        let meshes=build_meshes(&(0..100).map(sprite).collect::<Vec<_>>());
        assert_eq!(meshes.len(),1);
        assert_eq!(meshes[0].vertices.len(),400);
        assert_eq!(meshes[0].indices.len(),600);
        assert!(meshes[0].indices.iter().all(|i|
            (*i as usize)<meshes[0].vertices.len()));
    }
    #[test]
    fn segmented_batches_stay_under_u16_indices() {
        let meshes=build_meshes(&(0..10_001).map(sprite).collect::<Vec<_>>());
        assert_eq!(meshes.len(),3);
        assert_eq!(count_indexed_quads(&meshes),10_001);
        assert!(meshes.iter().all(|m|m.vertices.len()<=u16::MAX as usize));
        for m in &meshes {
            assert!(m.indices.iter().all(|i|(*i as usize)<m.vertices.len()));
        }
    }
    #[test]
    fn invalid_billboards_are_ignored() {
        let mut p=sprite(0);
        p.radius_px=f32::NAN;
        assert_eq!(count_indexed_quads(&build_meshes(&[p])),0);
    }
    #[test]
    fn empty_meshes_produce_no_draw_calls() {
        assert!(build_meshes(&[]).is_empty());
    }
}
