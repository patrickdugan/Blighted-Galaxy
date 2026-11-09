//! Three-dimensional Kerr photon camera using separable, non-equatorial
//! null geodesics in the stationary Kerr metric.
 //!
 //! Rays launch toward the scene from a local ZAMO tetrad, with past-directed
 //! four-momentum. They integrate (u=1/r, theta, phi, t) in Mino time, normalized
 //! to |E_infinity|=1, with nonzero Carter separation constant Q when needed.
 //! Boyer-Lindquist horizon boundary, opaque equatorial thin disc, and a
 //! synthetic star background are explicit visual choices. This is not a
 //! radiation-transfer or plasma-accretion simulation.
use crate::kerr::KerrGeometry;
use std::f64::consts::{PI,TAU};

const MAX_H:f64=0.004;
const HORIZON_MARGIN:f64=0.04;
const MAX_STEPS:usize=1600;
const TINY:f64=1e-12;

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
pub enum RayExit { Horizon, Sky, Disc, Unresolved }
#[derive(Debug,Clone,Copy)]
pub struct KerrPixelRay {
    pub u:f64,
    pub du:f64,
    pub theta:f64,
    pub dtheta:f64,
    pub phi:f64,
    pub t_geom:f64,
    pub b:f64,
    pub q:f64,
    pub steps:usize,
    pub exit:RayExit,
    pub disc_radius_rg:f64,
    pub disc_g_factor:f64,
    pub max_radial_error:f64,
    pub max_polar_error:f64,
}
impl KerrPixelRay {
    pub fn radius_rg(&self)->f64 {1.0/self.u}
    fn radial_polynomial(&self,a:f64)->f64 {
        let u=self.u;let k=(self.b-a).powi(2)+self.q;
        1.0+(a*a-self.b*self.b-self.q)*u*u+
            2.0*k*u.powi(3)-a*a*self.q*u.powi(4)
    }
    fn polar_polynomial(&self,a:f64)->f64 {
        let sin=self.theta.sin();
        if sin.abs()<1e-10 {return f64::NAN;}
        self.q+a*a*self.theta.cos().powi(2)-
            self.b*self.b*(self.theta.cos()/sin).powi(2)
    }
    fn state(&self)->[f64;6] {
        [self.u,self.du,self.theta,self.dtheta,self.phi,self.t_geom]
    }
    fn set_state(&mut self,s:[f64;6]) {
        self.u=s[0];self.du=s[1];self.theta=s[2];
        self.dtheta=s[3];self.phi=s[4];self.t_geom=s[5];
        if self.theta<0.0 {
            self.theta=-self.theta;self.dtheta=-self.dtheta;self.phi+=PI;
        } else if self.theta>PI {
            self.theta=2.0*PI-self.theta;
            self.dtheta=-self.dtheta;self.phi+=PI;
        }
    }
    pub fn invariant_errors(&self,a:f64)->(f64,f64) {
        (self.du*self.du-self.radial_polynomial(a),
         self.dtheta*self.dtheta-self.polar_polynomial(a))
    }
}
#[derive(Debug,Clone,Copy)]
pub struct KerrCamera {
    pub geometry:KerrGeometry,
    /// ZAMO observer's Boyer-Lindquist radius in gravitational radii.
    pub observer_rg:f64,
    pub inclination_rad:f64,
    /// Horizontal field of view expressed as tangent-plane length in rg.
    pub span_rg:f64,
    /// Thin disc cutoff radii. Negative spin chooses negative-L circular gas.
    pub disc_outer_rg:f64,
}
impl KerrCamera {
    pub fn new(geometry:KerrGeometry,observer_rg:f64,
               inclination_rad:f64,span_rg:f64)->Option<Self> {
        if !observer_rg.is_finite()||
            observer_rg<=geometry.outer_horizon_rg()+5.0||
            !inclination_rad.is_finite()||
            !(0.06..=PI-0.06).contains(&inclination_rad)||
            !span_rg.is_finite()||span_rg<=0.0 {return None;}
        Some(Self{geometry,observer_rg,inclination_rad,span_rg,
            disc_outer_rg:28.0})
    }
    pub fn disc_inner_rg(&self)->f64 {
        let direction=if self.geometry.spin<0.0 {-1} else {1};
        self.geometry.isco_rg(direction).unwrap()
    }
    /// Exact finite-radius local orthonormal (ZAMO) observer launch.
    /// A past-directed null ray is shot backward into the scene. The
    /// alpha/beta screen coordinates are pinhole tangent-plane lengths (rg).
    pub fn launch(&self,alpha_rg:f64,beta_rg:f64)->Option<KerrPixelRay> {
        if !alpha_rg.is_finite()||!beta_rg.is_finite() {return None;}
        let r=self.observer_rg;
        let a=self.geometry.spin;
        let theta=self.inclination_rad;
        let sin=theta.sin();let cos=theta.cos();
        let sigma=r*r+a*a*cos*cos;
        let delta=self.geometry.delta(r);
        let big_a=(r*r+a*a).powi(2)-a*a*delta*sin*sin;
        let lapse=(sigma*delta/big_a).sqrt();
        let omega=2.0*a*r/big_a;
        let root_g_phi=(big_a*sin*sin/sigma).sqrt();
        // Local photon direction in the ZAMO tetrad. Time component is
        // negative: ray points into the PAST of the observer.
        let norm=(r*r+alpha_rg*alpha_rg+beta_rg*beta_rg).sqrt();
        let nphi=alpha_rg/norm;
        let ntheta=-beta_rg/norm;
        let angular_momentum=nphi*root_g_phi;
        let e_signed=-lapse+omega*angular_momentum;
        if e_signed>=-1e-12 {return None;}
        let b=angular_momentum/e_signed;
        let q=sigma*ntheta*ntheta/(e_signed*e_signed)+
            cos*cos*(b*b/(sin*sin)-a*a);
        let u=1.0/r;
        let mut ray=KerrPixelRay {
            u,du:0.0,theta,dtheta:0.0,phi:0.0,t_geom:0.0,
            b,q,steps:0,exit:RayExit::Unresolved,
            disc_radius_rg:0.0,disc_g_factor:0.0,
            max_radial_error:0.0,max_polar_error:0.0,
        };
        ray.du=ray.radial_polynomial(a).max(0.0).sqrt();
        let pol=ray.polar_polynomial(a).max(0.0).sqrt();
        ray.dtheta=if beta_rg>=0.0 {-pol}else{pol};
        Some(ray)
    }
    fn deriv(&self,ray:&KerrPixelRay,state:[f64;6])->[f64;6] {
        let u=state[0];
        let theta=state[2];
        let a=self.geometry.spin;
        let b=ray.b;let q=ray.q;
        let k=(b-a).powi(2)+q;
        let u2=u*u;
        let rad_acc=(a*a-b*b-q)*u+
            3.0*k*u2-2.0*a*a*q*u2*u;
        let sin=theta.sin();
        let safe_sin=if sin.abs()<1e-8 {1e-8}else{sin};
        let cos=theta.cos();
        let polar_acc=-a*a*sin*cos+b*b*cos/safe_sin.powi(3);
        let r=1.0/u.max(1e-10);
        let delta=self.geometry.delta(r);
        let p=r*r+a*a-a*b;
        // Negative signs are essential for a past-directed ray.
        let phidot=-(a*p/delta+b/(safe_sin*safe_sin)-a);
        let tdot=-((r*r+a*a)*p/delta+a*(b-a*sin*sin));
        [state[1],rad_acc,state[3],polar_acc,phidot,tdot]
    }
    fn rk4(&self,ray:&KerrPixelRay,y:[f64;6],h:f64)->[f64;6] {
        let k1=self.deriv(ray,y);
        let y2=std::array::from_fn(|i|y[i]+0.5*h*k1[i]);
        let k2=self.deriv(ray,y2);
        let y3=std::array::from_fn(|i|y[i]+0.5*h*k2[i]);
        let k3=self.deriv(ray,y3);
        let y4=std::array::from_fn(|i|y[i]+h*k3[i]);
        let k4=self.deriv(ray,y4);
        std::array::from_fn(|i|y[i]+
            h*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0)
    }
    fn disc_redshift(&self,r:f64,b:f64)->f64 {
        let s=if self.geometry.spin<0.0 {-1}else{1};
        let orbit=match self.geometry.circular_timelike(r,s) {
            Some(v)=>v,None=>return 1.0,
        };
        let a=self.geometry.spin;
        let omega=orbit.omega_rad_s*self.geometry.r_g_m/
            crate::physics::C_SI;
        let gtt=-(1.0-2.0/r);
        let gtphi=-2.0*a/r;
        let gphiphi=r*r+a*a+2.0*a*a/r;
        let norm=-(gtt+2.0*gtphi*omega+gphiphi*omega*omega);
        if norm<=0.0 {return 1.0;}
        let ut=1.0/norm.sqrt();
        let denom=ut*(1.0-omega*b);
        if denom>0.0 {1.0/denom}else{1.0}
    }
    pub fn trace(&self,alpha_rg:f64,beta_rg:f64)->Option<KerrPixelRay> {
        let mut ray=self.launch(alpha_rg,beta_rg)?;
        let u_obs=1.0/self.observer_rg;
        let u_cut=1.0/(self.geometry.outer_horizon_rg()+HORIZON_MARGIN);
        for step in 0..MAX_STEPS {
            let mut h=MAX_H;
            if ray.du>0.0 {
                h=h.min(0.35*(u_cut-ray.u)/ray.du.max(1e-9));
            }
            if !h.is_finite()||h<1e-10 {
                ray.exit=RayExit::Horizon;break;
            }
            let old=ray.state();
            let next=self.rk4(&ray,old,h);
            if !next.iter().all(|v|v.is_finite()) {
                ray.exit=RayExit::Horizon;break;
            }
            if next[0]>=u_cut {
                ray.exit=RayExit::Horizon;break;
            }
            ray.set_state(next);
            ray.steps=step+1;
            let (er,ep)=ray.invariant_errors(self.geometry.spin);
            if ray.radius_rg()>self.geometry.outer_horizon_rg()+0.6 {
                if er.is_finite() {ray.max_radial_error=ray.max_radial_error.max(er.abs());}
                if ep.is_finite() {ray.max_polar_error=ray.max_polar_error.max(ep.abs());}
            }
            let previous=old[2]-PI*0.5;
            let following=ray.theta-PI*0.5;
            if (previous*following<0.0 || following.abs()<1e-9) &&
                (old[2]-ray.theta).abs()>1e-10 {
                let f=((PI*0.5-old[2])/(ray.theta-old[2])).clamp(0.0,1.0);
                let u_cross=old[0]+f*(ray.u-old[0]);
                let r_cross=1.0/u_cross;
                if r_cross>=self.disc_inner_rg() &&
                    r_cross<=self.disc_outer_rg {
                    ray.exit=RayExit::Disc;
                    ray.disc_radius_rg=r_cross;
                    ray.disc_g_factor=self.disc_redshift(r_cross,ray.b);
                    break;
                }
            }
            if ray.u<=u_obs && ray.du<0.0 {
                ray.exit=RayExit::Sky;break;
            }
        }
        Some(ray)
    }
    fn color(&self,ray:&KerrPixelRay)->[u8;4] {
        match ray.exit {
            RayExit::Horizon=>[3,5,12,255],
            RayExit::Unresolved=>[28,8,48,255],
            RayExit::Disc=>{
                let r=ray.disc_radius_rg;
                let scale=(self.disc_inner_rg()/r).powf(0.75);
                let heat=(scale*ray.disc_g_factor).clamp(0.0,2.2);
                let brightness=(0.25+0.75*(heat/1.6).min(1.0)).clamp(0.0,1.0);
                let red=(235.0*brightness).min(255.0);
                let green=(85.0+165.0*heat/1.8).min(255.0)*brightness;
                let blue=(45.0+185.0*heat/1.5).min(255.0)*brightness;
                [red as u8,green as u8,blue as u8,255]
            }
            RayExit::Sky=>{
                // Deterministic background field in the terminal angular
                // coordinates. The strong bending distorts this source map.
                let sx=ray.theta.sin()*ray.phi.cos();
                let sy=ray.theta.sin()*ray.phi.sin();
                let sz=ray.theta.cos();
                let band=(sz*5.0+0.35*ray.phi.sin()).abs();
                let haze=(1.0/(1.0+band*band*8.0))*24.0;
                let stars=((sx*181.37+sy*293.87+sz*117.21).sin()*
                    (sx*421.91-sy*97.13+sz*379.87).cos()).abs();
                let star=if stars>0.994 {225.0}else{0.0};
                [(7.0+haze*0.65+star).min(255.0) as u8,
                 (10.0+haze+star).min(255.0) as u8,
                 (24.0+haze*1.3+star).min(255.0) as u8,255]
            }
        }
    }
    pub fn render_rgba(&self,width:u16,height:u16)->CameraImage {
        assert!(width>0&&height>0);
        let mut pixels=Vec::with_capacity(width as usize*height as usize*4);
        let mut image=CameraImage{
            width,height,pixels:vec![],captured:0,sky:0,disc:0,unresolved:0,
            worst_radial_residual:0.0,worst_polar_residual:0.0,
            traced_steps:0,
        };
        let xscale=self.span_rg/width as f64;
        for y in 0..height {
            for x in 0..width {
                let alpha=(x as f64+0.5-width as f64/2.0)*xscale;
                let beta=(height as f64/2.0-y as f64-0.5)*xscale;
                let ray=self.trace(alpha,beta).expect("valid camera pixel");
                pixels.extend_from_slice(&self.color(&ray));
                match ray.exit {
                    RayExit::Horizon=>image.captured+=1,
                    RayExit::Sky=>image.sky+=1,
                    RayExit::Disc=>image.disc+=1,
                    RayExit::Unresolved=>image.unresolved+=1,
                }
                image.traced_steps+=ray.steps as u64;
                image.worst_radial_residual=image.worst_radial_residual.max(ray.max_radial_error);
                image.worst_polar_residual=image.worst_polar_residual.max(ray.max_polar_error);
            }
        }
        image.pixels=pixels;image
    }
}
#[derive(Debug)]
pub struct CameraImage {
    pub width:u16,
    pub height:u16,
    pub pixels:Vec<u8>,
    pub captured:usize,
    pub sky:usize,
    pub disc:usize,
    pub unresolved:usize,
    pub worst_radial_residual:f64,
    pub worst_polar_residual:f64,
    pub traced_steps:u64,
}
impl CameraImage {
    pub fn mean_steps_per_ray(&self)->f64 {
        self.traced_steps as f64/(self.width as f64*self.height as f64)
    }
    pub fn ppm_bytes(&self)->Vec<u8> {
        let mut out=format!("P6\n{} {}\n255\n",self.width,self.height).into_bytes();
        for p in self.pixels.chunks_exact(4) {out.extend_from_slice(&p[0..3]);}
        out
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::SOLAR_MASS_KG;
    fn camera(a:f64,i:f64)->KerrCamera {
        KerrCamera::new(
            KerrGeometry::new(10.0*SOLAR_MASS_KG,a).unwrap(),
            85.0,i,36.0).unwrap()
    }
    #[test]
    fn finite_radius_zamo_tetrad_satisfies_both_null_invariants() {
        let cam=camera(0.85,PI/3.0);
        for (x,y) in [(0.0,0.0),(5.0,7.0),(-10.0,4.0)] {
            let r=cam.launch(x,y).unwrap();
            let (rad,pol)=r.invariant_errors(cam.geometry.spin);
            assert!(rad.abs()<1e-10,"radial {rad}");
            assert!(pol.abs()<1e-10,"polar {pol}");
        }
    }
    #[test]
    fn central_ray_reaches_black_hole_silhouette() {
        let cam=camera(0.0,PI/3.0);
        let r=cam.trace(0.0,0.0).unwrap();
        assert_eq!(r.exit,RayExit::Horizon);
        assert!(r.t_geom<0.0);
    }
    #[test]
    fn distant_ray_escapes_with_correct_time_orientation() {
        let cam=camera(0.0,PI/3.0);
        let r=cam.trace(15.0,10.0).unwrap();
        assert_eq!(r.exit,RayExit::Sky);
        assert!(r.t_geom<0.0);
    }
    #[test]
    fn off_equatorial_pixel_traces_nonzero_carter_constant() {
        let cam=camera(0.8,PI/3.0);
        let r=cam.launch(3.0,8.0).unwrap();
        assert!(r.q.abs()>1.0);
        assert!(r.dtheta.abs()>0.01);
        let trace=cam.trace(3.0,8.0).unwrap();
        assert!((trace.theta-cam.inclination_rad).abs()>0.005);
    }
    #[test]
    fn equatorial_ray_remains_equatorial_in_schwarzschild() {
        let cam=camera(0.0,PI/2.0);
        let trace=cam.trace(7.0,0.0).unwrap();
        assert!((trace.theta-PI/2.0).abs()<1e-9);
    }
    #[test]
    fn render_camera_has_silhouette_disc_and_sky() {
        let cam=camera(0.85,PI/3.0);
        let image=cam.render_rgba(48,32);
        assert_eq!(image.pixels.len(),48*32*4);
        assert_eq!(image.captured+image.disc+image.sky+image.unresolved,48*32);
        assert!(image.captured>0);
        assert!(image.sky>0);
        assert!(image.disc>0);
        assert!(image.mean_steps_per_ray()<MAX_STEPS as f64);
        assert_eq!(image.ppm_bytes()[0..2],*b"P6");
    }
    #[test]
    fn spin_changes_silhouette_or_disc_ray_distribution() {
        let left=camera(0.0,PI/3.0).render_rgba(40,28);
        let right=camera(0.90,PI/3.0).render_rgba(40,28);
        assert_ne!(left.pixels,right.pixels);
    }
    #[test]
    fn disc_shift_is_positive_and_finite() {
        let cam=camera(0.85,PI/3.0);
        let r=cam.disc_inner_rg()+4.0;
        let g=cam.disc_redshift(r,2.0);
        assert!(g>0.0&&g.is_finite());
    }
}
