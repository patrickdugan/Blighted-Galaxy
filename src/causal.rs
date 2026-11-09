//! Bounded causal worldline histories for delayed sensor and projectile events.
//! Positions have SI units and shared inertial coordinate time. Observation is
//! solved strictly inside recorded history; extrapolating nonexistent past
//! events is prohibited. Samples are connected by piecewise-linear paths.
use crate::physics::C_SI;
use crate::V2;
#[derive(Debug,Clone,Copy)]
pub struct WorldlineSample {
    pub time_s:f64,
    pub pos_m:V2,
    pub vel_m_s:V2,
}
#[derive(Debug,Clone,Copy)]
pub struct CausalObservation {
    pub observed_at_s:f64,
    pub emitted_at_s:f64,
    pub delay_s:f64,
    pub pos_m:V2,
    /// Piecewise-interpolated historical estimate, not an instantaneous velocity.
    pub vel_m_s:V2,
}
#[derive(Debug,Clone)]
pub struct WorldlineHistory {
    samples:Vec<WorldlineSample>,
    limit:usize,
}
impl WorldlineHistory {
    pub fn new(max_samples:usize)->Self {
        assert!(max_samples>=2);
        Self{samples:Vec::with_capacity(max_samples.min(8192)),limit:max_samples}
    }
    pub fn len(&self)->usize {self.samples.len()}
    pub fn is_empty(&self)->bool {self.samples.is_empty()}
    pub fn oldest_time(&self)->Option<f64> {self.samples.first().map(|s|s.time_s)}
    pub fn latest_time(&self)->Option<f64> {self.samples.last().map(|s|s.time_s)}
    pub fn push(&mut self,sample:WorldlineSample)->Result<(),&'static str> {
        if !sample.time_s.is_finite()||
            !sample.pos_m.x.is_finite()||!sample.pos_m.y.is_finite()||
            !sample.vel_m_s.x.is_finite()||!sample.vel_m_s.y.is_finite()||
            sample.vel_m_s.len()>=C_SI {return Err("non-finite or noncausal sample");}
        if self.samples.last().is_some_and(|s|sample.time_s<=s.time_s) {
            return Err("worldline time must increase strictly");
        }
        self.samples.push(sample);
        if self.samples.len()>self.limit {self.samples.remove(0);}
        Ok(())
    }
    /// Interpolation is defined only between known sample times, inclusive.
    pub fn sample_at(&self,time_s:f64)->Option<WorldlineSample> {
        if !time_s.is_finite()||self.samples.is_empty(){return None;}
        let first=self.samples.first()?;
        let last=self.samples.last()?;
        if time_s<first.time_s||time_s>last.time_s {return None;}
        let idx=self.samples.partition_point(|s|s.time_s<time_s);
        if idx==0 {return Some(*first);}
        if idx==self.samples.len() {return Some(*last);}
        let right=self.samples[idx];
        if right.time_s==time_s {return Some(right);}
        let left=self.samples[idx-1];
        let fraction=(time_s-left.time_s)/(right.time_s-left.time_s);
        let pos=left.pos_m+(right.pos_m-left.pos_m)*fraction;
        let segment_velocity=(right.pos_m-left.pos_m)/(right.time_s-left.time_s);
        // The segment derivative ensures the interpolated worldline remains
        // causal if each recorded segment is causal.
        Some(WorldlineSample{time_s,pos_m:pos,vel_m_s:segment_velocity})
    }
    /// Solve |source(t_emit) - observer(t_obs)| = c*(t_obs-t_emit).
    /// For subluminal piecewise-linear motion, the causal intersection, if
    /// bracketted in the history, is unique. Return None when the source's
    /// corresponding emission event predates the stored history.
    pub fn observe(&self,observer_m:V2,observer_time_s:f64)->Option<CausalObservation> {
        if !observer_time_s.is_finite()||
            !observer_m.x.is_finite()||!observer_m.y.is_finite()||
            self.samples.len()<2 {return None;}
        let first=self.samples.first()?.time_s;
        let last=self.samples.last()?.time_s.min(observer_time_s);
        if last<first{return None;}
        let earliest=self.sample_at(first)?;
        let latest=self.sample_at(last)?;
        let light_margin=|s:WorldlineSample|->
            f64 {C_SI*(observer_time_s-s.time_s)-(s.pos_m-observer_m).len()};
        let early=light_margin(earliest);
        let late=light_margin(latest);
        // No history sufficient to determine the actual emission event.
        if early<0.0||late>0.0 {return None;}
        let mut lo=first;let mut hi=last;
        // 75 bisections yield sub-ns resolution for ordinary encounter spans.
        for _ in 0..75 {
            let mid=lo+(hi-lo)*0.5;
            if mid==lo||mid==hi {break;}
            let event=self.sample_at(mid)?;
            if light_margin(event)>=0.0 {lo=mid;} else {hi=mid;}
        }
        let emitted_at=lo+(hi-lo)*0.5;
        let event=self.sample_at(emitted_at)?;
        Some(CausalObservation{
            observed_at_s:observer_time_s,
            emitted_at_s:emitted_at,
            delay_s:(observer_time_s-emitted_at).max(0.0),
            pos_m:event.pos_m,
            vel_m_s:event.vel_m_s,
        })
    }
    pub fn validate_segments(&self)->bool {
        self.samples.windows(2).all(|pair|{
            let a=pair[0];let b=pair[1];
            (b.pos_m-a.pos_m).len()<(b.time_s-a.time_s)*C_SI
        })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn check(a:f64,b:f64,absolute:f64) {
        assert!((a-b).abs()<=absolute,
            "actual={a:.16e} expected={b:.16e}");
    }
    fn sample(t:f64,x:f64)->WorldlineSample {
        WorldlineSample{time_s:t,pos_m:V2::new(x,0.0),vel_m_s:V2::default()}
    }
    #[test]
    fn stationary_source_at_one_light_second() {
        let mut history=WorldlineHistory::new(50);
        for t in 0..=20 {history.push(sample(t as f64,C_SI)).unwrap();}
        let obs=history.observe(V2::default(),9.0).unwrap();
        check(obs.emitted_at_s,8.0,1e-12);
        check(obs.delay_s,1.0,1e-12);
        check(obs.pos_m.x,C_SI,1e-4);
    }
    #[test]
    fn moving_source_requires_retarded_position() {
        let mut history=WorldlineHistory::new(64);
        for t in 0..=20 {
            let t=t as f64/2.0;
            history.push(WorldlineSample{time_s:t,
                pos_m:V2::new(C_SI*(1.0+0.5*t),0.0),
                vel_m_s:V2::new(0.5*C_SI,0.0)}).unwrap();
        }
        let obs=history.observe(V2::default(),3.0).unwrap();
        check(obs.emitted_at_s,4.0/3.0,1e-12);
        check(obs.delay_s,5.0/3.0,1e-12);
        check(obs.pos_m.x,(5.0/3.0)*C_SI,1e-3);
    }
    #[test]
    fn accelerated_history_interpolates_measured_trajectory() {
        let mut history=WorldlineHistory::new(20);
        history.push(sample(0.0,C_SI)).unwrap();
        history.push(sample(1.0,1.5*C_SI)).unwrap();
        history.push(sample(2.0,1.5*C_SI)).unwrap();
        history.push(sample(3.0,2.0*C_SI)).unwrap();
        let mid=history.sample_at(1.5).unwrap();
        check(mid.pos_m.x,1.5*C_SI,1e-5);
        assert!(history.validate_segments());
        let obs=history.observe(V2::default(),3.0).unwrap();
        check(obs.emitted_at_s,1.5,1e-12);
    }
    #[test]
    fn missing_earlier_emission_returns_none() {
        let mut history=WorldlineHistory::new(10);
        history.push(sample(0.0,10.0*C_SI)).unwrap();
        history.push(sample(1.0,10.0*C_SI)).unwrap();
        assert!(history.observe(V2::default(),2.0).is_none());
    }
    #[test]
    fn history_is_strictly_bounded() {
        let mut history=WorldlineHistory::new(3);
        for i in 0..12 {history.push(sample(i as f64,0.0)).unwrap();}
        assert_eq!(history.len(),3);
        assert_eq!(history.oldest_time(),Some(9.0));
        assert!(history.sample_at(3.0).is_none());
        assert!(history.observe(V2::default(),12.0).is_none());
    }
    #[test]
    fn invalid_order_and_superluminal_velocity_rejected() {
        let mut history=WorldlineHistory::new(3);
        history.push(sample(1.0,0.0)).unwrap();
        assert!(history.push(sample(1.0,0.0)).is_err());
        assert!(history.push(sample(0.5,0.0)).is_err());
        assert!(history.push(WorldlineSample{time_s:2.0,pos_m:V2::default(),
            vel_m_s:V2::new(C_SI,0.0)}).is_err());
    }
    #[test]
    fn same_location_observation_is_zero_delay() {
        let mut history=WorldlineHistory::new(3);
        history.push(sample(0.0,0.0)).unwrap();
        history.push(sample(1.0,0.0)).unwrap();
        let o=history.observe(V2::default(),1.0).unwrap();
        check(o.emitted_at_s,1.0,1e-12);
    }
    #[test]
    fn emitter_in_future_is_unavailable() {
        let mut history=WorldlineHistory::new(4);
        history.push(sample(10.0,0.0)).unwrap();
        history.push(sample(11.0,0.0)).unwrap();
        assert!(history.observe(V2::default(),2.0).is_none());
    }
}
