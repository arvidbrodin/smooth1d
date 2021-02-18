/*
 * Copyright (c) 2021 Arvid Brodin
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

pub struct Segment {
	initvals: Vec<f64>,
	duration: f64,
	padto: usize,		// Return state of at least this length
}

impl Segment {
	pub fn new(initvals: &[f64], duration: f64, padto: usize) -> Self {
		let seg = Self {
			initvals: initvals.to_vec(),
			duration: duration,
			padto: padto,
		};
		seg.print();
		seg
	}

	/*
	 * Treat initvals as coefficients of the terms in integrals of some
	 * derivative of position - that is, given
	 *
	 * s0 = initvals[2] (position)
	 * v0 = initvals[1] (velocity)
	 * a0 = initvals[0] (acceleration)
	 *
	 * return state at time t according to
	 *
	 * s = 1/6*j0*t³ + 1/2*a0*t² + v0*t + s0
	 * v = 1/2*j0*t² + a0*t + v0
	 * a = j0*t + a0
	 * j = j0
	 */
	pub fn get_state_at(&self, t: f64) -> Vec<f64> {
		assert!(t >= 0.0);
		assert!(t <= self.duration);

		let mut terms = Vec::new();
		let mut state = vec![0.0; self.padto - self.initvals.len()];
		for initval in &self.initvals {
			let mut val = 0.0;
			let degree = terms.len();
			for n in 0..degree {
				terms[n] *= t/(degree - n) as f64;
				val += terms[n];
			}
			terms.push(*initval);
			val += *initval;
			state.push(val);
		}

		state
	}

	pub fn get_end_state(&self) -> Vec<f64> {
		self.get_state_at(self.duration)
	}

	pub fn get_duration(&self) -> f64 {
		self.duration
	}

	pub fn print(&self) {
		eprintln!("Segment: duration {}", self.duration);
		eprintln!("   Initvals: {:?}", self.initvals);
		eprintln!("   Endstate: {:?}", self.get_end_state());
	}
}
