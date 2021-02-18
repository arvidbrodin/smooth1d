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

use std::collections::VecDeque;
use crate::segment::Segment;
use crate::poly::Poly;

const CLOSE_ENOUGH: f64 = 1e-12;
const TINY_DURATION: f64 = 1e-12;

pub struct Path {
	limits: Vec<f64>,	// Highest-derivative first: limits[0] is jerk (or acc).
	segments: VecDeque<Segment>,
	time: f64,
	state: Vec<f64>,	// Highest-derivative first: state[0] is jerk (or acc).
	target: f64,		// Position target. Used to zero inaccuracies at end of move.
}

impl Path {
	/*
	 * Parameter limits are from highest order to lowest, excluding velocity limit
	 * (since velocity target is specific to each move): (MAX_)ACC, (JERK)
	 */
	pub fn new(mut limits: Vec<f64>) -> Self {
		// Assure get_state() can return pos, vel and acc
		assert!(limits.len() >= 1);
		// We can only handle 2nd- and 3rd-degree paths
		assert!(limits.len() <= 2);

		let degree = limits.len() + 1;
		limits.reverse();
		Self {
			limits: limits,
			segments: VecDeque::new(),
			time: 0.0,
			state: vec![0.0; degree + 1],
			target: 0.0,
		}
	}

	pub fn replan(&mut self, s_target: f64, v_limit: f64) {
		let mut limits = self.limits.clone();
		limits.push(v_limit);

		eprintln!("Path::replan(), state {:?}, s_target {}, limits {:?}", self.state, s_target, limits);
		assert!(v_limit > 0.0);

		self.time = 0.0;
		self.segments.clear();

		if self.limits.len() == 1 {
			// Acc-limited path
			self.calc_path_2(&limits, s_target);
		} else /* self.limits.len() == 2 */ {
			// Jerk-limited path
			self.calc_path_3(&limits, s_target);
		}

		self.target = s_target;
	}

	pub fn stop(&mut self) {
		eprintln!("Path::stop(), state {:?}", self.state);

		self.time = 0.0;
		self.segments.clear();

		if self.limits.len() == 1 {
			// Acc-limited path
			self.calc_path_1(0.0);
		} else /* self.limits.len() == 2 */ {
			// Jerk-limited path
			self.calc_path_2(&self.limits.clone(), 0.0);
		}

		if !self.segments.is_empty() {
			let end_state = self.segments.back().unwrap().get_end_state();

			// Accept whatever position we end up at as target
			self.target = *end_state.last().unwrap();
		} else {
			self.target = *self.state.last().unwrap();
		}
	}

	pub fn update(&mut self, dt: f64) {
		assert!(dt >= 0.0);
		if self.segments.is_empty() {
			// No movement planned - nothing to do
			return;
		}

		self.time += dt;
		while self.time > self.segments[0].get_duration() {
			let seg = self.segments.pop_front().unwrap();
			self.time -= seg.get_duration();
			if self.segments.is_empty() {
				// Zero out any accumulated inaccuracies
				self.state = vec![0.0; self.limits.len() + 1];
				self.state.push(self.target);

				break;
			}
		}

		if !self.segments.is_empty() {
			self.state = self.segments[0].get_state_at(self.time);
		}
	}

	pub fn get_state(&self) -> (f64, f64, f64) {
		let pos_index = self.state.len() - 1;
		(self.state[pos_index], self.state[pos_index - 1], self.state[pos_index - 2])
	}

	pub fn is_active(&self) -> bool {
		!self.segments.is_empty()
	}

	fn get_end_state(&self) -> Vec<f64> {
		if self.segments.is_empty() {
			return self.state.clone();
		}
		self.segments.back().unwrap().get_end_state()
	}

	fn calc_path_1(&mut self, v_target: f64) {
		let mut state = self.get_end_state();
		let v_diff = v_target - state[1];
		let a0 = v_diff.signum()*self.limits[0];
		let t0 = v_diff/a0;

		let degree = self.limits.len() + 1;
		if t0.abs() > TINY_DURATION {
			state[0] = a0;
			self.segments.push_back(Segment::new(&state[..], t0, degree + 1));
		}
	}

	fn calc_path_2(&mut self, limits: &Vec<f64>, s_target: f64) {
		let mut state = self.get_end_state();
		let s_diff = s_target - state[2];
		let v0 = state[1];

		let v1_target = s_diff.signum()*limits[1];
		let v1_diff = v1_target - v0;

//		println!("calc_path_2(): s_diff = {}; v1_target = {}", s_diff, v1_target);

		let a0 = v1_diff.signum()*limits[0];
		let mut t0 = v1_diff/a0;

		let mut a2 = -v1_target.signum()*limits[0];
		let mut t2 = -v1_target/a2;

		println!("t0 = {}; a0 = {}; t2 = {}; a2 = {}", t0, a0, t2, a2);
		let mut t1 = s_diff/v1_target + 0.5*v0.powi(2)/(a0*v1_target) - 0.5*v1_target/a0 + 0.5*v1_target/a2;

		if t1 < 0.0 {
			// Solve for t0 with t1 = 0 (v_target never reached)
			let x = v0/a0;
			let roots = Poly::solve_quadratic(1.0, 2.0*x, 0.5*v0*x/a0 - s_diff/a0);
			eprintln!("Roots: {:?}", roots);
			t0 = roots[1];
			t1 = 0.0;
			t2 = t0 + x;
			if roots[0] > 0.0 {
				/* If both roots are positive, that means this
				 is a shortened move that overshoots the new
				 target, i.e. it crosses the target during t0.
				 Thus s_diff will have effectively changed sign
				 during t0, so change sign of a2. */
				a2 = -a2;
			}
		}

		let degree = self.limits.len() + 1;
		if t0.abs() > TINY_DURATION {
			state[0] = a0;
			self.segments.push_back(Segment::new(&state[..], t0, degree + 1));
			state = self.segments.back().unwrap().get_end_state();
		}

		if t1.abs() > TINY_DURATION {
			state[1] = v1_target;
			self.segments.push_back(Segment::new(&state[1..], t1, degree + 1));
			state = self.segments.back().unwrap().get_end_state();
		}

		if t2.abs() > TINY_DURATION {
			state[0] = a2;
			self.segments.push_back(Segment::new(&state[..], t2, degree + 1));
		}

		// Check result
		if self.segments.is_empty() {
			assert!(self.state[1].abs() < CLOSE_ENOUGH);
			assert!((self.state[2] - s_target).abs() < CLOSE_ENOUGH, "s_target = {}; self.target = {}", s_target, self.state[2]);
			return;
		}

		let state = self.segments.back().unwrap().get_end_state();
		assert!(state[1].abs() < CLOSE_ENOUGH);
		assert!((s_target - state[2]).abs() < CLOSE_ENOUGH);
	}

	fn calc_path_3(&mut self, limits: &Vec<f64>, s_target: f64) {
		let s_diff = s_target - self.state[3];
		let v3_target = s_diff.signum()*limits[2];

		self.calc_path_2(limits, v3_target);
		let coast_index = self.segments.len();
		self.calc_path_2(limits, 0.0);

		let mut state = self.get_end_state();
		let t3 = (s_target - state[3])/v3_target;

		let degree = self.limits.len() + 1;
		if t3 >= 0.0 {
			while self.segments.len() > coast_index {
				self.segments.pop_back();
			}
			state = self.get_end_state();
			state[2] = v3_target;
			self.segments.push_back(Segment::new(&state[2..], t3, degree + 1));
			self.calc_path_2(limits, 0.0);

			return;
		}

		todo!();

		// let v = <something>;
		// self.calc_path_2(limits, v);
		// self.calc_path_2(limits, 0.0);
	}

	pub fn print(&self) {
		for seg in &self.segments {
			seg.print();
		}
	}
}


#[cfg(test)]
mod tests {
	use super::Path;
	use super::CLOSE_ENOUGH;
	use std::io::Write;

	enum ActionType {
		MoveTo((f64, f64)),
		CheckAcc(f64),
		CheckVel(f64),
		CheckPos(f64),
		CheckState((f64, f64, f64)),
		Stop,
		Done,
	}

	struct Action {
		t: f64,
		action: ActionType,
	}

	fn val_eq(val1: f64, val2: f64) -> bool {
		(val1 - val2).abs() < CLOSE_ENOUGH
	}

	fn check_eq(time: f64, val: f64, goal: f64) -> Result<(), String> {
		if val_eq(val, goal) {
			return Ok(());
		}

		Err(format!("Time {:.3}: value {} differs from {}", time, val, goal))
	}

	fn check_states_eq(time: f64, state: (f64, f64, f64), goal: (f64, f64, f64)) -> Result<(), String> {
		if val_eq(state.0, goal.0) && val_eq(state.1, goal.1) && val_eq(state.2, goal.2) {
			return Ok(());
		}

		Err(format!("Time {:.3}: state {:?} differs from {:?}", time, state, goal))
	}

	fn write_state(orig_state: &Vec<f64>, file: &mut std::fs::File) {
		for val in orig_state.iter().rev() {
			write!(file, "{:.6} ", *val).expect("Could not write to file");
		}

		// Print dummy jerk to keep gnuplot happy
		if orig_state.len() < 4 {
			write!(file, "{:.6} ", 0.0).expect("Could not write to file");
		}
	}

	fn run_test(limits: Vec<f64>, actions: &[Action], test_name: &str) -> Result<(), String> {
		let mut file = std::fs::File::create(test_name.to_owned() + ".data").expect("Cannot create data file");

		let mut replans = Vec::new();

		let mut result = Ok(());

		let dt = 0.001;
		let tolerance_fact = 1.01;

		let mut path = Path::new(limits.clone());
		let mut t = 0.0;
		let mut s_prev = 0.0;
		let mut v_prev = 0.0;
		let mut a_prev = 0.0;
		let mut action_index = 0;
		loop {
			// Check limits
			let state = path.get_state();
			let v = (state.0 - s_prev)/dt;
			let a = (v - v_prev)/dt;
			let j = (a - a_prev)/dt;

			// Velocity is variable so we don't have a reference to check against

			if result.is_ok() {
				if a.abs() > limits[0]*tolerance_fact {
					result = Err(format!("Time {}: acceleration ({}) over limit ({})", t, a, limits[0]));
				}
			}

			if limits.len() > 1 {
				if result.is_ok() {
					if j.abs() > limits[1]*tolerance_fact {
						result = Err(format!("Time {}: jerk ({}) over limit ({})", t, j, limits[1]));
					}
				}
			}

			s_prev = state.0;
			v_prev = v;
			a_prev = a;

			// Check specific values at specific times
			let action = &actions[action_index];
			if t + dt >= action.t {
				let t_frac = action.t - t;
				path.update(t_frac);
				match action.action {
					ActionType::MoveTo((x, v)) => {
						path.replan(x, v);
						replans.push(action.t);
					},
					ActionType::Stop => {
						path.stop();
						replans.push(action.t);
					},
					ActionType::CheckAcc(acc) => {
						if result.is_ok() {
							let state = path.get_state();
							result = check_eq(action.t, state.2, acc);
						}
					},
					ActionType::CheckVel(vel) => {
						if result.is_ok() {
							let state = path.get_state();
							result = check_eq(action.t, state.1, vel);
						}
					},
					ActionType::CheckPos(pos) => {
						if result.is_ok() {
							let state = path.get_state();
							result = check_eq(action.t, state.0, pos);
						}
					},
					ActionType::CheckState(goal) => {
						if result.is_ok() {
							let state = path.get_state();
							result = check_states_eq(action.t, state, goal);
						}
					},
					ActionType::Done => {
						if result.is_ok() {
							if path.is_active() {
								result = Err(format!("Time {:.3}: Planner still active", action.t));
							}
						}
						break;
					},
				}
				let t_frac = t + dt - action.t;
				path.update(t_frac);
				action_index += 1;
			} else {
				path.update(dt);
			}

			// Plot state to file
			write!(file, "{:.6} ", t).expect("Could not write to file");
			write_state(&path.state, &mut file);
			writeln!(file, "").expect("Could not write to file");

			t += dt;
		}

		writeln!(file, "").expect("Could not write to file");
		writeln!(file, "").expect("Could not write to file");
		for t in replans {
			write!(file, "{:.6} ", t).expect("Could not write to file");
		}
		writeln!(file, "").expect("Could not write to file");

		result
	}

	// ***
	// *** Test acceleration-limited motion ***
	// ***

	// Zero-distance move
	#[test]
	fn alim_null_move() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.00, MAX_VEL)) },
			Action { t: 0.01, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_null_move")
	}

	// Long move that reach all the limits
	#[test]
	fn alim_move_all_limits_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckAcc(MAX_ACC) },
			Action { t: 0.30, action: ActionType::CheckState((0.02, MAX_VEL, 0.0)) },
			Action { t: 0.60, action: ActionType::CheckState((0.04, 0.0, 0.0)) },
			Action { t: 0.61, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_move_all_limits_pos")
	}
	#[test]
	fn alim_move_all_limits_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckAcc(-MAX_ACC) },
			Action { t: 0.30, action: ActionType::CheckState((-0.02, -MAX_VEL, 0.0)) },
			Action { t: 0.60, action: ActionType::CheckState((-0.04, 0.0, 0.0)) },
			Action { t: 0.61, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_move_all_limits_neg")
	}

	// Short, separate moves that reach a_max but not v_max
	#[test]
	fn alim_move_no_vmax_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.02, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::CheckState((0.02, 0.0, 0.0)) },
			Action { t: 0.31, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_move_no_vmax_pos")
	}
	#[test]
	fn alim_move_no_vmax_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.02, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::CheckState((-0.02, 0.0, 0.0)) },
			Action { t: 0.31, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_move_no_vmax_neg")
	}

	// Plan a move, then plan another one that continues in exactly the same way
	#[test]
	fn alim_continued_move_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.60, action: ActionType::CheckState((0.04, 0.0, 0.0)) },
			Action { t: 0.61, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_continued_move_pos")
	}
	#[test]
	fn alim_continued_move_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.60, action: ActionType::CheckState((-0.04, 0.0, 0.0)) },
			Action { t: 0.61, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_continued_move_neg")
	}

	// Interrupted move, higher velocity
	#[test]
	fn alim_inc_vel_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.06, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.06, MAX_VEL*1.5)) },
			Action { t: 0.41, action: ActionType::CheckVel(MAX_VEL*1.5) },
			Action { t: 0.80, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_inc_vel_pos")
	}

	// Interrupted move, higher velocity, doesn't reach vmax
	#[test]
	fn alim_inc_vel_novmax_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.05, MAX_VEL*1.5)) },
			Action { t: 0.70, action: ActionType::CheckPos(0.05) },
			Action { t: 0.71, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_inc_vel_novmax_pos")
	}

	// Interrupted move, lower velocity
	#[test]
	fn alim_dec_vel_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.05, MAX_VEL*0.5)) },
			Action { t: 1.00, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_dec_vel_pos")
	}

	// Shortened move - new command at s = 0.020 overshoots to 0.030 before going back to 0.025
	#[test]
	fn alim_shortened_move_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.025, MAX_VEL*0.5)) },
			Action { t: 0.80, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_shortened_move_pos")
	}
	#[test]
	fn alim_shortened_move_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((-0.025, MAX_VEL*0.5)) },
			Action { t: 0.80, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_shortened_move_neg")
	}

	// Reversed move - new command at s = 0.020, v > 0, with target that needs v < 0.
	#[test]
	fn alim_reversed_move_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((0.010, MAX_VEL*0.5)) },
			Action { t: 1.01, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_reversed_move_pos")
	}
	#[test]
	fn alim_reversed_move_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.05, MAX_VEL)) },
			Action { t: 0.30, action: ActionType::MoveTo((-0.010, MAX_VEL*0.5)) },
			Action { t: 1.01, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_reversed_move_neg")
	}

	// Stop during positive velocity
	#[test]
	fn alim_stop_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.15, action: ActionType::Stop },
			Action { t: 0.35, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_stop_pos")
	}
	// Stop during negative velocity
	#[test]
	fn alim_stop_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.15, action: ActionType::Stop },
			Action { t: 0.35, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_stop_neg")
	}
	#[test]
	fn alim_stop_interrupted() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 0.15;
		let limits = vec![MAX_ACC];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.25, MAX_VEL)) },
			Action { t: 2.00, action: ActionType::MoveTo((-0.10, MAX_VEL)) },
			Action { t: 3.00, action: ActionType::Stop },
			Action { t: 4.00, action: ActionType::Done },
		];
		run_test(limits, &actions, "alim_stop_interrupted")
	}

	// ***
	// *** Test jerk-limited motion ***
	// ***

	// Long move that reach all the limits
	#[test]
	fn jlim_move_all_limits_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.15, action: ActionType::CheckAcc(MAX_ACC) },
			Action { t: 0.35, action: ActionType::CheckState((0.02, MAX_VEL, 0.0)) },
			Action { t: 0.70, action: ActionType::CheckState((0.04, 0.0, 0.0)) },
			Action { t: 0.71, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_all_limits_pos")
	}
	#[test]
	fn jlim_move_all_limits_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.15, action: ActionType::CheckAcc(-MAX_ACC) },
			Action { t: 0.35, action: ActionType::CheckState((-0.02, -MAX_VEL, 0.0)) },
			Action { t: 0.70, action: ActionType::CheckState((-0.04, 0.0, 0.0)) },
			Action { t: 0.71, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_all_limits_neg")
	}

	// Long, separate moves that reach v_max but not a_max
	#[test]
	fn jlim_move_no_amax_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 10.0;
		const JERK: f64 = 62.5;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.02, MAX_VEL)) },
			Action { t: 0.14, action: ActionType::CheckState((0.01, MAX_VEL, 0.0)) },
			Action { t: 0.28, action: ActionType::CheckState((0.02, 0.0, 0.0)) },
			Action { t: 0.29, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_amax_pos")
	}
	#[test]
	fn jlim_move_no_amax_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 10.0;
		const JERK: f64 = 62.5;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.02, MAX_VEL)) },
			Action { t: 0.14, action: ActionType::CheckState((0.01, MAX_VEL, 0.0)) },
			Action { t: 0.28, action: ActionType::CheckState((0.02, 0.0, 0.0)) },
			Action { t: 0.29, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_amax_neg")
	}

	// Short, separate moves that reach a_max but not v_max
	#[test]
	fn jlim_move_no_vmax_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		const JERK: f64 = 10.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.02, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckAcc(1.0) },
			Action { t: 0.20, action: ActionType::CheckState((0.01, 0.1, 0.0)) },
			Action { t: 0.40, action: ActionType::CheckState((0.02, 0.0, 0.0)) },
			Action { t: 0.41, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_vmax_pos")
	}
	#[test]
//	#[ignore]
	fn jlim_move_no_vmax_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		const JERK: f64 = 10.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.02, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckAcc(-1.0) },
			Action { t: 0.20, action: ActionType::CheckState((-0.01, -0.1, 0.0)) },
			Action { t: 0.40, action: ActionType::CheckState((-0.02, 0.0, 0.0)) },
			Action { t: 0.41, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_vmax_neg")
	}

	// Short, separate moves that reach neither v_max nor a_max
	#[test]
//	#[ignore]
	fn jlim_move_no_vmax_nor_amax_pos() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		const JERK: f64 = 10.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.0025, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckState((0.00125, 0.025, 0.0)) },
			Action { t: 0.20, action: ActionType::CheckState((0.0025, 0.0, 0.0)) },
			Action { t: 0.21, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_vmax_nor_amax_pos")
	}

	#[test]
//	#[ignore]
	fn jlim_move_no_vmax_nor_amax_neg() -> Result<(), String> {
		const MAX_VEL: f64 = 0.2;
		const MAX_ACC: f64 = 1.0;
		const JERK: f64 = 10.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((-0.0025, MAX_VEL)) },
			Action { t: 0.10, action: ActionType::CheckState((-0.00125, -0.025, 0.0)) },
			Action { t: 0.20, action: ActionType::CheckState((-0.0025, 0.0, 0.0)) },
			Action { t: 0.21, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_move_no_vmax_nor_amax_neg")
	}

	// Plan a move, then plan another one that continues in exactly the same way
	#[test]
	fn jlim_continued_move() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.35, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.70, action: ActionType::CheckState((0.04, 0.0, 0.0)) },
			Action { t: 0.71, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_continued_move")
	}

	// Stop from v_max
	#[test]
	fn jlim_stop_at_vmax() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.35, action: ActionType::Stop },
			Action { t: 0.65, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_stop_at_vmax")
	}

	// Stop from a_max (not yet at v_max)
	#[test]
	fn jlim_stop_at_amax() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.15, action: ActionType::Stop },
			Action { t: 0.51, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_stop_at_amax")
	}

	// Stop from positive a, v (neither yet at max)
	#[test]
	fn jlim_stop_not_at_limits() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.01, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.25, action: ActionType::Stop },
			Action { t: 0.62, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_stop_not_at_limits")
	}

	/*
	 * Stop from a = -a_max and v > 0 but small enough to pass v = 0
	 * before a = 0. This requires heading into negative v before stopping.
	 */
	#[test]
	fn jlim_stop_switch_v() -> Result<(), String> {
		const MAX_VEL: f64 = 0.1;
		const MAX_ACC: f64 = 0.5;
		const JERK: f64 = 5.0;
		let limits = vec![MAX_ACC, JERK];
		let actions = [
			Action { t: 0.00, action: ActionType::MoveTo((0.04, MAX_VEL)) },
			Action { t: 0.35, action: ActionType::MoveTo((-0.04, MAX_VEL)) },
			Action { t: 0.58, action: ActionType::Stop },
			Action { t: 0.80, action: ActionType::Done },
		];
		run_test(limits, &actions, "jlim_stop_switch_v")
	}


	// Interrupted moves at v_max: same v_max, same direction
	// Interrupted moves below v_max: same v_max, same direction
	// Interrupted moves at v_max: same v_max, other direction
	// Interrupted moves below v_max: same v_max, other direction
	// Interrupted moves at v_max: target too close (position overshoot)
	// Interrupted moves below v_max: target too close (position overshoot)

}
