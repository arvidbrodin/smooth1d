/*
 * This code ported from GNU Science Library, poly/solve_quadratic.c and
 * poly/solve_cubic.c, copyright (C) 1996, 1997, 1998, 1999, 2000
 * Brian Gough and released under GPLv3 or later.
 *
 * Ported code copyright (C) 2021 Arvid Brodin.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#![allow(non_snake_case)]

pub struct Poly;

impl Poly {
	/* GSL code written to return +1.0 if num == -0.0. Rust's num.signum()
	 * returns -1.0 if num == -0.0, so we can't use that. */
	fn sgn(num: f64) -> f64 {
		if num >= -0.0 {
			return 1.0;
		}
		-1.0
	}

	pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
		let mut res = Vec::new();

		// Handle linear case
		if a == 0.0 {
			if b == 0.0 {
				return res;
			} else {
				res.push(-c/b);
				return res;
			}
		}

		let disc = b.powi(2) - 4.0*a*c;

		if disc > 0.0 {
			if b == 0.0 {
				let r = (-c/a).sqrt();
				res.push(-r);
				res.push(r);
				return res;
			}

			let temp = -0.5*(b + Self::sgn(b)*disc.sqrt());
			let r1 = temp/a;
			let r2 = c/temp;

			if r1 < r2 {
				res.push(r1);
				res.push(r2);
			} else {
				res.push(r2);
				res.push(r1);
			}

			return res;
		}

		if disc == 0.0 {
			res.push(-0.5*b/a);
			res.push(-0.5*b/a);
			return res;
		}

		// Discriminant < 0.0; no roots
		return res;
	}

	fn gsl_poly_solve_cubic(a: f64, b: f64, c: f64) -> Vec<f64> {
		let q = a.powi(2) - 3.0*b;
		let r = 2.0*a.powi(3) - 9.0*a*b + 27.0*c;

		let Q = q/9.0;
		let R = r/54.0;

		let Q3 = Q.powi(3);
		let R2 = R.powi(2);

		let CR2 = 729.0*r.powi(2);
		let CQ3 = 2916.0*q.powi(3);

		let mut res = Vec::new();

		if R == 0.0 && Q == 0.0 {
			res.push(-a/3.0);
			res.push(-a/3.0);
			res.push(-a/3.0);
			return res;
		}

		if CR2 == CQ3 {
			/* this test is actually R2 == Q3, written in a form suitable
			   for exact computation with integers */

			/* Due to finite precision some double roots may be missed, and
			   considered to be a pair of complex roots z = x +/- epsilon i
			   close to the real axis. */

			let sqrtQ = Q.sqrt();

			if R > 0.0 {
				res.push(-2.0*sqrtQ - a/3.0);
				res.push(sqrtQ - a/3.0);
				res.push(sqrtQ - a/3.0);
			} else {
				res.push(-sqrtQ - a/3.0);
				res.push(-sqrtQ - a/3.0);
				res.push(2.0*sqrtQ - a/3.0);
			}
			return res;
		}

		if R2 < Q3 {
			let ratio = Self::sgn(R)*(R2/Q3).sqrt();
			let theta = ratio.acos();
			let norm = -2.0*Q.sqrt();
			res.push(norm*(theta/3.0).cos() - a/3.0);
			res.push(norm*((theta + 2.0*std::f64::consts::PI)/3.0).cos() - a/3.0);
			res.push(norm*((theta - 2.0*std::f64::consts::PI)/3.0).cos() - a/3.0);

			// Sort roots into increasing order
			res.sort_by(|a, b| a.partial_cmp(b).unwrap());

			return res;
		}

		let A = -Self::sgn(R)*(R.abs() + (R2 - Q3).sqrt()).powf(1.0/3.0);
		let B = Q/A;
		res.push(A + B - a/3.0);
		return res;
	}

	pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
		Self::gsl_poly_solve_cubic(b/a, c/a, d/a)
	}
}
