module main

import vsl.la
import vsl.vlas
import os

const xmin = -6.0
const xmax = 6.0

fn potential(x f64) f64 {
	return x * x
}

fn flatten(m &la.Matrix[f64]) []f64 {
	mut flat := []f64{}
	for i in 0 .. m.m {
		flat << m.get_row(i)
	}
	return flat
}

enum Method {
	three_point
	five_point
	invalid
}

/*
	* Create a matrix for the 1D Schrodinger equation
	@param[in] n The number of points
	@param[in] ve The potential energy function

	@return corresponding matrix
*/
fn create_matrix(n int, ve fn (f64) f64, m Method) &la.Matrix[f64] {
	/*
		The equation:latex
    		\[
			1/2* ( -d^2/dz^2 + v(z) ) psi(z) = E / (hbar omega) psi(z)
		\]	
	*/
	mut hamiltonian := la.Matrix.new[f64](n, n)

	dx := (xmax - xmin) / f64(n - 1)
	println('dx: ${dx}')

	mut x := xmin
	if m == .three_point {
		c_1 := -1.0 / (dx * dx) / 2.0
		c1 := -1.0 / (dx * dx) / 2.0
		for i in 0 .. n {
			if i == 0 {
				// We set [0, 1] seperatlly
				hamiltonian.set(i, i + 1, c1)
			} else if i == n - 1 {
				// We set [n-1, n-2] seperatlly
				hamiltonian.set(i, i - 1, c_1)
			} else {
				// we set the off-set diagonal elements
				hamiltonian.set(i, i - 1, c_1)
				hamiltonian.set(i, i + 1, c1)
			}

			hamiltonian.set(i, i, (2.0 / (dx * dx) + ve(x)) / 2.0)
			x += dx
		}
	} else if m == .five_point {
		c_2 := 1.0 / (12 * dx * dx) / 2.0
		c_1 := -16.0 / (12 * dx * dx) / 2.0
		c1 := 1.0 / (12 * dx * dx) / 2.0
		c2 := -16.0 / (12 * dx * dx) / 2.0
		for i in 0 .. n {
			if i == 0 {
				hamiltonian.set(i, i + 1, c1)
				hamiltonian.set(i, i + 2, c2)
			} else if i == 1 {
				hamiltonian.set(i, i - 1, c_1)
				hamiltonian.set(i, i + 1, c1)
				hamiltonian.set(i, i + 2, c2)
			} else if i == n - 2 {
				hamiltonian.set(i, i - 2, c_2)
				hamiltonian.set(i, i - 1, c_1)
				hamiltonian.set(i, i + 1, c1)
			} else if i == n - 1 {
				hamiltonian.set(i, i - 2, c_2)
				hamiltonian.set(i, i - 1, c_1)
			} else {
				hamiltonian.set(i, i - 2, c_2)
				hamiltonian.set(i, i - 1, c_1)
				hamiltonian.set(i, i + 1, c1)
				hamiltonian.set(i, i + 2, c2)
			}

			hamiltonian.set(i, i, (30 / (12 * dx * dx) + ve(x)) / 2.0)
			x += dx
		}
	} else {
		panic('Invalid method')
	}
	return hamiltonian
}

fn main() {
	// v -d vsl_vlas_lapacke run main.v to compile and run with lapacke backend use the following command
	// else; v run main.v
	n := 20

	mut matrix := create_matrix(n, potential, Method.five_point)

	mut eigen_vectors := la.Matrix.new[f64](n, n)

	mut values := []f64{len: n}

	la.jacobi(mut eigen_vectors, mut values, mut matrix) or { panic(err) } // TODO: This method does not work?

	println(values)

	// Below is for trying another method
	// Now we flatten the matrix
	/*
	mut hamiltonian := flatten(matrix)
	mut wr := []f64{len: n} // real part of the eigenvalues
	mut wi := []f64{len: n} // imaginary part of the eigenvalues (is zero in our case)

	ldvl := 1 // no left eigenvalues
	mut vl := []f64{len: n * ldvl} // storage

	ldvr := n // n right eigenvalues
	mut vr := []f64{len: n * ldvr} // storage
	vlas.dgeev(false, true, n, mut hamiltonian, n, wr, wi, vl, ldvl, vr, ldvr) // compute

	println(hamiltonian)
	*/
}
