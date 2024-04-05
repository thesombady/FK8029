module main

import vsl.la
import vsl.vlas

const xmin = -6.0
const xmax = 6.0


fn potential(x f64) f64 {
	return x * x
}

/*
	* Create a matrix for the 1D Schrodinger equation
	@param[in] n The number of points
	@param[in] ve The potential energy function

	@return corresponding matrix
*/
fn create_matrix(n int, ve fn ( f64 ) f64 ) &la.Matrix[f64] {
	/*
		The equation:
		\[
			-( d^2/dz^2 + z^2)psi(z) = 2 E / (hbar omega) psi(z)
		\]	
	*/	
	mut hamiltonian := la.Matrix.new[f64](n, n)

	dx := (xmax - xmin) / f64(n - 1)
	println("dx: ${dx}")

	mut x := xmin

	for i in 0..n {
		if i == 0 {
			// We set [0, 1] seperatlly
			hamiltonian.set(i, 1, -1 / (dx * dx))
		} else if i == n - 1 {
			// We set [n-1, n-2] seperatlly
			hamiltonian.set(i, i - 1, -1 / (dx * dx))
		} else {
			// we set the off-set diagonal elements
			hamiltonian.set(i, i - 1, -1 / (dx * dx))
			hamiltonian.set(i, i + 1, -1 / (dx * dx))
		}

		hamiltonian.set(i, i, 2 / (dx * dx) + ve(x))
		x += dx

	}

	return hamiltonian
}

fn main() {
	// v -d vsl_vlas_lapacke run main.v to compile and run with lapacke backend use the following command
	// else; v run main.v
	n := 10

	matrix := create_matrix(n, potential)

	println(matrix)
}