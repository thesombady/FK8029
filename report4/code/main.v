import math
import rand
import rand.seed
import rand.pcg32
// import vtikz // Module that needs to be imported, v install --git www.github.com/thesombady/vtikz
import os
import vector { Vec } // Vector struct that is used to represent the position of the particles, includes new, norm, norm_squared
// and is operator overloaded


const xmin = -10
const xmax = 10
const ymin = -10
const ymax = 10
const sep = 0.250
const res = 100

// --------------- QMC ----------------

/*
	generate the trail wave-function
	: position []Vec - an array of all positions
	: alpha f64 - the current alpha value
	: lambda f64 - the current lambda value
	; f64 returns the wave-funciton at the two coordinates
*/
fn wave_function(v Vec, u Vec, alpha f64, lambda f64) f64 {
	mut value := math.exp(-(v.norm_squared() + u.norm_squared()) / 2.0)

	f := math.sqrt(math.pow(v.x - u.x, 2.0) + math.pow(v.y - u.y, 2.0))

	value *= math.exp((lambda * f / (1 + alpha * f)))

	return value
}

/*
	laplacian computes the laplacian operator on the wave-function, given a specific set of coordinates
	Uses the three point stencil to compute the laplacian
	: pos []Vec - The position of the two particles
	: alpha f64 - The current alpha value
	: lambda f64 - The current lambda value
	; f64 - returns the laplacian on psi
*/
fn laplacian(pos []Vec, alpha f64, lambda f64) f64 {
	v := pos[0]
	u := pos[1]

 	/*
	// Three point stencil
	h := 0.00001

	wf := wave_function(v, u, alpha, lambda)

	dx1 := (wave_function(v + Vec.new(h, 0), u, alpha, lambda) - 2 * wf +
		wave_function(v - Vec.new(h, 0), u, alpha, lambda))

	dy1 := (wave_function(v + Vec.new(0, h), u, alpha, lambda) - 2 * wf +
		wave_function(v - Vec.new(0, h), u, alpha, lambda))

	dx2 := (wave_function(v, u + Vec.new(h, 0), alpha, lambda) - 2 * wf + 
		wave_function(v, u + Vec.new(h, 0), alpha, lambda))

	dy2 := (wave_function(v, u + Vec.new(0, h), alpha, lambda) - 2 * wf +
	 	wave_function(v, u + Vec.new(0, h), alpha, lambda))

	return (dx1 + dy1 + dx2 + dy2) / (wf * h * h) // We skip by We skip by multiplying by the factor since we divide by it too
	*/
	r_ij := math.sqrt(math.pow(v.x - u.x, 2.0) + math.pow(v.y - u.y, 2.0))

	beta := lambda / (1 + alpha * r_ij)
	x12 := v.x - u.x
	y12 := v.y - u.y

	mut t := math.pow(beta * x12 / r_ij - v.x, 2.0)
	t += math.pow(beta * y12 / r_ij - v.y, 2.0)
	t += math.pow(beta * (-x12) / r_ij - u.x, 2.0)
	t += math.pow(beta * (-y12) / r_ij - u.y, 2.0)

	t += 2 * (beta * beta / (r_ij * lambda)) - 4 * alpha * beta * beta * beta / (lambda * lambda) -4

	return t / wave_function(v, u, alpha, lambda)
}

fn hamiltonian_on_psi(pos []Vec, alpha f64, lambda f64) f64 {
	mut value := 0.0

	value = -1.0 / 2.0 * laplacian(pos, alpha, lambda)

	// for i in 0 .. 2 {
	// 	value += (pos[i].dot(pos[i])) / 2 // x_i^2 / 2 + y_i^2/2
	// }
	v := pos[0]
	u := pos[1]

	value += ( v.norm_squared() + u.norm_squared() ) / 2.0

	value += lambda / (math.sqrt(math.pow(v.x - u.x, 2.0) + math.pow(v.y - u.y, 2.0)))

	return value
}

fn random_pos(v Vec, mut rng rand.PRNG) !Vec {
	x := rng.f64_in_range(-sep, sep)!
	y := rng.f64_in_range(-sep, sep)!
	if v.x + x < xmin || v.x + x > xmax || v.y + y < ymin || v.y + y > ymax {
		println('Performing new position calculation')
		return random_pos(v, mut rng) // Calling again to ensure we're within boundaries
	}

	return Vec.new(v.x + x, v.y + y)
}

fn sample(pos []Vec, mut sim QMC, iteration_number int) {
	mut energy := hamiltonian_on_psi(pos, sim.alpha, sim.lambda) // Local energy
	sim.energy << energy
}

fn run(positions []Vec, mut sim QMC, mut rng rand.PRNG) {
	simulation_time := 1_500_000 // Iterations (250_000)
	dx := (xmax - xmin) / f64(res - 1)
	dy := (ymax - ymin) / f64(res - 1)

	mut pos_bin := [][]int{len: res, init: []int{len: res, init: 0}} // Two dimensional bin
	mut acc_rate := 0

	mut pos := positions.clone()

	for ts := 0; ts < simulation_time; ts++ {
		idx := rng.intn(2) or { panic('Could not generator number') }
		ndx := if idx == 1 { 0 } else { 1 } // The other index, of the atom not moved

		new_pos := random_pos(pos[idx], mut rng) or { panic('Could not generate random position') }
		// Should we generate totally random positions or pertubate the positions from the current location

		mut prob := wave_function(new_pos, pos[ndx], sim.alpha, sim.lambda)
		prob /= wave_function(pos[0], pos[1], sim.alpha, sim.lambda)
		prob = math.pow(prob, 2.0) // More or less |psi_new/psi_old|^2

		accepting := (rng.f64() < math.min(1.0, prob))

		if accepting {

			pos[idx] = new_pos
			if ts > 5000 && ts % 10 == 0 {
				sample(pos, mut sim, ts)
			}

			mut idx_x := int(new_pos.x / dx)
			idx_x = if idx_x < 0 { res / 2 - idx_x } else { idx_x + res / 2 }

			mut idx_y := int(new_pos.y / dy)
			idx_y = if idx_y < 0 { res / 2 - idx_y } else { idx_y + res / 2 }

			pos_bin[idx_x][idx_y] += 1
			acc_rate++
		}

	}
	println('Acceptance rate: ${acc_rate / f64(simulation_time)}')
	sim.pos = pos_bin
	sim.average_energy = average(sim.energy)
	sim.variance_energy = variance(sim.energy, sim.average_energy)
}

struct QMC {
	lambda f64
mut:
	alpha    f64
	average_energy   f64
	variance_energy f64
	pos      [][]int
	energy []f64
	var []f64
}

fn initialize() ![]Vec {
	mut pos := []Vec{len: 2}
	pos[0] = Vec{
		x: rand.f64_in_range(xmin, xmax)!
		y: rand.f64_in_range(ymin, ymax)!
	}

	pos[1] = Vec{
		x: rand.f64_in_range(xmin, xmax)!
		y: rand.f64_in_range(ymin, ymax)!
	}

	return pos
}

fn write_energy(mut file os.File, sim QMC) {
	file.writeln('${sim.alpha}\t${sim.average_energy}') or { panic('Coult not write to file') }
}

fn write_variance(mut file os.File, sim QMC) {
	file.writeln('${sim.alpha}\t${sim.variance_energy}') or { panic('Coult not write to file') }
}

fn average(v []f64) f64 {
	mut t := 0.0
	for i in 0..v.len {
		t += v[i]
	}
	return t / f64(v.len)
}

fn variance(v []f64, average f64) f64 {
	mut t := 0.0
	for i in 0..v.len {
		t += math.pow(v[i] - average, 2.0)
	}
	return t / f64(v.len)
}
fn main() {
	// rand.seed([u32(3223878742), 1732001562]) // Setting a seed
	mut rng := &rand.PRNG(pcg32.PCG32RNG{})

	rng.seed(seed.time_seed_array(pcg32.seed_len))

	mut pos := initialize() or { panic(err) }

	mut lambda := 2.0

	delta := 1.618033988749//2 / (3 - math.sqrt(5))

	mut alpha1 := 0.1
	mut alpha4 := 1.1

	/*
	mut alpha2 := alpha1 + (alpha4 - alpha1) / delta

	mut alpha3 := alpha4 - (alpha4 - alpha1) / delta
	*/

	mut alpha2 := alpha4 - (alpha4 - alpha1) / delta

	mut alpha3 := alpha1 + (alpha4 - alpha1) / delta

	mut sim_a := QMC{
		alpha: alpha2 // Random value
		lambda: lambda
	}

	mut sim_d := QMC{
		alpha: alpha3 // initial guess
		lambda: lambda
	}

	run(pos, mut sim_a, mut rng) //

	run(pos, mut sim_d, mut rng)

	mut energy_alpha := os.create('energy_alpha.dat') or { panic('Could not generate file') }

	mut variance_alpha := os.create('variance_alpha.dat') or { panic('Could not generate file') }

	energy_alpha.writeln('alpha\tenergy')!
	variance_alpha.writeln('alpha\tvariance')!
	defer {
		energy_alpha.close()
		variance_alpha.close()
	}

	// Here is where we do our goldi-search
	for math.abs(alpha2 - alpha3) > 1e-3 {
		// Check which way to go

		if sim_a.average_energy > sim_d.average_energy {
			alpha1 = alpha2
			write_energy(mut energy_alpha, sim_a)
			write_variance(mut variance_alpha, sim_a)
		} else {
			alpha4 = alpha3
			write_energy(mut energy_alpha, sim_d)
			write_variance(mut variance_alpha, sim_d)
		}

		/*
		alpha2 = alpha1 + (alpha4 - alpha1) / delta

		alpha3 = alpha4 - (alpha4 - alpha1) / delta
		*/
		
		alpha2 = alpha4 - (alpha4 - alpha1) / delta

		alpha3 = alpha1 + (alpha4 - alpha1) / delta
			
		sim_a.alpha = alpha2
		sim_d.alpha = alpha3


		println('Alpha2: ${alpha2}')
		println('Alpha3: ${alpha3}')


		// Reset
		sim_a.average_energy = 0
		sim_d.average_energy = 0
		sim_a.energy = []
		sim_d.energy = []
		run(pos, mut sim_a, mut rng)

		run(pos, mut sim_d, mut rng)

	}
	// They have converged to a point
	println(sim_a.average_energy)
}
