import math
import rand
import rand.seed
import rand.pcg32
import vtikz // Module that needs to be imported, v install --git www.github.com/thesombady/vtikz
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
	wave_function: generate the trail wave-function
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
	laplacian: computes the laplacian operator on the wave-function, given a specific set of coordinates
	Uses the three point stencil to compute the laplacian
	: pos []Vec - The position of the two particles
	: alpha f64 - The current alpha value
	: lambda f64 - The current lambda value
	; f64 - returns the laplacian on psi
*/
fn laplacian(pos []Vec, alpha f64, lambda f64) f64 {
	v := pos[0]
	u := pos[1]

	// Nine point stencil, central approximations
	h := 0.001
	wf := wave_function(v, u, alpha, lambda)
	x_1 := wave_function(v - Vec.new(h, 0), u, alpha, lambda)
	x1 := wave_function(v + Vec.new(h, 0), u, alpha, lambda)
	y_1 := wave_function(v - Vec.new(0, h), u, alpha, lambda)
	y1 := wave_function(v + Vec.new(0, h), u, alpha, lambda)
	x_2 := wave_function(v, u - Vec.new(h, 0), alpha, lambda)
	x2 := wave_function(v, u + Vec.new(h, 0), alpha, lambda)
	y_2 := wave_function(v, u - Vec.new(0, h), alpha, lambda)
	y2 := wave_function(v, u + Vec.new(0, h), alpha, lambda)

	return (x_1 + x1 + y_1 + y1 + x_2 + x2 + y_2 + y2 - 8 * wf) / (h * h * wf)
}

/*
	loc_energy: computes the local energy at the current positions
	: pos []Vec - The position of the two particles
	: alpha f64 - The current alpha value
	: lambda f64 - The current lambda value
	; f64 - returns the local energy 
*/
fn loc_energy(pos []Vec, alpha f64, lambda f64) f64 {
	mut value := 0.0

	value = -1.0 / 2.0 * laplacian(pos, alpha, lambda)

	v := pos[0]
	u := pos[1]

	value += (v.norm_squared() + u.norm_squared()) / 2.0

	value += lambda / (math.sqrt(math.pow(v.x - u.x, 2.0) + math.pow(v.y - u.y, 2.0)))

	return value
}

/*
	random_pos: generates a random position within the boundaries of the simulation
	The function is recurive and will call itself if the new position is outside the boundaries
	: v Vec - The current position
	: rng rand.PRNG - The random number generator
	; Vec - returns a new position
*/
fn random_pos(v Vec, mut rng rand.PRNG) !Vec {
	x := rng.f64_in_range(-sep, sep)!
	y := rng.f64_in_range(-sep, sep)!
	if v.x + x < xmin || v.x + x > xmax || v.y + y < ymin || v.y + y > ymax {
		println('Performing new position calculation')
		return random_pos(v, mut rng) // Calling again to ensure we're within boundaries
	}

	return Vec.new(v.x + x, v.y + y)
}

/*
	sample: samples the local energy, the probability distribution and the distance between the two particles
	: pos []Vec - The position of the two particles
	: sim QMC - The current simulation
*/
fn sample(pos []Vec, mut sim QMC) {
	mut energy := loc_energy(pos, sim.alpha, sim.lambda) // Local energy
	sim.energy << energy
	wf := wave_function(pos[0], pos[1], sim.alpha, sim.lambda)
	sim.prob_dist << wf * wf
	sim.r_12 << (pos[0] - pos[1]).norm()
}

/*
	run: Simulates via Monte Carlo (metroplis algorithm) of a system with two particles
	: positions []Vec - The initial positions of the two particles
	: mut sim QMC - The current simulation (Containing information such as alpha, lambda, average energy, variance energy, position and energy)
	: mut rng rand.PRNG - The random number generator
*/
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
			mut idx_x := int(new_pos.x / dx)
			idx_x = if idx_x <= res / 2 { res / 2 - idx_x } else { idx_x + res / 2 }

			mut idx_y := int(new_pos.y / dy)
			idx_y = if idx_y <= res / 2 { res / 2 - idx_y } else { idx_y + res / 2 }

			pos_bin[idx_x][idx_y] += 1
			acc_rate++
		}

		if ts > 5000 && ts % 10 == 0 {
			sample(pos, mut sim)
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
	alpha           f64
	average_energy  f64
	variance_energy f64
	pos             [][]int
	energy          []f64
	var             []f64
	prob_dist       []f64
	r_12            []f64
}

/*
	initialize: Initializes the positions of the two particles
	; []Vec - returns the positions of the two particles
*/
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

/*
	write_energy: Writes the energy to a file
	: file os.File - The file to write to
	: sim QMC - The current simulation
*/
fn write_energy(mut file os.File, sim QMC) {
	file.writeln('${sim.alpha}\t${sim.average_energy}') or { panic('Coult not write to file') }
}

/*
	write_variance: Writes the variance to a file
	: file os.File - The file to write to
	: sim QMC - The current simulation
*/
fn write_variance(mut file os.File, sim QMC) {
	file.writeln('${sim.alpha}\t${sim.variance_energy}') or { panic('Coult not write to file') }
}

/*
	write_dist: Writes the probability distribution to a file given a distance r_12
	: file os.File - The file to write to
	: sim QMC - The current simulation
*/
fn write_dist(mut file os.File, sim QMC) {
	for i := 0; i < sim.prob_dist.len / 10; i += 10 {
		file.writeln('${sim.r_12[i]}\t${sim.prob_dist[i]}') or { panic('Could not write to file') }
	}
}

/*
	average: Computes the average of a vector
	: v []f64 - The vector to compute the average of
	; f64 - returns the average
*/
fn average(v []f64) f64 {
	mut t := 0.0
	for i in 0 .. v.len {
		t += v[i]
	}
	return t / f64(v.len)
}

/*
	variance: Computes the variance of a vector
	: v []f64 - The vector to compute the variance of
	: average f64 - The average of the vector
	; f64 - returns the variance
*/
fn variance(v []f64, average f64) f64 {
	mut t := 0.0
	for i in 0 .. v.len {
		t += math.pow(v[i] - average, 2.0)
	}
	return t / f64(v.len)
}

/*
	main: The entry point of the program
*/
fn main() {
	mut rng := &rand.PRNG(pcg32.PCG32RNG{}) // random number generator

	rng.seed(seed.time_seed_array(pcg32.seed_len)) // Seed the generator

	mut pos := initialize() or { panic(err) } // Initialize the positions of the two particles

	mut lambda := 2.0 // The lambda parameter

	delta := 2 / (3 - math.sqrt(5)) // Golden ratio

	mut alpha1 := 0.1
	mut alpha4 := 1.1

	mut alpha2 := alpha1 + (alpha4 - alpha1) / delta

	mut alpha3 := alpha4 - (alpha4 - alpha1) / delta

	mut sim_a := QMC{
		alpha: alpha2 // Random value
		lambda: lambda
	}

	mut sim_d := QMC{
		alpha: alpha3 // initial guess
		lambda: lambda
	}

	run(pos, mut sim_a, mut rng) // run the simulation to th left

	run(pos, mut sim_d, mut rng) // run the simulation to the right

	// Create the files to store information
	mut energy_alpha := os.create('energy_alpha${int(lambda)}.dat') or {
		panic('Could not generate file')
	}

	mut variance_alpha := os.create('variance_alpha${int(lambda)}.dat') or {
		panic('Could not generate file')
	}

	mut prob_dist := os.create('prod_dist${int(lambda)}.dat') or {
		panic('Could not generate file')
	}

	energy_alpha.writeln('alpha\tenergy')!
	variance_alpha.writeln('alpha\tvariance')!
	prob_dist.writeln('r_12\tp')!

	defer {
		energy_alpha.close()
		variance_alpha.close()
		prob_dist.close()
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

		alpha2 = alpha1 + (alpha4 - alpha1) / delta

		alpha3 = alpha4 - (alpha4 - alpha1) / delta

		sim_a.alpha = alpha2
		sim_d.alpha = alpha3

		println('Alpha2: ${alpha2}')
		println('Alpha3: ${alpha3}')

		// Reset
		sim_a.average_energy = 0
		sim_d.average_energy = 0
		sim_a.prob_dist = []
		sim_d.prob_dist = []
		sim_a.r_12 = []
		sim_d.r_12 = []
		sim_a.energy = []
		sim_d.energy = []

		run(pos, mut sim_a, mut rng)

		run(pos, mut sim_d, mut rng)
	}
	// They have converged to a point
	write_dist(mut prob_dist, sim_a)
	println(sim_a.average_energy)

	// The below code is for plotting the wave-function, the package vtikz is needed,
	// but it's currently not up to date on github. If you want to run this code, you need to comment out the below code
	mut total := 0
	for i in 0..res {
		for j in 0..res {
			total += sim_a.pos[i][j]
		}
	}

	mut norm_bin := [][]f64{len: res, init: []f64{len: res, init: 0.0}}
	for i in 0..res {
		for j in 0..res {
			norm_bin[i][j] = f64(sim_a.pos[i][j]) / f64(total)
		}
	}

	mut tikz := vtikz.Tikz.new('x', 'y', '\\psi_T')
	tikz.scatter3[f64](data: norm_bin, xlim: [f64(xmin), xmax]!, ylim: [f64(ymin), ymax]!, type_: .surface)
	//tikz.set_compiler(.lualatex)
	tikz.set_compiler(.@none)
	//tikz.set_grid(true)
	tikz.plot('wave2.tex')

}
