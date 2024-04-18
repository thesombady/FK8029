import math
import rand
import rand.seed
import rand.pcg32
import vtikz // Module that needs to be imported, v install --git www.github.com/thesombady/vtikz
import os
import vector { Vec }


const xmin = -10
const xmax = 10
const sep = 0.250
const res = 500

// --------------- QMC ----------------

/*
	generate the trail wave-function
	: v f64 - the curret position
	: alpha f64 - the current alpha value
	: lambda f64 - the current lambda value
	; f64 returns the wave-funciton at the two coordinates
*/
fn wave_function(v f64, alpha f64, lambda f64) f64 {
	return math.exp(-alpha * math.pow(v, 2.0))
}

/*
	laplacian computes the laplacian operator on the wave-function, given a specific set of coordinates
	Uses the three point stencil to compute the laplacian
	: v f64 - The the current position 
	: alpha f64 - The current alpha value
	: lambda f64 - The current lambda value
	; f64 - returns the laplacian on psi
*/
fn laplacian(v f64, alpha f64, lambda f64) f64 {

	// Three point stencil
	h := 0.00001

	wf := wave_function(v, alpha, lambda)

	dx1 := (wave_function(v + h, alpha, lambda) - 2 * wf +
		wave_function(v - h, alpha, lambda)) / (math.pow(h, 2.0))

	return (dx1) / wf // We skip by We skip by multiplying by the factor since we divide by it too
}

fn hamiltonian_on_psi(pos f64, alpha f64, lambda f64) f64 {
	return (alpha + math.pow(pos, 2.0) * (0.5 - 2 * math.pow(alpha, 2.0)))
}

fn random_pos(v f64, mut rng rand.PRNG) !f64 {
	x := rng.f64_in_range(-sep, sep)!
	if v + x < xmin || v + x > xmax  {
		println('Performing new position calculation')
		return random_pos(v, mut rng) // Calling again to ensure we're within boundaries
	}

	return v + x 
}

/*
	Sample the monte carlo simulation

*/
fn sample(pos f64, mut sim QMC) {
	mut energy := hamiltonian_on_psi(pos, sim.alpha, sim.lambda) // Local energy
	sim.energy << energy
	//sim.average_energy = (sim.average_energy* (iteration_number - 1) + energy ) / f64(iteration_number)
	//sim.variance_energy = (sim.variance_energy * (iteration_number - 1) + math.pow(energy, 2.0)) / f64(iteration_number) - math.pow(sim.average_energy, 2.0)
}

fn run(v f64, mut sim QMC, mut rng rand.PRNG) {
	simulation_time := 1_500_000 // Iterations (250_000)
	dx := (xmax - xmin) / f64(res - 1)

	mut pos_bin := []int{len: res, init: 0} // Two dimensional bin
	mut acc_rate := 0

	mut pos := v

	for ts := 0; ts < simulation_time; ts++ {
		new_pos := random_pos(pos, mut rng) or { panic('Could not generate random position') }
		// Should we generate totally random positions or pertubate the positions from the current location

		mut prob := wave_function(new_pos, sim.alpha, sim.lambda)
		prob /= wave_function(pos, sim.alpha, sim.lambda)
		prob = math.pow(prob, 2.0) // More or less |psi_new/psi_old|^2

		accepting := (rand.f64() < math.min(1.0, prob))

		if accepting {

			pos = new_pos

			

			mut idx_x := int(new_pos / dx)
			idx_x = if idx_x <= res / 2 { res / 2 - idx_x } else { idx_x + res / 2 }

			pos_bin[idx_x] += 1
			acc_rate++
		}	

		if ts > 5000 {
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
	alpha    f64
	average_energy   f64
	variance_energy f64
	pos      []int
	energy []f64
	var []f64
}

fn initialize() !f64 {
	return rand.f64_in_range(xmin, xmax)!
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

	delta := 2 / (3 - math.sqrt(5))

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

	run(pos, mut sim_a, mut rng) //

	run(pos, mut sim_d, mut rng)

	mut energy_alpha := os.create('harmonic_energy_alpha.dat') or { panic('Could not generate file') }

	mut variance_alpha := os.create('harmonic_variance_alpha.dat') or { panic('Could not generate file') }

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
			
		alpha2 = alpha1 + (alpha4 - alpha1) / delta

		alpha3 = alpha4 - (alpha4 - alpha1) / delta

		sim_a.alpha = alpha2
		sim_d.alpha = alpha3


		// Reset
		sim_a.average_energy = 0
		sim_d.average_energy = 0
		sim_a.energy = []
		sim_d.energy = []

		run(pos, mut sim_a, mut rng)

		run(pos, mut sim_d, mut rng)

	}
	// They have converged to a point
	//println(sim_a.pos)
	println(average(sim_a.energy))
	println(average(sim_d.energy))

	///*
	mut total := 0
	for i in 0..res {
		total += sim_a.pos[i]
	}
	norm_bin := []f64{len: res, init: sim_a.pos[index] / f64(total)}
	x := []f64{len: res, init: xmin + index * (xmax - xmin) / f64(res - 1)}

	mut tikz := vtikz.Tikz.new('Position', '$\\psi(x)$', 'Harmonic oscillator')
	tikz.scatter(x: x, y: norm_bin, color: .blue )
	//tikz.set_compiler(.lualatex)
	tikz.set_compiler(.@none)
	tikz.set_grid(true)
	tikz.plot('harmonic_oscillator.tex')
	println(norm_bin)
	//*/
}
