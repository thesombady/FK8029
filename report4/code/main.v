import math
import rand
// import vtikz // Module that needs to be imported, v install --git www.github.com/thesombady/vtikz

// ------------- Vector 2d implementation -------------

const xmin = -10
const xmax = 10
const ymin = -10
const ymax = 10
const sep = 1.0

struct Vec {
	x f64
	y f64
}

@[inline]
fn (v Vec) + (u Vec) Vec {
	return Vec{
		x: v.x + u.x
		y: v.y + u.y
	}
}

@[inline]
fn (v Vec) - (u Vec) Vec {
	return Vec{
		x: v.x - u.x
		y: v.y - u.y
	}
}

fn (v Vec) str() string {
	return '[${v.x}, ${v.y}]'
}

@[inline]
fn (v Vec) mul(s f64) Vec {
	return Vec{v.x * s, v.y * s}
}

@[inline]
fn Vec.new(x f64, y f64) Vec {
	return Vec{
		x: x
		y: y
	}
}

//  TODO: Fine to have as inline?
@[inline]
fn (v Vec) norm() f64 {
	return math.sqrt(v.x * v.x + v.y * v.y)
}

@[inline]
fn (v Vec) normalize() Vec {
	return v.mul(v.norm())
}

//  TODO: Fine to have as inline?
@[inline]
fn (v Vec) dot(u Vec) f64 {
	return v.x * u.x + v.y * u.y
}

// --------------- QMC ----------------

/*
	generate the trail wave-function
	: position []Vec - an array of all positions
	: alpha f64 - the current alpha value
	: lambda f64 - the current lambda value
	; f64 returns the wave-funciton at the two coordinates
*/
fn wave_function(v Vec, u Vec, alpha f64, lambda f64) f64 {
	mut value := 1.0

	value *= math.exp(-(v.dot(v)) / 2)
	value *= math.exp(-(u.dot(u)) / 2)

	f := math.sqrt(math.pow(v.x * v.x - u.x * u.x, 2) + math.pow(v.y * v.y - u.y * u.y, 2))

	value *= math.exp((lambda * f / (1 + alpha * f)))

	return value
}

/*
	laplacian computes the laplacian operator on the wave-function, given a specific set of coordinates
	: pos []Vec - The position of the two particles
	: alpha f64 - The current alpha value
	: lambda f64 - The current lambda value
	; f64 - returns the laplacian on psi
*/
fn laplacian(pos []Vec, alpha f64, lambda f64) f64 {
	v := pos[0]
	u := pos[1]
	f := math.sqrt(math.pow(v.x * v.x - u.x * u.x, 2) + math.pow(v.y * v.y - u.y * u.y, 2))
	mut value := 0.0

	const_f1 := lambda / (math.pow(1 + alpha * f, 2) * f)

	const_f2 := (1 + 3 * alpha * f) / (math.pow(1 + alpha * f, 3) * math.pow(f, 3))

	// x_1

	value += -1 + const_f1 - const_f2 * math.pow(v.x - u.x, 2) + math.pow(-v.x +
		const_f1 * (v.x - u.x), 2)

	// x_2
	value += -1 + const_f1 - const_f2 * math.pow(u.x - v.x, 2) + math.pow(-u.x +
		const_f1 * (u.x - v.x), 2)

	// y_1

	value += -1 + const_f1 - const_f2 * math.pow(v.y - u.y, 2) + math.pow(-v.y +
		const_f1 * (v.x - u.x), 2)
	// y_2

	value += -1 + const_f1 - const_f2 * math.pow(u.y - v.y, 2) + math.pow(-u.y +
		const_f1 * (u.x - v.x), 2)

	return value // We skip by We skip by multiplying by the factor since we divide by it too
}

fn hamiltonian_on_psi(pos []Vec, alpha f64, lambda f64) f64 {
	mut value := 0.0

	value = -1.0 / 2.0 * laplacian(pos, alpha, lambda)

	for i in 0 .. 2 {
		value += (pos[i].dot(pos[i])) / 2 // x_i^2 / 2 + y_i^2/2
	}

	v := pos[0]
	u := pos[1]
	value += lambda / (math.sqrt(math.pow(v.x * v.x - u.x * u.x, 2) +
		math.pow(v.y * v.y - u.y * u.y, 2)))

	return value
}

fn random_pos(v Vec) !Vec {
	x := rand.f64_in_range(-sep, sep)!
	y := rand.f64_in_range(-sep, sep)!
	if v.x + x < xmin || v.x + x > xmax || v.y + y < ymin || v.y + y > ymax {
		println('Performing new position calculation')
		return random_pos(v) // Calling again to ensure we're within boundaries
	}

	return Vec.new(v.x + x, v.y + y)
}

fn sample(pos []Vec, mut sim QMC, iteration_number int) {
	sim.particle1 << pos[0]
	sim.particle2 << pos[1]
	mut energy := hamiltonian_on_psi(pos, sim.alpha, sim.lambda)
	sim.energy = energy
	// dump(energy)
}

fn run(mut pos []Vec, mut sim QMC) {
	simulation_time := 10000 // Iterations (250_000)
	res := 100 // Resolution for bin (Used for plotting)
	dx := (xmax - xmin) / f64(res - 1)
	dy := (ymax - ymin) / f64(res - 1)

	mut pos_bin := [][]int{len: res, init: []int{len: res, init: 0}} // Two dimensional bin

	for ts := 0; ts < simulation_time; ts++ {
		idx := rand.intn(2) or { panic('Could not generator number') }
		ndx := if idx == 1 { 0 } else { 1 } // The other index, of the atom not moved

		new_pos := random_pos(pos[idx]) or { panic('Could not generate random position') }
		// Should we generate totally random positions or pertubate the positions from the current location

		mut prob := wave_function(new_pos, pos[ndx], sim.alpha, sim.lambda)
		prob /= wave_function(pos[0], pos[1], sim.alpha, sim.lambda)
		prob = math.pow(prob, 2.0) // More or less |psi_new/psi_old|^2

		accepting := (rand.f64() < math.min(1.0, prob))

		if accepting {
			mut idx_x := int(new_pos.x / dx)
			idx_x = if idx_x < 0 { res / 2 - idx_x } else { idx_x + res / 2 }

			mut idx_y := int(new_pos.y / dy)
			idx_y = if idx_y < 0 { res / 2 - idx_y } else { idx_y + res / 2 }

			pos_bin[idx_x][idx_y] += 1
			// println('x: ${idx_x}, y: ${idx_y}')
			pos[idx] = new_pos
		}

		if ts % 2500 == 0 && ts > 10 {
			// Save data
			// println('Saving data')
			sample(pos, mut sim, ts)
			// We want so save E_{alpha, L}
			// Variance and the bins visited
			// Here we want to sample
		}
	}
	sim.pos = pos_bin
}

struct QMC {
	lambda f64
mut:
	alpha     f64
	energy    f64
	variance  f64
	pos       [][]int
	particle1 []Vec
	particle2 []Vec
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

fn main() {
	rand.seed([u32(3223878742), 1732001562]) // Setting a seed

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

	mut pos1 := pos.clone()
	run(mut pos1, mut sim_a)

	mut pos2 := pos.clone()
	run(mut pos2, mut sim_d)

	// Here is where we do our goldi-search
	for math.abs(alpha2 - alpha3) > 1e-7 {
		if sim_a.energy > sim_d.energy {
			alpha1 = alpha2
		} else {
			alpha4 = alpha2
		}

		alpha2 = alpha4 - (alpha4 - alpha1) / delta

		alpha3 = alpha1 + (alpha4 - alpha1) / delta

		println(alpha2)
		println(alpha3)

		pos1 = pos.clone()
		pos2 = pos.clone()
		sim_a.alpha = alpha2
		sim_d.alpha = alpha3
		run(mut pos1, mut sim_a)

		run(mut pos2, mut sim_d)
	}
	// They have converged to a point
}
