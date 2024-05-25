import vtikz
import os

// import os

const alpha = 1.0
const beta = 0.1
const delta = 0.8
const gamma = 1.4

const methods = {
	'rk4':   runge_kutta_4
	'euler': euler
}

fn euler(x f64, y f64, dt f64) (f64, f64) {
	return x + dt * dx(x, y), y + dt * dy(x, y)
}

fn runge_kutta_4(x f64, y f64, dt f64) (f64, f64) {
	// k1
	kx_1 := dt * dx(x, y)
	ky_1 := dt * dy(x, y)
	// k2
	kx_2 := dt * dx(x + 0.5 * kx_1, y + 0.5 * ky_1)
	ky_2 := dt * dy(x + 0.5 * kx_1, y + 0.5 * ky_1)
	// k3
	kx_3 := dt * dx(x + 0.5 * kx_2, y + 0.5 * ky_2)
	ky_3 := dt * dy(x + 0.5 * kx_2, y + 0.5 * ky_2)
	// k4
	kx_4 := dt * dx(x + kx_3, y + ky_3)
	ky_4 := dt * dy(x + kx_3, y + ky_3)

	return x + (kx_1 + 2 * (kx_2 + kx_3) + kx_4) / 6, y + (ky_1 + 2 * (ky_2 + ky_3) + ky_4) / 6
}

@[inline]
fn dx(x f64, y f64) f64 {
	return alpha * x - beta * x * y
}

@[inline]
fn dy(x f64, y f64) f64 {
	return delta * x * y - gamma * y
}

fn solver(f fn (f64, f64, f64) (f64, f64)) ([]f64, []f64, []f64) {
	mut x, mut y := 8.0, 5.0

	dt := 0.001
	mut t := 0.0
	t_end := 20.0

	mut x_array := []f64{}
	mut y_array := []f64{}
	mut t_array := []f64{len: int(t_end / dt), init: index * dt}

	for (t < t_end) {
		x, y = f(x, y, dt)
		x_array << x
		y_array << y
		t += dt
	}

	return x_array, y_array, t_array
}

fn main() {
	args := os.args[1..].clone()
	mut x_ := []f64{}
	mut y_ := []f64{}
	mut t_ := []f64{}

	match args[0] {
		'rk4' {
			x_, y_, t_ = solver(methods['rk4'])
		}
		'euler' {
			x_, y_, t_ = solver(methods['euler'])
		}
		else {
			eprintln('Method not implemented')
			exit(0)
		}
	}

	mut pop := vtikz.Tikz.new('Time', 'Population', 'Population of the two species')

	pop.scatter(x: t_, y: x_, color: .red, legend: 'A')
	pop.scatter(x: t_, y: y_, color: .blue, legend: 'B')
	// pop.scatter(x: x_array, y: y_array, color: .magenta, legend: 'Phase')
	// pop.set_show_legends(false)

	pop.set_grid(true)
	pop.set_legend_pos(.north_east)
	pop.set_compiler(.lualatex)
	pop.plot('${args[0]}_pop.tex')

	mut phase := vtikz.Tikz.new('A', 'B', 'Phase diagram')

	phase.scatter(x: x_, y: y_, color: .red, legend: 'Pray')
	// phase.scatter(x: x_array, y: y_array, color: .magenta, legend: 'Phase')
	phase.set_show_legends(false)

	phase.set_grid(true)
	phase.set_legend_pos(.north_east)
	phase.set_compiler(.lualatex)
	phase.plot('${args[0]}_phase.tex')
}
