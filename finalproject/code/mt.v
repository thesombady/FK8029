import vtikz
import os

const alpha = 1.1
const beta = 0.5
const delta = 0.3
const gamma = 0.9
const eta = 0.1
const chi = 0.1

const methods = {
	'rk4':   runge_kutta_4
	'euler': euler
}

fn euler(x f64, y f64, z f64, dt f64) (f64, f64, f64) {
	return x + dt * dx(x, y, z), y + dt * dy(x, y, z), z + dt * dz(x, y, z)
}

fn runge_kutta_4(x f64, y f64, z f64, dt f64) (f64, f64, f64) {
	// k1
	kx_1 := dt * dx(x, y, z)
	ky_1 := dt * dy(x, y, z)
	kz_1 := dt * dz(x, y, z)
	// k2
	kx_2 := dt * dx(x + 0.5 * kx_1, y + 0.5 * ky_1, z + 0.5 * kz_1)
	ky_2 := dt * dy(x + 0.5 * kx_1, y + 0.5 * ky_1, z + 0.5 * kz_1)
	kz_2 := dt * dz(x + 0.5 * kx_1, y + 0.5 * ky_1, z + 0.5 * kz_1)
	// k3
	kx_3 := dt * dx(x + 0.5 * kx_2, y + 0.5 * ky_2, z + 0.5 * kz_2)
	ky_3 := dt * dy(x + 0.5 * kx_2, y + 0.5 * ky_2, z + 0.5 * kz_2)
	kz_3 := dt * dz(x + 0.5 * kx_2, y + 0.5 * ky_2, z + 0.5 * kz_2)
	// k4
	kx_4 := dt * dx(x + kx_3, y + ky_3, z + kz_3)
	ky_4 := dt * dy(x + kx_3, y + ky_3, z + kz_3)
	kz_4 := dt * dz(x + kx_3, y + ky_3, z + kz_3)

	return x + (kx_1 + 2 * (kx_2 + kx_3) + kx_4) / 6, y + (ky_1 + 2 * (ky_2 + ky_3) + ky_4) / 6,
		z + (kz_1 + 2 * (kz_2 + kz_3) + kz_4) / 6
}

@[inline]
fn dx(x f64, y f64, z f64) f64 {
	return alpha * x - beta * x * z
}

@[inline]
fn dy(x f64, y f64, z f64) f64 {
	return eta * x * y - chi * y * z
}

@[inline]
fn dz(x f64, y f64, z f64) f64 {
	return delta * x * z - gamma * z
}

fn solver(f fn (f64, f64, f64, f64) (f64, f64, f64)) ([]f64, []f64, []f64, []f64) {
	mut x, mut y, mut z := 8.0, 1.0, 5.0

	dt := 0.001
	mut t := 0.0
	t_end := 20.0

	mut x_array := []f64{}
	mut y_array := []f64{}
	mut z_array := []f64{}
	mut t_array := []f64{len: int(t_end / dt), init: index * dt}

	for (t < t_end) {
		x, y, z = f(x, y, z, dt)
		x_array << x
		y_array << y
		z_array << z
		t += dt
	}

	return x_array, y_array, z_array, t_array
}

fn main() {
	args := os.args[1..].clone()
	mut x_ := []f64{}
	mut y_ := []f64{}
	mut z_ := []f64{}
	mut t_ := []f64{}

	match args[0] {
		'rk4' {
			x_, y_, z_, t_ = solver(methods['rk4'])
		}
		'euler' {
			x_, y_, z_, t_ = solver(methods['euler'])
		}
		else {
			eprintln('Method not implemented')
			exit(0)
		}
	}
	mut file := os.create('3comp.dat') or { panic(err) }
	mut file2 := os.create('3phase.dat') or { panic(err) }
	defer {
		file.close()
		file2.close()
	}

	file.writeln('t\tx\ty\tz') or { panic(err) }
	file2.writeln('x\ty\tz') or { panic(err) }
	for i in 0 .. x_.len {
		file.writeln('${t_[i]}\t${x_[i]}\t${y_[i]}\t${z_[i]}') or { panic(err) }
		if (i % 10 == 0) {
			file2.writeln('${x_[i]}\t${y_[i]}\t${z_[i]}') or { panic(err) }
		}
	}
}
