module main

import math

struct Complex {
	real f64
	imag f64
}

fn (c Complex) abs() f64 {
	return math.sqrt(c.real * c.real + c.imag * c.imag)
}

[@inline]
fn (c Complex) + (o Complex) Complex {
	return Complex{c.real + o.real, c.imag + o.imag}
}

[@inline]
fn (c Complex) - (o Complex) Complex {
	return Complex{c.real - o.real, c.imag - o.imag}
}

[@inline]
fn (c Complex) * (o Complex) Complex {
	return Complex{c.real * o.real - c.imag * o.imag, c.real * o.imag + c.imag * o.real}
}

[@inline]
fn (c Complex) conjugate() Complex {
	return Complex{c.real, -c.imag}
}

[@inline]
fn (c Complex) str () string {
	return '${c.real} + ${c.imag}i'
}

fn (c Complex) / (o Complex) Complex {
	denom := o.real * o.real + o.imag * o.imag
	if denom == 0 {
		panic('Division by zero')
	}
	return Complex{(c.real * o.real + c.imag * o.imag) / denom, (c.imag * o.real - c.real * o.imag) / denom}
}

[@inline]
fn (c Complex) == (o Complex) bool {
	return c.real == o.real && c.imag == o.imag
}

fn (c Complex) angle() f64 {
	return math.atan2(c.imag, c.real)
}

fn (c Complex) pow(n int) Complex {
	if n == 0 {
		return Complex{1, 0}
	}
	mut res := c
	for i := 1; i < n; i++ {
		res = res * c
	}
	return res
}

fn (c Complex) angle_deg() f64 {
	return c.angle() * 180.0 / math.pi
}


[@inline]
fn Complex.new(real f64, imag f64) Complex {
	return Complex{real, imag}
}



fn main() {
	c1 := Complex{1, 2}

	c2 := Complex{3, 4}

	println((c1 + c2).angle() * 180 / math.pi)
}