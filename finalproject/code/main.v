// module main
import vsl.vcl
import os
import stbi
// v -d vsl_vcl_dlopencl run main.v

const height = 1020
const width = 720

const kernel_path = os.join_path(os.dir(@FILE), 'kernel/slime.cl')
const output = os.join_path(os.dir(@FILE), 'output')

/**
	run_kernel: Runs the kernel specified as const kernel_path
*/
fn run_kernel() {
	kernel_cont := os.read_file(kernel_path) or { exit(-1) }
	println(kernel_cont)

	mut device := vcl.get_default_device() or { panic('Could not generate device') }
	defer {
		device.release() or { panic(err) }
	}

	// rbg does not exist, intensity does...
	mut img := device.image_2d(.rgba, width: width, height: height) or {
		panic('Could not generate image')
	}

	defer {
		img.release() or { panic(err) }
	}

	device.add_program(kernel_cont) or { panic(err) }

	k := device.kernel('slime') or { panic(err) }

	kernel_err := <-k.global(int(img.bounds.width), int(img.bounds.height))
		.local(1, 1).run(img)

	println(typeof(kernel_err))

	if kernel_err !is none {
		panic(kernel_err)
	}

	buffer := img.data_2d() or { panic(err) }

	stbi.stbi_write_png(os.join_path(output, 'slime.png'), width, height, 4, buffer.data,
		0) or { panic(err) }
}

fn main() {
	println('Hello world!')
	println(kernel_path)
	run_kernel()
}
