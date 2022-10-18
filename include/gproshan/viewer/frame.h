#ifndef FRAME_H
#define FRAME_H


#include <gproshan/geometry/vec.h>
#include <gproshan/viewer/shader.h>
#include <gproshan/viewer/include_opengl.h>

#ifdef GPROSHAN_CUDA
	#include <cuda_runtime.h>
	#include <cuda_gl_interop.h>
#endif // GPROSHAN_CUDA


// geometry processing and shape analysis framework
namespace gproshan {


class frame
{
	private:
		GLuint render_tex = 0;
		GLuint vao = 0;
		GLuint vbo = 0;
		GLuint pbo = 0;

		size_t width = 0;
		size_t height = 0;

		shader program;

	#ifdef GPROSHAN_CUDA
		cudaGraphicsResource_t pbo_cuda;
	#endif // GPROSHAN_CUDA

	public:
		frame();
		~frame();

		operator const GLuint & () const;

		vec4 * map_pbo(bool cuda = false);
		void unmap_pbo(bool cuda = false);

		bool resize(const size_t & w, const size_t & h);
		void display();
};


} // namespace gproshan

#endif // FRAME_H

