#ifndef FRAME_H
#define FRAME_H


#include <gproshan/viewer/shader.h>

#include <gproshan/viewer/include_opengl.h>


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

	public:
		frame();
		~frame();

		operator const GLuint & () const;

		bool resize(const size_t & w, const size_t & h);
		void display();
};


} // namespace gproshan

#endif // FRAME_H

