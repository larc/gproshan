#ifndef FRAME_H
#define FRAME_H


#include "viewer/shader.h"

#include "viewer/include_opengl.h"


// geometry processing and shape analysis framework
namespace gproshan {


class frame
{
	private:
		GLuint render_tex = 0;
		GLuint vao = 0;
		GLuint vbo = 0;
		GLuint pbo = 0;

		size_t pbo_size = 0;

		shader program;

	public:
		frame(const size_t & width, const size_t & height);
		~frame();

		operator const GLuint & () const;

		void display(const int & width, const int & height);
};


} // namespace gproshan

#endif // FRAME_H

