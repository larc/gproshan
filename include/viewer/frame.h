#ifndef VIEWER_FRAME_H
#define VIEWER_FRAME_H


#include "shader.h"

#include "include_opengl.h"


// geometry processing and shape analysis framework
namespace gproshan {


class frame
{
	private:
		GLuint render_tex;
		GLuint vao, vbo;
		
		shader program;

	public:
		frame();
		~frame();

		void display(const int & width, const int & height, void * buffer);
};


} // namespace gproshan

#endif // VIEWER_FRAME_H

