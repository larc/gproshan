#ifndef FRAME_H
#define FRAME_H


#include "viewer/shader.h"

#include "viewer/include_opengl.h"


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

#endif // FRAME_H

