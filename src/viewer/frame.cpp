#include "frame.h"


// geometry processing and shape analysis framework
namespace gproshan {


frame::frame()
{
	render_tex = 0;
	vbo = 0;

	program.load_vertex("../shaders/vertex_frame.glsl");	
	program.load_fragment("../shaders/fragment_frame.glsl");
	

	const GLfloat vertices[] = {
								-1, -1, 0,
								1, -1, 0,
								-1, 1, 0,
								-1, 1, 0,
								1, -1, 0,
								1, 1, 0,
								};

	glGenTextures(1, &render_tex);
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &pbo);
	

	glBindVertexArray(vao);

	glBindTexture(GL_TEXTURE_2D, render_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	glBindVertexArray(0);
}

frame::~frame()
{
	glDeleteBuffers(1, &pbo);
	glDeleteBuffers(1, &vbo);
	glDeleteVertexArrays(1, &vao);
	glDeleteTextures(1, &render_tex);
}

void frame::display(const int & width, const int & height, void * buffer)
{
	glBindBuffer(GL_ARRAY_BUFFER, pbo);
	glBufferData(GL_ARRAY_BUFFER, width * height * sizeof(float) * 4, buffer, GL_STREAM_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	program.enable();

	glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, render_tex);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
	glUniform1i(program("render_tex"), 0);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	
	program.disable();
}


} // namespace gproshan

