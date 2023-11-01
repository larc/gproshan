#include <gproshan/viewer/frame.h>

#include <gproshan/include.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan {


frame::frame()
{
	program.load_vertex(shaders_path("vertex_frame.glsl"));
	program.load_fragment(shaders_path("fragment_frame.glsl"));


	const GLfloat vertices[] = {
								-1, -1, 0,
								1, -1, 0,
								-1, 1, 0,
								-1, 1, 0,
								1, -1, 0,
								1, 1, 0,
								};


	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glGenBuffers(1, &pbo);

	glGenTextures(1, &render_tex);
	glBindTexture(GL_TEXTURE_2D, render_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);
}

frame::~frame()
{
#ifdef GPROSHAN_CUDA
	if(width && height)
		cudaGraphicsUnregisterResource(pbo_cuda);
#endif // GPROSHAN_CUDA

	glDeleteBuffers(1, &vbo);
	glDeleteVertexArrays(1, &vao);
	glDeleteTextures(1, &render_tex);
	glDeleteBuffers(1, &pbo);
}

frame::operator const GLuint & () const
{
	return pbo;
}

vec4 * frame::map_pbo(bool cuda)
{
	if(!cuda)
	{
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
		return (vec4 *) glMapBuffer(GL_PIXEL_UNPACK_BUFFER, GL_READ_WRITE);
	}

	vec4 * img = nullptr;
	#ifdef GPROSHAN_CUDA
		size_t num_bytes = 0;
		cudaGraphicsMapResources(1, &pbo_cuda, 0);
		cudaGraphicsResourceGetMappedPointer((void **) &img, &num_bytes, pbo_cuda);
	#endif // GPROSHAN_CUDA
	return img;
}

void frame::unmap_pbo(bool cuda)
{
	if(!cuda)
	{
		glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
		return;
	}

	#ifdef GPROSHAN_CUDA
		cudaGraphicsUnmapResources(1, &pbo_cuda, 0);
	#endif // GPROSHAN_CUDA
}

bool frame::resize(const size_t & w, const size_t & h)
{
	if(w == width && height == h)
		return false;

	if(w * h > width * height)
	{
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
		glBufferData(GL_PIXEL_UNPACK_BUFFER, 4 * sizeof(float) * w * h, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
	#ifdef GPROSHAN_CUDA
		if(width && height)
			cudaGraphicsUnregisterResource(pbo_cuda);
		cudaGraphicsGLRegisterBuffer(&pbo_cuda, pbo, cudaGraphicsMapFlagsNone);
	#endif // GPROSHAN_CUDA
	}

	width = w;
	height = h;

	return true;
}

void frame::display()
{
	program.enable();

	glBindVertexArray(vao);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, render_tex);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, 0);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glDrawArrays(GL_TRIANGLES, 0, 6);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	glBindVertexArray(0);

	program.disable();
}


} // namespace gproshan

