#include <gproshan/viewer/scene_viewer.h>


// geometry processing and shape analysis framework
namespace gproshan {


scene_viewer::scene_viewer(scene * p_sc): che_viewer(p_sc), sc(p_sc)
{
	gltextures = new GLuint[sc->textures.size()];

	glGenTextures(sc->textures.size(), gltextures);
	for(index_t i = 0; i < sc->textures.size(); ++i)
	{
		gproshan_log_var(sc->texture_name[i]);
		init_texture(gltextures[i], sc->textures[i]);
	}
}

scene_viewer::~scene_viewer()
{
	glDeleteTextures(sc->textures.size(), gltextures);
	delete [] gltextures;
}

void scene_viewer::init_texture(const GLuint & gltex, const scene::texture & tex)
{
	glBindTexture(GL_TEXTURE_2D, gltex);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.cols, tex.rows, 0, GL_RGB, GL_FLOAT, tex.data);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);
}

void scene_viewer::draw(shader & program)
{
	glProgramUniformMatrix4fv(program, program("model_mat"), 1, true, &model_mat[0][0]);
	glProgramUniform1ui(program, program("idx_colormap"), idx_colormap);
	glProgramUniform1i(program, program("render_flat"), render_flat);
	glProgramUniform1i(program, program("render_lines"), render_lines);
	glProgramUniform1i(program, program("render_wireframe"), render_triangles);

	glPolygonMode(GL_FRONT_AND_BACK, render_wireframe ? GL_LINE : GL_FILL);

	program.enable();

	glBindVertexArray(vao);
	if(sc->objects.size() == 1)
	{
		glDrawArrays(GL_TRIANGLES, 0, mesh->n_vertices);
	}
	else
	{
		for(index_t i = 0; i < sc->objects.size() - 1; ++i)
		{
			const scene::object & obj = sc->objects[i];
			const scene::material & mat = sc->materials[obj.material_id];
			glProgramUniform3f(program, program("mat.Kd"), mat.Kd.x(), mat.Kd.y(), mat.Kd.z());
			glDrawArrays(GL_TRIANGLES, obj.begin, sc->objects[i + 1].begin - obj.begin);
		}
	}
	glBindVertexArray(0);

	program.disable();
}


} // namespace gproshan

