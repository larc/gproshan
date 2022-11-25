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

	glGenBuffers(1, &tex_vbo);
	update_vbo_texcoords();
}

scene_viewer::~scene_viewer()
{
	glDeleteBuffers(1, &tex_vbo);
	glDeleteTextures(sc->textures.size(), gltextures);
	delete [] gltextures;
}

void scene_viewer::init_texture(const GLuint & gltex, const scene::texture & tex)
{
	glBindTexture(GL_TEXTURE_2D, gltex);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.cols, tex.rows, 0, GL_RGB, GL_UNSIGNED_BYTE, tex.data);
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

			gl_uniform_material(program, mat);
			glDrawArrays(GL_TRIANGLES, obj.begin, sc->objects[i + 1].begin - obj.begin);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, 0);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, 0);
		}
	}
	glBindVertexArray(0);

	program.disable();
}

void scene_viewer::gl_uniform_material(shader & program, const scene::material & mat)
{
	glProgramUniform3f(program, program("mat.Ka"), mat.Ka.x(), mat.Ka.y(), mat.Ka.z());
	glProgramUniform3f(program, program("mat.Kd"), mat.Kd.x(), mat.Kd.y(), mat.Kd.z());
	glProgramUniform3f(program, program("mat.Ks"), mat.Ks.x(), mat.Ks.y(), mat.Ks.z());
	glProgramUniform1f(program, program("mat.d"), mat.d);
	glProgramUniform1f(program, program("mat.Ns"), mat.Ns);
	glProgramUniform1f(program, program("mat.Ni"), mat.Ni);
	glProgramUniform1i(program, program("mat.illum"), mat.illum);
	glProgramUniform1i(program, program("mat.map_Ka"), mat.map_Ka);
	glProgramUniform1i(program, program("mat.map_Kd"), mat.map_Kd);
	glProgramUniform1i(program, program("mat.map_Ks"), mat.map_Ks);

	if(mat.map_Ka > -1)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, gltextures[mat.map_Ka]);
		glProgramUniform1i(program, program("tex_Ka"), 0);
	}

	if(mat.map_Kd > -1)
	{
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, gltextures[mat.map_Kd]);
		glProgramUniform1i(program, program("tex_Kd"), 1);
	}
}

void scene_viewer::update_vbo_texcoords()
{
	if(!sc->texcoords) return;
	gproshan_error(texcoords);
	glBindVertexArray(vao);

	// 6 TEXTURE COORDS
	glBindBuffer(GL_ARRAY_BUFFER, tex_vbo);
	glBufferData(GL_ARRAY_BUFFER, mesh->n_vertices * sizeof(vec2), sc->texcoords, GL_STATIC_DRAW);
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 2, GL_REAL, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}


} // namespace gproshan

