#include <gproshan/viewer/scene_viewer.h>


// geometry processing and shape analysis framework
namespace gproshan {


static const int gltex_nums[] = {GL_TEXTURE0, GL_TEXTURE1, GL_TEXTURE2, GL_TEXTURE3, GL_TEXTURE4};

scene_viewer::scene_viewer(scene * p_sc): che_viewer(p_sc), sc(p_sc)
{
	gltextures = new GLuint[sc->textures.size()];

	glGenTextures(sc->textures.size(), gltextures);
	for(index_t i = 0; i < sc->textures.size(); ++i)
		init_texture(gltextures[i], sc->textures[i]);

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
	tex.spectrum == 3 	? glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex.width, tex.height, 0, GL_RGB, GL_UNSIGNED_BYTE, tex.data)
	: tex.spectrum == 4	? glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.width, tex.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.data)
	: glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, tex.width, tex.height, 0, GL_RED, GL_UNSIGNED_BYTE, tex.data);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);
}

void scene_viewer::draw(shader & program)
{
	program.uniform("model_mat", model_mat);
	program.uniform("idx_colormap", idx_colormap);
	program.uniform("render_lines", render_lines);
	program.uniform("render_flat", render_flat);
	program.uniform("render_wireframe", render_triangles);

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

			for(auto & gltex_num: gltex_nums)
			{
				glActiveTexture(gltex_num);
				glBindTexture(GL_TEXTURE_2D, 0);
			}
		}
	}
	glBindVertexArray(0);

	program.disable();
}

void scene_viewer::draw_pointcloud(shader & program)
{
	program.uniform("model_mat", model_mat);
	program.uniform("idx_colormap", idx_colormap);
	program.uniform("render_lines", render_lines);
	program.uniform("point_normals", point_normals);
	program.uniform("point_size", point_size);

	program.enable();

	glBindVertexArray(vao);
	if(sc->objects.size() == 1)
	{
		glDrawArrays(GL_POINTS, 0, mesh->n_vertices);
	}
	else
	{
		for(index_t i = 0; i < sc->objects.size() - 1; ++i)
		{
			const scene::object & obj = sc->objects[i];
			const scene::material & mat = sc->materials[obj.material_id];

			gl_uniform_material(program, mat);
			glDrawArrays(GL_POINTS, obj.begin, sc->objects[i + 1].begin - obj.begin);

			for(auto & gltex_num: gltex_nums)
			{
				glActiveTexture(gltex_num);
				glBindTexture(GL_TEXTURE_2D, 0);
			}
		}
	}
	glBindVertexArray(0);

	program.disable();
}

void scene_viewer::gl_uniform_material(shader & program, const scene::material & mat)
{
	program.uniform("mat.Ka", mat.Ka);
	program.uniform("mat.Kd", mat.Kd);
	program.uniform("mat.Ks", mat.Ks);
	program.uniform("mat.d", mat.d);
	program.uniform("mat.Ns", mat.Ns);
	program.uniform("mat.Ni", mat.Ni);
	program.uniform("mat.illum", mat.illum);
	program.uniform("mat.map_Ka", mat.map_Ka);
	program.uniform("mat.map_Kd", mat.map_Kd);
	program.uniform("mat.map_Ks", mat.map_Ks);
	program.uniform("mat.map_d", mat.map_d);
	program.uniform("mat.map_bump", mat.map_bump);

	static auto bind_texture = [&](const int & i, const int & map, const char * tex)
	{
		if(map < 0) return;
		glActiveTexture(gltex_nums[i]);
		glBindTexture(GL_TEXTURE_2D, gltextures[map]);
		program.uniform(tex, i);
	};

	bind_texture(0, mat.map_Ka, "tex_Ka");
	bind_texture(1, mat.map_Kd, "tex_Kd");
	bind_texture(2, mat.map_Ks, "tex_Ks");
	bind_texture(3, mat.map_d, "tex_d");
	bind_texture(4, mat.map_bump, "tex_bump");
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

