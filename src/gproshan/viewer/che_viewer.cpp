#include <gproshan/viewer/che_viewer.h>

#include <cassert>
#include <cmath>
#include <numeric>
#include <algorithm>

#include <glm/gtc/matrix_transform.hpp>

#include <gproshan/raytracing/rt_embree.h>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_viewer::~che_viewer()
{
	glDeleteBuffers(5, vbo);
	glDeleteVertexArrays(1, &vao);
}

che *& che_viewer::operator -> ()
{
	return mesh;
}

che *const & che_viewer::operator -> () const
{
	return mesh;
}

che_viewer::operator che *& ()
{
	return mesh;
}

void che_viewer::init(che * m, const bool & center)
{
	glGenVertexArrays(1, &vao);
	glGenBuffers(6, vbo);

	mesh = m;
	center_mesh = center;

	update();
}

void che_viewer::update()
{
	model_mat = glm::mat4(1);
	if(center_mesh)
	{
		vertex pmin(INFINITY, INFINITY, INFINITY);
		vertex pmax(0, 0, 0);

		for(index_t v = 0; v < mesh->n_vertices; ++v)
		{
			const vertex & p = mesh->gt(v);

			pmin.x = min(pmin.x, p.x);
			pmin.y = min(pmin.y, p.y);
			pmin.z = min(pmin.z, p.z);

			pmax.x = max(pmax.x, p.x);
			pmax.y = max(pmax.y, p.y);
			pmax.z = max(pmax.z, p.z);
		}

		scale(2.0 / std::max({pmax.x - pmin.x, pmax.y - pmin.y, pmax.z - pmin.z}));
		translate(- (pmax + pmin) / 2);
	}

	render_pointcloud = mesh->is_pointcloud();
	selected_xyz.clear();
	update_vbo();

	delete pick_vertex;
	pick_vertex = new rt::embree({mesh}, {model_mat});
}

void che_viewer::update_vbo()
{
	update_vbo_geometry();
	update_vbo_normal();
	update_vbo_color();
	update_vbo_heatmap();
}

void che_viewer::update_vbo_geometry()
{
	glBindVertexArray(vao);

	// 0 VERTEX
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, mesh->n_vertices * sizeof(vertex), &mesh->gt(0), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 4 INDEXES
	if(!mesh->is_pointcloud())
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[4]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->n_half_edges * sizeof(index_t), &mesh->vt(0), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	glBindVertexArray(0);
}

void che_viewer::update_vbo_normal(const vertex * vnormal)
{
	if(!vnormal) vnormal = &mesh->normal(0);

	glBindVertexArray(vao);

	// 1 NORMAL
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, mesh->n_vertices * sizeof(vertex), vnormal, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_REAL, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void che_viewer::update_vbo_color(const che::rgb_t * vcolor)
{
	if(!vcolor) vcolor = &mesh->rgb(0);

	glBindVertexArray(vao);

	// 2 COLOR
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, mesh->n_vertices * sizeof(che::rgb_t), vcolor, GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void che_viewer::update_vbo_heatmap(const real_t * vheatmap)
{
	if(!vheatmap) vheatmap = &mesh->heatmap(0);

	glBindVertexArray(vao);

	// 3 HEAT MAP
	glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
	glBufferData(GL_ARRAY_BUFFER, mesh->n_vertices * sizeof(real_t), vheatmap, GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 1, GL_REAL, GL_TRUE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void che_viewer::update_instances_positions(const vector<vertex> & translations)
{
	n_instances = translations.size();
	if(!n_instances) return;

	glBindVertexArray(vao);

	// 4 INSTANCES (translations)
	glBindBuffer(GL_ARRAY_BUFFER, vbo[5]);
	glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(vertex), translations.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_REAL, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glVertexAttribDivisor(3, 1);

	glBindVertexArray(0);
}

void che_viewer::draw(shader & program)
{
	glProgramUniformMatrix4fv(program, program("model_mat"), 1, 0, &model_mat[0][0]);
	glProgramUniform1ui(program, program("idx_colormap"), idx_colormap);
	glProgramUniform1i(program, program("render_flat"), render_flat);
	glProgramUniform1i(program, program("render_lines"), render_lines);
	glProgramUniform1i(program, program("render_wireframe"), render_triangles);

	glPolygonMode(GL_FRONT_AND_BACK, render_wireframe ? GL_LINE : GL_FILL);

	program.enable();

	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[4]);

	if(n_instances)
		glDrawElementsInstanced(GL_TRIANGLES, mesh->n_half_edges, GL_UNSIGNED_INT, 0, n_instances);
	else
		glDrawElements(GL_TRIANGLES, mesh->n_half_edges, GL_UNSIGNED_INT, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	program.disable();
}

void che_viewer::draw_point_cloud(shader & program)
{
	glProgramUniformMatrix4fv(program, program("model_mat"), 1, 0, &model_mat[0][0]);
	glProgramUniform1ui(program, program("idx_colormap"), idx_colormap);
	glProgramUniform1i(program, program("render_lines"), render_lines);
	glProgramUniform1i(program, program("point_normals"), point_normals);
	glProgramUniform1ui(program, program("point_size"), point_size);

	program.enable();

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);

	glDrawArrays(GL_POINTS, 0, mesh->n_vertices);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	program.disable();
}

void che_viewer::draw_selected_vertices(che_viewer & sphere, shader & program)
{
	if(selected_xyz.size() != selected.size())
	{
		selected_xyz.clear();
		selected_xyz.reserve(selected.size());

		for(const index_t & v: selected)
			selected_xyz.push_back(mesh->gt(v));
	}

	if(selected_xyz.size())
	{
		sphere.model_mat = model_mat;
		sphere.update_instances_positions(selected_xyz);
		sphere.draw(program);
	}
}

void che_viewer::translate(const vertex & p)
{
	model_mat = glm::translate(model_mat, glm_vec3(p));
}

void che_viewer::scale(const real_t & s)
{
	model_mat = glm::scale(model_mat, {s, s, s});
}

void che_viewer::select(const index_t & x, const index_t & y, const glm::uvec2 & windows_size, const glm::mat4 & proj_view_mat, const glm::vec3 & cam_pos)
{
	glm::vec3 dir = pick_vertex->ray_view_dir(x, windows_size.y - y, windows_size, glm::inverse(proj_view_mat), cam_pos);
	index_t v = pick_vertex->closest_vertex(cam_pos, dir);
	if(v != NIL) selected.push_back(v);
}

void che_viewer::log_info()
{
	if(!mesh) return;

	gproshan_log_var(mesh->filename);
	gproshan_log_var(mesh->n_vertices);
	gproshan_log_var(mesh->n_faces);
	gproshan_log_var(mesh->n_half_edges);
	gproshan_log_var(mesh->n_edges);
	gproshan_log_var(mesh->area_surface());
	gproshan_log_var(mesh->is_manifold());
	gproshan_log_var(mesh->bounds().size());
	gproshan_log_var(mesh->memory() / 1E6);
	gproshan_log_var(mesh->quality());
	gproshan_log_var(mesh->genus());
}


} // namespace gproshan

