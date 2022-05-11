#include "viewer/che_viewer.h"

#include <cassert>
#include <cmath>
#include <numeric>

#include <glm/gtc/matrix_transform.hpp>

#include "raytracing/rt_embree.h"


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
	factor = 1;

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

		factor = pmax.x - pmin.x;
		factor = std::max(factor, pmax.y - pmin.y);
		factor = std::max(factor, pmax.z - pmin.z);
		factor = 2.0 / factor;

		translate(- factor * (pmax + pmin) / 2);
		scale(factor);
	}

	factor *= mesh->mean_edge();

	render_pointcloud = mesh->is_pointcloud();

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

void che_viewer::translate(const vertex & p)
{
	model_mat = glm::translate(model_mat, glm_vec3(p));
}

void che_viewer::scale(const real_t & s)
{
	model_mat = glm::scale(model_mat, {s, s, s});
}

void che_viewer::select(const real_t & x, const real_t & y, const glm::uvec2 & windows_size, const glm::mat4 & proj_view_mat, const glm::vec3 & cam_pos)
{
	if(!pick_vertex) return;

	glm::mat4 inv_proj_view = glm::inverse(proj_view_mat);
	glm::vec2 screen = glm::vec2(float(x) / windows_size.x, float(windows_size.y - y) / windows_size.y);
	glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
	glm::vec4 q = inv_proj_view * view;
	glm::vec3 p = glm::vec3(q * (1.f / q.w));

	index_t v = pick_vertex->cast_ray(cam_pos, glm::normalize(p - cam_pos));
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

