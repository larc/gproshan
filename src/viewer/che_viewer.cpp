#include "che_viewer.h"

#include <cassert>
#include <cmath>
#include <numeric>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_viewer::che_viewer(const size_t & n): n_instances(n)
{
	mesh = nullptr;
	n_vertices = 0;

	normals = nullptr;
	colors = nullptr;
	
	vx = vy = 0;
}

che_viewer::~che_viewer()
{
	if(n_vertices) return;

	glDeleteBuffers(5, vbo);
	glDeleteVertexArrays(1, &vao);

	if(normals) delete [] normals;
	if(colors) delete [] colors;
}

che *& che_viewer::operator -> ()
{
	return mesh;
}

che_viewer::operator che *& ()
{
	return mesh;
}

void che_viewer::init(che * _mesh, const bool & normalize)
{
	n_vertices = 0;
	mesh = _mesh;
	
	if(normalize)
		mesh->normalize();

	invert_normals = false;

	glGenVertexArrays(1, &vao);
	glGenBuffers(5, vbo);

	update();
}

void che_viewer::reload()
{
	n_vertices = 0;
	mesh->reload();
	mesh->normalize();
	update();

	translate(v_translate);
	update();
}

void che_viewer::update()
{
	assert(mesh != nullptr);

	if(n_vertices != mesh->n_vertices())
	{
		if(normals) delete [] normals;
		if(colors) delete [] colors;

		n_vertices = mesh->n_vertices();
		normals = new vertex[n_vertices];
		colors = new color_t[n_vertices];

		update_normals();
		update_colors();
	}

	factor = mesh->mean_edge();

	update_vbo();
}

void che_viewer::update_vbo()
{
	glBindVertexArray(vao);

	// 0 VERTEX
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertex), &mesh->gt(0), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 1 NORMAL
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertex), normals, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 2 COLOR
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(real_t), colors, GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 1, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 3 INDEXES
	if(mesh->n_faces())
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[3]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->n_half_edges() * sizeof(index_t), &mesh->vt(0), GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
	else
	{
		index_t * indices = new index_t[n_vertices];
		iota(indices, indices + n_vertices, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[3]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_vertices * sizeof(index_t), indices, GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		delete [] indices;
	}
	
	glBindVertexArray(0);
}

void che_viewer::update_normals()
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		if(invert_normals)
			normals[v] = -mesh->normal(v);
		else
			normals[v] = mesh->normal(v);
}

void che_viewer::update_colors(const color_t *const c)
{
	if(!c)
	{
		#pragma omp parallel for
		for(index_t v = 0; v < n_vertices; v++)
			colors[v] = COLOR;

		return;
	}

	distance_t max_c = 0;
	
	#pragma omp parallel for reduction(max: max_c)
	for(index_t v = 0; v < n_vertices; v++)
		if(c[v] < INFINITY)
			max_c = max(c[v], max_c);

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; v++)
		colors[v] = c[v] / max_c;
}

void che_viewer::update_instances_translations(const vector<vertex> & translations)
{	
	n_instances = translations.size();
	if(!n_instances) return;

	glBindVertexArray(vao);

	// 4 INSTANCES (translations)
	glBindBuffer(GL_ARRAY_BUFFER, vbo[4]);
	glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(vertex), translations.data(), GL_STATIC_DRAW);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glVertexAttribDivisor(3, 1);

	glBindVertexArray(0);
}

void che_viewer::draw()
{
	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[3]);

	if(n_instances)
		glDrawElementsInstanced(GL_TRIANGLES, mesh->n_half_edges(), GL_UNSIGNED_INT, 0, n_instances);
	else if(mesh->n_faces())
		glDrawElements(GL_TRIANGLES, mesh->n_half_edges(), GL_UNSIGNED_INT, 0);
	else
		glDrawElements(GL_POINTS, mesh->n_vertices(), GL_UNSIGNED_INT, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

color_t & che_viewer::color(const index_t & v)
{
	return colors[v];
}

vertex & che_viewer::normal(const index_t & v)
{
	return normals[v];
}

vertex *& che_viewer::normals_ptr()
{
	return normals;
}

void che_viewer::translate(const vertex & p)
{
	v_translate = p;

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh->get_vertex(v) += v_translate;
}

void che_viewer::invert_orientation()
{
	invert_normals = !invert_normals;
}

void che_viewer::log_info()
{
	if(!mesh) return;

	gproshan_log_var(mesh->n_vertices());
	gproshan_log_var(mesh->n_faces());
	gproshan_log_var(mesh->n_half_edges());
	gproshan_log_var(mesh->n_edges());
	gproshan_log_var(mesh->area_surface());
	gproshan_log_var(mesh->is_manifold());
	gproshan_log_var(mesh->n_borders());
	gproshan_log_var(mesh->memory() / 1E6);
	gproshan_log_var(mesh->quality());
	gproshan_log_var(mesh->genus());
}


} // namespace gproshan

