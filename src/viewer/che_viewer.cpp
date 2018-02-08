#include "che_viewer.h"
#include "viewer.h"

#include <GLES3/gl3.h>

#include <cassert>

che_viewer::che_viewer()
{
	mesh = NULL;
	_n_vertices = 0;

	normals = NULL;
	colors = NULL;
}

che_viewer::~che_viewer()
{
	if(!mesh) return;

	glDeleteBuffers(4, vbo);
	glDeleteVertexArrays(1, &vao);

	if(normals) delete [] normals;
	if(colors) delete [] colors;
	delete mesh;
}

che *& che_viewer::operator -> ()
{
	return mesh;
}

che_viewer::operator che *& ()
{
	return mesh;
}

void che_viewer::init(che * _mesh)
{
	_n_vertices = 0;
	mesh = _mesh;
	mesh->normalize();

	_invert_orientation = false;

	glGenVertexArrays(1, &vao);
	glGenBuffers(4, vbo);

	update();
}

void che_viewer::reload()
{
	_n_vertices = 0;
	mesh->reload();
	mesh->normalize();
	update();

	translate(v_translate);
	update();
}

void che_viewer::update()
{
	assert(mesh != NULL);

	if(_n_vertices != mesh->n_vertices())
	{
		if(normals) delete [] normals;
		if(colors) delete [] colors;

		_n_vertices = mesh->n_vertices();
		normals = new vertex[_n_vertices];
		colors = new color_t[_n_vertices];

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
	glBufferData(GL_ARRAY_BUFFER, _n_vertices * sizeof(vertex), &mesh->gt(0), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 1 NORMAL
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	glBufferData(GL_ARRAY_BUFFER, _n_vertices * sizeof(vertex), normals, GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 2 COLOR
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
	glBufferData(GL_ARRAY_BUFFER, _n_vertices * sizeof(vertex_t), colors, GL_STATIC_DRAW);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 1, GL_VERTEX_T, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// INDEXES
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[3]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->n_half_edges() * sizeof(index_t), &mesh->vt(0), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void che_viewer::update_normals()
{
	#pragma omp parallel for
	for(index_t v = 0; v < _n_vertices; v++)
	{
		normals[v] = mesh->normal(v);
		if(_invert_orientation) normals[v] = -normals[v];
	}
}

void che_viewer::update_colors(const color_t *const c)
{
	#pragma omp parallel for
	for(index_t v = 0; v < _n_vertices; v++)
		colors[v] = c ? c[v] : COLOR;
}

void che_viewer::draw()
{
	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[3]);
	glDrawElements(GL_TRIANGLES, mesh->n_half_edges(), GL_UNSIGNED_INT, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void che_viewer::draw_wireframe()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor4f(0., 0., 0., 0.5);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glBegin(GL_LINES);

	for(index_t e = 0; e < mesh->n_edges(); e++)
	{
		glVertex3v(&mesh->gt_vt(mesh->et(e)).x);
		glVertex3v(&mesh->gt_vt(next(mesh->et(e))).x);
	}

	glEnd();
	glPopAttrib();
}

void che_viewer::draw_normal_field()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor3f(.8, .8, 1.0);
	glLineWidth(2.0);

	glBegin(GL_LINES);
	for(index_t v = 0; v < _n_vertices; v++)
	{
		vertex n = factor * normals[v];
		vertex a = mesh->get_vertex(v);
		vertex b = a + n;

		glVertex3v(&a[0]);
		glVertex3v(&b[0]);
	}
	glEnd();

	glPopAttrib();
}

void che_viewer::draw_gradient_field()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor3f(.8, 1.0, .8);
	glLineWidth(1.2);

	double h = 0.3 * factor;

	for(index_t f = 0; f < mesh->n_faces(); f++)
	{
		vertex g = h * mesh->gradient_he(f * P, colors);
		vertex a = mesh->barycenter(f);
		vertex b = a + g;
		vertex n = mesh->normal_he(f * P);

		vertex v = b - a;
		vertex v90 = n * v;
		vertex p0 = b;
		vertex p1 = p0 - 0.25 * v - 0.15 * v90;
		vertex p2 = p0 - 0.25 * v + 0.15 * v90;

		glBegin(GL_LINES);
		glVertex3v(&a[0]);
		glVertex3v(&b[0]);
		glEnd();

		glBegin(GL_TRIANGLES);
		glVertex3v(&p0[0]);
		glVertex3v(&p1[0]);
		glVertex3v(&p2[0]);
		glEnd();
	}

	glPopAttrib();
}

void che_viewer::draw_mesh_info()
{
	if(!mesh) return;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, viewer::window_width(), 0, viewer::window_height());
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	char str[256];
	float color[4] = {1, .75, .25, 1};
	int h = 2, dh = 16;

	sprintf(str, "%9lu n_vertices", mesh->n_vertices());
	draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_9_BY_15);

	sprintf(str, "%9lu n_faces", mesh->n_faces());
	draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_9_BY_15);

	sprintf(str, "%9lu n_edges", mesh->n_edges());
	draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_9_BY_15);

	sprintf(str, "%9lu n_half_edges", mesh->n_half_edges());
	draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_9_BY_15);

	sprintf(str, "%9lu n_borders", mesh->n_borders());
	draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_9_BY_15);

	//sprintf(str, "%9.3lf quality", mesh->quality());
	//draw_str(str, 10, viewer::window_height() - (h += 18), color, GLUT_BITMAP_8_BY_13);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

const size_t & che_viewer::n_vertices() const
{
	return _n_vertices;
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
	_invert_orientation = !_invert_orientation;
}

void che_viewer::debug_info()
{
	if(!mesh) return;

	debug(mesh->n_vertices())
	debug(mesh->n_faces())
	debug(mesh->n_half_edges())
	debug(mesh->n_edges())
	debug(mesh->area_surface())
	debug(mesh->is_manifold())
	debug(mesh->n_borders())
	debug(mesh->memory() / 1E6)
	debug(mesh->quality())
	debug(mesh->genus())
}

