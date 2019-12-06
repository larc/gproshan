#ifndef CHE_VIEWER_H
#define CHE_VIEWER_H

#include "che.h"

#include "include_opengl.h"

#define COLOR 0.5

#ifdef SINGLE_P
	#define glVertex3v(x) glVertex3fv(x)
	#define GL_VERTEX_T GL_FLOAT
#else
	#define glVertex3v(x) glVertex3dv(x)
	#define GL_VERTEX_T GL_DOUBLE
#endif


// geometry processing and shape analysis framework
namespace gproshan {


typedef real_t color_t;

class che_viewer
{
	protected:
		che * mesh;
		size_t _n_vertices; // current number of vertices
		bool _invert_orientation;
		real_t factor;
		vertex v_translate;

		vertex * normals;
		color_t * colors;

		GLuint vao;
		GLuint vbo[4];
	
	public:
		int vx, vy;					///< viewport positions.

	public:
		che_viewer();
		virtual ~che_viewer();
		che *& operator -> ();
		operator che *& ();
		void init(che * _mesh);
		void reload();
		void update();
		void update_vbo();
		void update_normals();
		void update_colors(const color_t *const c = nullptr);
		void draw();
		void draw_normal_field();
		void draw_gradient_field();
		void draw_mesh_info();

		const size_t & n_vertices() const;
		color_t & color(const index_t & v);
		vertex & normal(const index_t & v);
		vertex *& normals_ptr();
		void translate(const vertex & p);
		void invert_orientation();

		void log_info();
};


} // namespace gproshan

#endif // CHE_VIEWER_H

