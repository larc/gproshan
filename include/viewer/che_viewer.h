#ifndef CHE_VIEWER_H
#define CHE_VIEWER_H

#include "che.h"

#include "include_opengl.h"

#define COLOR 0.4

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

		size_t n_instances;
		size_t n_vertices; // current number of vertices
		bool invert_normals;
		vertex v_translate;

		vertex * normals;
		color_t * colors;

		GLuint vao;
		GLuint vbo[5];
	
	public:
		int vx, vy;					///< viewport positions.
		real_t factor;
		std::vector<index_t> selected;

	public:
		che_viewer(const size_t & n = 0);
		virtual ~che_viewer();
		che *& operator -> ();
		operator che *& ();
		void init(che * _mesh, const bool & normalize = true);
		void reload();
		void update();
		void update_vbo();
		void update_normals();
		void update_colors(const color_t *const c = nullptr);
		void update_instances_translations(const std::vector<vertex> & translations);
		void draw();

		color_t & color(const index_t & v);
		vertex & normal(const index_t & v);
		vertex *& normals_ptr();
		void translate(const vertex & p);
		void invert_orientation();

		void log_info();
};


} // namespace gproshan

#endif // CHE_VIEWER_H

