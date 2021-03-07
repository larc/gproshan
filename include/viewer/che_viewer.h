#ifndef CHE_VIEWER_H
#define CHE_VIEWER_H

#include "mesh/che.h"
#include "viewer/shader.h"

#include "viewer/include_opengl.h"


#ifdef SINGLE_P
	#define glVertex3v(x) glVertex3fv(x)
	#define GL_REAL GL_FLOAT
#else
	#define glVertex3v(x) glVertex3dv(x)
	#define GL_REAL GL_DOUBLE
#endif


// geometry processing and shape analysis framework
namespace gproshan {


class che_viewer
{
	protected:
		che * mesh = nullptr;

		size_t n_instances = 0;
		bool normalize = false;
		vertex v_translate;

		GLuint vao;
		GLuint vbo[6];

	public:
		int vx, vy;					///< viewport positions.
		real_t factor;
		std::vector<index_t> selected;

	public:
		che_viewer() = default;
		virtual ~che_viewer();
		che *& operator -> ();
		operator che *& ();
		void init(che * mesh, const bool & normalize = true);
		void reload();
		void update();
		void update_vbo();
		void update_instances_translations(const std::vector<vertex> & translations);
		void draw(shader & program);
		void draw_point_cloud(shader & program);

		void translate(const vertex & p);
		void invert_orientation();

		void log_info();
};


} // namespace gproshan

#endif // CHE_VIEWER_H

