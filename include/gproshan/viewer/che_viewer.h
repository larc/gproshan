#ifndef CHE_VIEWER_H
#define CHE_VIEWER_H

#include "mesh/che.h"
#include "viewer/shader.h"
#include "raytracing/raytracing.h"

#include "viewer/include_opengl.h"


#ifdef GPROSHAN_FLOAT
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
		rt::raytracing * pick_vertex = nullptr;

		size_t n_instances = 0;
		bool normalize = false;
		vertex v_translate;

		GLuint vao;
		GLuint vbo[6];

	public:
		int vx, vy;							///< viewport positions.
		real_t factor;
		std::vector<index_t> selected;

		index_t idx_colormap	= 1;		// colormap index defined in shaders/colormap.glsl
		index_t point_size		= 1;
		bool point_normals		= true;
		bool render_pointcloud	= false;
		bool render_wireframe	= false;
		bool render_triangles	= false;
		bool render_gradients	= false;
		bool render_normals		= false;
		bool render_border		= false;
		bool render_lines		= false;
		bool render_flat		= false;

	public:
		che_viewer() = default;
		virtual ~che_viewer();
		che *& operator -> ();
		che *const & operator -> () const;
		operator che *& ();
		void init(che * mesh, const bool & normalize = true);
		void reload();
		void update();
		void update_vbo();
		void update_vbo_geometry();
		void update_vbo_normal(const vertex * vnormal = nullptr);
		void update_vbo_color(const che::rgb_t * vcolor = nullptr);
		void update_vbo_heatmap(const real_t * vheatmap = nullptr);
		void update_instances_translations(const std::vector<vertex> & translations);
		void draw(shader & program);
		void draw_point_cloud(shader & program);

		void translate(const vertex & p);
		void invert_orientation();
		void select(const real_t & x, const real_t & y, const glm::uvec2 & windows_size, const glm::mat4 & view_mat, const glm::mat4 & proj_mat);

		void log_info();
};


} // namespace gproshan

#endif // CHE_VIEWER_H

