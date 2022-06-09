#ifndef CHE_VIEWER_H
#define CHE_VIEWER_H

#include <gproshan/mesh/che.h>
#include <gproshan/viewer/shader.h>
#include <gproshan/raytracing/raytracing.h>

#include <gproshan/viewer/include_opengl.h>


#ifdef GPROSHAN_FLOAT
	#define glVertex3v(x) glVertex3fv(x)
	#define GL_REAL GL_FLOAT
#else
	#define glVertex3v(x) glVertex3dv(x)
	#define GL_REAL GL_DOUBLE
#endif


// geometry processing and shape analysis framework
namespace gproshan {


enum render_type: index_t { R_GL, R_EMBREE, R_OPTIX };

class che_viewer
{
	protected:
		che * mesh = nullptr;

		size_t n_instances = 0;
		bool center_mesh = false;
		vertex v_translate;

		GLuint vao;
		GLuint vbo[6];

	public:
		int vx, vy;							///< viewport positions.
		std::vector<index_t> selected;
		std::vector<vertex> selected_xyz;
		rt::raytracing * rt_embree	= nullptr;
		rt::raytracing * rt_optix	= nullptr;

		mat4 model_mat = mat4::identity();

		index_t idx_colormap	= 1;		// colormap index defined in shaders/colormap.glsl
		index_t point_size		= 1;
		index_t render_opt		= R_GL;
		bool point_normals		= true;
		bool render_pointcloud	= false;
		bool render_wireframe	= false;
		bool render_triangles	= false;
		bool render_gradients	= false;
		bool render_normals		= false;
		bool render_lines		= false;
		bool render_flat		= false;

	public:
		che_viewer() = default;
		virtual ~che_viewer();

		che *& operator -> ();
		che *const & operator -> () const;
		operator che *& ();

		void init(che * m, const bool & center = true);
		void update();
		void update_vbo();
		void update_vbo_geometry();
		void update_vbo_normal(const vertex * vnormal = nullptr);
		void update_vbo_color(const che::rgb_t * vcolor = nullptr);
		void update_vbo_heatmap(const real_t * vheatmap = nullptr);
		void update_instances_positions(const std::vector<vertex> & translations);

		void draw(shader & program);
		void draw_point_cloud(shader & program);
		void draw_selected_vertices(che_viewer & sphere, shader & program);

		void translate(const vertex & p);
		void scale(const real_t & s);

		void select(const index_t & x, const index_t & y, const uvec2 & windows_size, const mat4 & proj_view_mat, const vertex & cam_pos);

		void log_info();
};


} // namespace gproshan

#endif // CHE_VIEWER_H

