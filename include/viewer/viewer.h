#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <map>
#include <cstring>

#include <glm/glm.hpp>

#include "camera.h"
#include "shader.h"
#include "frame.h"
#include "che_viewer.h"

#include "include_opengl.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


#ifdef GPROSHAN_EMBREE
	#include "rt_embree.h"
#endif // GPROSHAN_EMBREE

#ifdef GPROSHAN_OPTIX
	#include "rt_optix.h"
#endif // GPROSHAN_OPTIX


#define N_MESHES 12


// geometry processing and shape analysis framework
namespace gproshan {


class viewer;


struct vcorr_t
{
	index_t mesh_i;
	corr_t * corr;

	vcorr_t()
	{
		mesh_i = NIL;
		corr = nullptr;
	}

	void init(const size_t & n, const index_t & _mesh_i, const corr_t * _corr)
	{
		if(corr) delete [] corr;
		corr = new corr_t[n];

		mesh_i = _mesh_i;
		memcpy(corr, _corr, n * sizeof(corr_t));
	}

	operator corr_t *& ()
	{
		return corr;
	}

	bool is_loaded()
	{
		return mesh_i != NIL && corr != nullptr;
	}
};


class viewer
{
	protected:

		using function_t = void (*) (viewer *);
		
		struct process_t
		{
			std::string key;
			std::string name;
			function_t function;
			index_t sub_menu;
			
			process_t() = default;
			process_t(const std::string & k, const std::string & n, function_t f, const index_t & sm = NIL): key(k), name(n), function(f), sub_menu(sm) {};
		};

		static const int m_window_size[N_MESHES][2];


		GLFWwindow * window;
		int viewport_width;
		int viewport_height;

		shader shader_program;
		shader shader_normals;
		shader shader_gradient;
		shader shader_frame;
		
		camera cam;
		
		quaternion eye;
		quaternion center;
		quaternion up;
		
		quaternion light;

		glm::mat4 view_mat;
		glm::mat4 proj_mat;

		che_viewer meshes[N_MESHES];
		vcorr_t corr_mesh[N_MESHES];
		size_t n_meshes;
		index_t current; // current mesh

		index_t render_opt;

		frame * render_frame;

	#ifdef GPROSHAN_EMBREE
		rt::embree * rt_embree;
	#endif // GPROSHAN_EMBREE
	
	#ifdef GPROSHAN_OPTIX
		rt::optix * rt_optix;
	#endif // GPROSHAN_OPTIX

		bool action;

		bool render_wireframe;
		bool render_wireframe_fill;
		bool render_gradient_field;
		bool render_normal_field;
		bool render_border;
		bool render_lines;
		bool render_flat;
		float bgc;

		std::map<int, process_t> processes;

		che_viewer sphere;
		std::vector<vertex> sphere_translations;
		shader shader_sphere;
	
	public:

		std::vector<index_t> select_vertices;
		std::vector<vertex> other_vertices;
		std::vector<vertex> vectors;
		std::vector<std::string> sub_menus;

	public:

		viewer();
		virtual ~viewer();
		
		bool run();
		
		che_viewer & mesh();
		void add_process(const int & key, const process_t & process);
		void add_mesh(const std::vector<che *> & _meshes);
		
	private:

		// init
		void info_gl();
		void init_gl();
		void init_imgui();
		void init_menus();
		void init_glsl();
		void update_vbo();

		void render_gl();
		void render_embree();
		void render_optix();

		// callbacks
		static void idle();
		static void keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods);
		static void mouse_callback(GLFWwindow * window, int button, int action, int mods);
		static void cursor_callback(GLFWwindow * window, double x, double y);
		static void scroll_callback(GLFWwindow * window, double xoffset, double yoffset);

		// menu functions
		void menu_meshes(int value);
		void menu_process(int value);
		void menu_process(function_t pro);

		static void menu_help(viewer * view);
		static void menu_reset_mesh(viewer * view);
		static void menu_save_mesh(viewer * view);
		static void menu_zoom_in(viewer * view);
		static void menu_zoom_out(viewer * view);
		static void menu_bgc_inc(viewer * view);
		static void menu_bgc_dec(viewer * view);
		static void menu_bgc_white(viewer * view);
		static void menu_bgc_black(viewer * view);

		// render options
		static void invert_orientation(viewer * view);
		static void set_render_gl(viewer * view);
		static void set_render_embree(viewer * view);
		static void set_render_optix(viewer * view);
		static void set_render_wireframe(viewer * view);
		static void set_render_wireframe_fill(viewer * view);
		static void set_render_gradient_field(viewer * view);
		static void set_render_normal_field(viewer * view);
		static void set_render_border(viewer * view);
		static void set_render_lines(viewer * view);
		static void set_render_flat(viewer * view);
		
		static void raycasting(viewer * view);

		// draw routines
		void draw_scene();
		void draw_polygons();
		void draw_border();
		void draw_selected_vertices();
		
		void pick_vertex(int x, int y);
};


} // namespace gproshan

#endif // VIEWER_VIEWER_H

