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


class viewer
{
	protected:

		using function_t = bool (*) (viewer *);
		
		struct process_t
		{
			std::string key;
			std::string name;
			function_t function;
			index_t sub_menu;
			bool selected = false;
			
			process_t() = default;
			process_t(const std::string & k, const std::string & n, function_t f, const index_t & sm = NIL): key(k), name(n), function(f), sub_menu(sm) {};
		};

		static const int m_window_size[N_MESHES + 1][2];


		GLFWwindow * window = nullptr;
		int viewport_width;
		int viewport_height;

		shader shader_program;
		shader shader_normals;
		shader shader_gradient;
		shader shader_pointcloud;
		
		camera cam;
		
		quaternion eye;
		quaternion center;
		quaternion up;
		
		quaternion light;

		glm::mat4 view_mat;
		glm::mat4 proj_mat;

		che_viewer meshes[N_MESHES];
		size_t n_meshes	= 0;
		index_t idx_active_mesh = 0; // idx_active_mesh mesh

		index_t render_opt = 0;
		
		frame * render_frame = nullptr;

	#ifdef GPROSHAN_EMBREE
		rt::embree * rt_embree = nullptr;
	#endif // GPROSHAN_EMBREE
	
	#ifdef GPROSHAN_OPTIX
		rt::optix * rt_optix = nullptr;
	#endif // GPROSHAN_OPTIX

		bool action = false;

		bool render_wireframe = false;
		bool render_wireframe_fill = false;
		bool render_gradient_field = false;
		bool render_normal_field = false;
		bool render_border = false;
		bool render_lines = false;
		bool render_flat = false;
		float bgc = 0;

		std::map<int, process_t> processes;

		che_viewer sphere;
		std::vector<vertex> sphere_translations;
		shader shader_sphere;
	
	public:

		std::vector<vertex> other_vertices;
		std::vector<vertex> vectors;
		std::vector<std::string> sub_menus;

	public:

		viewer();
		virtual ~viewer();
		
		bool run();
		
		che_viewer & active_mesh();
		void add_process(const int & key, const process_t & process);
		void add_mesh(che * p_mesh);
		
	private:

		// init
		void info_gl();
		void init_gl();
		void init_imgui();
		void init_menus();
		void init_glsl();

		void render_gl();
		void render_embree();
		void render_optix();

		// callbacks
		static void keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods);
		static void mouse_callback(GLFWwindow * window, int button, int action, int mods);
		static void cursor_callback(GLFWwindow * window, double x, double y);
		static void scroll_callback(GLFWwindow * window, double xoffset, double yoffset);

		// menu functions
		static bool menu_help(viewer * view);
		static bool menu_save_load_view(viewer * view);
		static bool menu_reset_mesh(viewer * view);
		static bool menu_save_mesh(viewer * view);
		static bool menu_zoom_in(viewer * view);
		static bool menu_zoom_out(viewer * view);
		static bool menu_bgc_inc(viewer * view);
		static bool menu_bgc_dec(viewer * view);
		static bool menu_bgc_white(viewer * view);
		static bool menu_bgc_black(viewer * view);

		// render options
		static bool invert_orientation(viewer * view);
		static bool set_render_gl(viewer * view);
		static bool set_render_embree(viewer * view);
		static bool set_render_optix(viewer * view);
		static bool set_render_wireframe(viewer * view);
		static bool set_render_wireframe_fill(viewer * view);
		static bool set_render_gradient_field(viewer * view);
		static bool set_render_normal_field(viewer * view);
		static bool set_render_border(viewer * view);
		static bool set_render_lines(viewer * view);
		static bool set_render_flat(viewer * view);
		
		static bool raycasting(viewer * view);

		// draw routines
		void draw_meshes(shader & program);
		void draw_selected_vertices(shader & program);
		
		void select_border_vertices();
		void pick_vertex(int x, int y);
};


} // namespace gproshan

#endif // VIEWER_VIEWER_H

