#ifndef VIEWER_H
#define VIEWER_H

#include <map>
#include <cstring>

#include <glm/glm.hpp>

#include "viewer/camera.h"
#include "viewer/shader.h"
#include "viewer/frame.h"
#include "viewer/che_viewer.h"

#include "raytracing/raytracing.h"

#include "viewer/include_opengl.h"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>


#ifdef GPROSHAN_FLOAT
	#define ImGui_InputReal ImGui::InputFloat
	#define ImGuiDataType_Real ImGuiDataType_Float
#else
	#define ImGui_InputReal ImGui::InputDouble
	#define ImGuiDataType_Real ImGuiDataType_Double
#endif // GPROSHAN_FLOAT

#define N_MESHES 12


// geometry processing and shape analysis framework
namespace gproshan {


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
		static const std::vector<std::string> colormap;


		GLFWwindow * window = nullptr;
		int window_width, window_height;
		int viewport_width, viewport_height;

		shader shader_triangles;
		shader shader_normals;
		shader shader_gradient;
		shader shader_pointcloud;

		camera cam;

		quaternion light;

		glm::mat4 view_mat;
		glm::mat4 proj_mat;

		che_viewer meshes[N_MESHES];
		size_t n_meshes	= 0;
		index_t idx_active_mesh = 0;

		enum render_type: index_t { R_GL, R_EMBREE, R_OPTIX };
		index_t render_opt = R_GL;

		frame * render_frame = nullptr;

		rt::raytracing * rt_embree = nullptr;
		rt::raytracing * rt_optix = nullptr;

		bool action = false;

		float bgc = 0;

		std::map<int, process_t> processes;

		che_viewer sphere;
		std::vector<vertex> sphere_translations;
		shader shader_sphere;

	public:
		std::vector<vertex> other_vertices;
		std::vector<vertex> vectors;
		std::vector<std::string> sub_menus;

		char status_message[1024] = {};

	public:
		viewer(const int & width = 1920, const int & height = 1080);
		virtual ~viewer();

		bool run();

		che_viewer & active_mesh();
		void add_process(const int & key, const process_t & process);
		void add_mesh(che * p_mesh);

	private:
		void info_gl();
		void init_gl();
		void init_imgui();
		void init_menus();
		void init_glsl();

		void render_gl();
		void render_rt(rt::raytracing * rt);

		static void framebuffer_size_callback(GLFWwindow * window, int width, int height);
		static void window_size_callback(GLFWwindow * window, int width, int height);
		static void keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods);
		static void mouse_callback(GLFWwindow * window, int button, int action, int mods);
		static void cursor_callback(GLFWwindow * window, double x, double y);
		static void scroll_callback(GLFWwindow * window, double xoffset, double yoffset);

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

		static bool setup_raytracing(viewer * view);
		static bool set_render_gl(viewer * view);
		static bool set_render_embree(viewer * view);
		static bool set_render_optix(viewer * view);

		static bool invert_orientation(viewer * view);
		static bool set_render_pointcloud(viewer * view);
		static bool set_render_wireframe(viewer * view);
		static bool set_render_triangles(viewer * view);
		static bool set_render_gradients(viewer * view);
		static bool set_render_normals(viewer * view);
		static bool set_render_border(viewer * view);
		static bool set_render_lines(viewer * view);
		static bool set_render_flat(viewer * view);

	#ifdef GPROSHAN_EMBREE
		static bool raycasting(viewer * view);
	#endif // GPROSHAN_EMBREE

		// draw routines
		void draw_selected_vertices(shader & program, const che_viewer & mesh);

		void select_border_vertices(che_viewer & mesh);
		void pick_vertex(const real_t & x, const real_t & y);
};


} // namespace gproshan

#endif // VIEWER_H

