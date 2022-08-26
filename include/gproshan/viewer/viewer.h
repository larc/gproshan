#ifndef VIEWER_H
#define VIEWER_H

#include <cstring>
#include <functional>
#include <map>

#include <gproshan/viewer/camera.h>
#include <gproshan/viewer/shader.h>
#include <gproshan/viewer/frame.h>
#include <gproshan/viewer/che_viewer.h>
#include <gproshan/viewer/include_opengl.h>
#include <gproshan/raytracing/render_params.h>

#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>


#ifdef GPROSHAN_FLOAT
	#define ImGui_InputReal ImGui::InputFloat
	#define ImGuiDataType_Real ImGuiDataType_Float
#else
	#define ImGui_InputReal ImGui::InputDouble
	#define ImGuiDataType_Real ImGuiDataType_Double
#endif // GPROSHAN_FLOAT


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

		static const std::vector<ivec2> m_window_split;
		static const size_t max_n_meshes;
		static const std::vector<std::string> colormap;

		bool apply_all_meshes = false;


		GLFWwindow * window = nullptr;
		rt::render_params render_params;
		int & window_width = render_params.window_width;
		int & window_height = render_params.window_height;
		int & viewport_width = render_params.viewport_width;
		int & viewport_height = render_params.viewport_height;
		mat4 proj_view_mat;

		bool hide_imgui = false;

		shader shader_triangles;
		shader shader_normals;
		shader shader_gradient;
		shader shader_pointcloud;

		camera cam;
		quaternion cam_light;

		double render_time = 0;

		che_viewer * meshes = nullptr;
		size_t n_meshes	= 0;
		index_t idx_active_mesh = 0;

		frame * frames = nullptr;

		float bgc = 0;

		std::map<int, process_t> processes;

		che_viewer sphere;
		shader shader_sphere;
		std::vector<vertex> sphere_points;

		std::vector<vertex> vectors;
		std::vector<std::string> sub_menus;

		char status_message[1024] = {};

	public:
		viewer(const int & width = 1920, const int & height = 1080);
		virtual ~viewer();

		che_viewer & active_mesh();
		void add_process(const int & key, const std::string & skey, const std::string & name, const function_t & f);
		bool add_mesh(che * p_mesh);

	protected:
		virtual bool run();

		void info_gl();
		void init_gl();
		void init_imgui();
		void init_menus();
		void init_glsl();

		void imgui();

		void render_gl();
		void render_rt(che_viewer & mesh, frame & rt_frame);

		static void framebuffer_size_callback(GLFWwindow * window, int width, int height);
		static void window_size_callback(GLFWwindow * window, int width, int height);
		static void keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods);
		static void mouse_callback(GLFWwindow * window, int button, int action, int mods);
		static void cursor_callback(GLFWwindow * window, double x, double y);
		static void scroll_callback(GLFWwindow * window, double xoffset, double yoffset);

		static bool m_help(viewer * view);
		static bool m_close(viewer * view);
		static bool m_hide_show_imgui(viewer * view);

		static bool m_save_load_view(viewer * view);
		static bool m_reset_mesh(viewer * view);
		static bool m_save_mesh(viewer * view);
		static bool m_normalize_mesh(viewer * view);
		static bool m_zoom_in(viewer * view);
		static bool m_zoom_out(viewer * view);
		static bool m_bgc_inc(viewer * view);
		static bool m_bgc_dec(viewer * view);
		static bool m_bgc_white(viewer * view);
		static bool m_bgc_black(viewer * view);

		static bool m_setup_raytracing(viewer * view);
		static bool m_render_gl(viewer * view);
		static bool m_render_embree(viewer * view);
		static bool m_render_optix(viewer * view);

		static bool m_invert_normals(viewer * view);
		static bool m_select_border_vertices(viewer * view);
		static bool m_clean_selected_vertices(viewer * view);
		static bool m_render_pointcloud(viewer * view);
		static bool m_render_wireframe(viewer * view);
		static bool m_render_triangles(viewer * view);
		static bool m_render_gradients(viewer * view);
		static bool m_render_normals(viewer * view);
		static bool m_render_lines(viewer * view);
		static bool m_render_flat(viewer * view);

		static bool m_raycasting(viewer * view);

		void pick_vertex(const int & x, const int & y);
		void check_apply_all_meshes(const std::function<void(che_viewer &)> & fun);
};


} // namespace gproshan

#endif // VIEWER_H

