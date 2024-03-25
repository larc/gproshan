#ifndef VIEWER_H
#define VIEWER_H

#include <gproshan/mesh/che_sphere.h>
#include <gproshan/viewer/camera.h>
#include <gproshan/viewer/shader.h>
#include <gproshan/viewer/frame.h>
#include <gproshan/viewer/che_viewer.h>
#include <gproshan/viewer/include_opengl.h>
#include <gproshan/scenes/scene.h>
#include <gproshan/raytracing/render_params.h>

#include <cstring>
#include <functional>
#include <map>

#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>


// geometry processing and shape analysis framework
namespace gproshan {


class viewer
{
	protected:

		using function_t = bool (*) (viewer *);

		struct process_t
		{
			const char * key = nullptr;
			const char * name = nullptr;
			function_t function = nullptr;
			index_t id_menu = NIL;
			bool selected = false;
		};

		static const std::vector<ivec2> m_window_split;
		static const std::vector<std::string> colormap;
		static const size_t max_meshes;
		static const size_t max_nframes = 1000;

		static che_sphere sphere_data;


		bool apply_all_meshes = false;


		GLFWwindow * window = nullptr;
		rt::render_params render_params;
		unsigned int & window_width = render_params.window_size.x();
		unsigned int & window_height = render_params.window_size.y();
		unsigned int & viewport_width = render_params.viewport_size.x();
		unsigned int & viewport_height = render_params.viewport_size.y();
		mat4 proj_mat;
		mat4 proj_view_mat;

		bool hide_imgui = false;

		shader shader_triangles;
		shader shader_normals;
		shader shader_gradient;
		shader shader_pointcloud;
		shader shader_depth;
		scene::material mat;

		camera cam;
		quaternion cam_light;

		double render_time = 0;
		double frametime[max_nframes] = {};
		index_t nframes = 0;

		std::vector<che_viewer *> meshes;
		std::vector<che_viewer *> removed_meshes;
		index_t idx_selected_mesh = 0;

		frame * frames = nullptr;

		float bgc = 0;

		std::vector<std::string> menus;
		std::vector<std::vector<int> > menu_processes;
		std::unordered_map<int, process_t> processes;

		che_viewer * sphere = nullptr;
		shader shader_sphere;
		std::vector<vertex> sphere_points;

		std::vector<vertex> vectors;

		char status_message[1024] = {};

	public:
		viewer(const char * title = "gproshan", const int width = 1600, const int height = 900);
		virtual ~viewer();

		che_viewer & selected_mesh();
		void add_menu(const std::string & str, const std::vector<int> & vprocesses);
		int add_process(const char * name, const function_t & f, const int key = -1);
		bool add_mesh(che * p_mesh, const bool reset_normals = true);
		bool remove_mesh(const index_t idx);
		bool pop_mesh();
		void update_viewport_meshes();
		void update_status_message(const char * format, ...);

	protected:
		virtual bool run();

		void info_gl();
		void init_gl(const char * title);
		void init_imgui();
		void init_menus();
		void init_glsl();

		void imgui();

		void render_gl();
		void render_rt(che_viewer & mesh, frame & rt_frame);

		void save_history(const std::string & file);
		void save_frametime(const std::string & file);

		static void framebuffer_size_callback(GLFWwindow * window, int width, int height);
		static void window_size_callback(GLFWwindow * window, int width, int height);
		static void keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods);
		static void mouse_callback(GLFWwindow * window, int button, int action, int mods);
		static void cursor_callback(GLFWwindow * window, double x, double y);
		static void scroll_callback(GLFWwindow * window, double xoffset, double yoffset);

		static bool m_help(viewer * view);
		static bool m_close(viewer * view);
		static bool m_maximize(viewer * view);
		static bool m_hide_show_imgui(viewer * view);

		static bool m_save_load_view(viewer * view);
		static bool m_reset_mesh(viewer * view);
		static bool m_save_mesh(viewer * view);
		static bool m_remove_mesh(viewer * view);
		static bool m_pop_mesh(viewer * view);
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

		void pick_vertex(const uvec2 & pos);
		void check_apply_all_meshes(const std::function<void(che_viewer &)> & fun);
};


} // namespace gproshan

#endif // VIEWER_H

