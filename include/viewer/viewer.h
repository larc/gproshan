#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <map>
#include <cstring>

#include "camera.h"
#include "shader.h"
#include "che_viewer.h"

#include "include_opengl.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#ifdef GPROSHAN_EMBREE
	#include "embree.h"
#endif // GPROSHAN_EMBREE


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

		typedef void (*function_t) (viewer *);
		
		struct process_t
		{
			std::string key;
			std::string name;
			function_t function;
			index_t sub_menu;
		};

		static const int m_window_size[N_MESHES][2];


		GLFWwindow * window;

		shader shader_program;
		camera cam;
		
		int viewport_width;
		int viewport_height;
		
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

	#ifdef GPROSHAN_EMBREE
		embree * r_embree;
	#endif // GPROSHAN_EMBREE

		bool action;

		bool render_wireframe;
		bool render_gradient_field;
		bool render_normal_field;
		bool render_border;
		bool render_corr;
		bool render_lines;
		bool render_flat;
		float bgc;

		std::map<int, process_t> processes;
	
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
		static void set_render_gradient_field(viewer * view);
		static void set_render_normal_field(viewer * view);
		static void set_render_border(viewer * view);
		static void set_render_lines(viewer * view);
		static void set_render_corr(viewer * view);
		static void set_is_flat(viewer * view);
		
		static void raycasting(viewer * view);

		// draw routines
		void draw_scene();
		void draw_corr();
		void draw_polygons();
		void draw_wireframe();
		void draw_gradient_field();
		void draw_normal_field();
		void draw_vertices();
		void draw_border();
		void draw_selected_vertices();
		void draw_vectors();
		void draw_isolated_vertices();
		
		void pick_vertex(int x, int y);
};

void draw_str(const char * str, int x, int y, float color[4], void * font);


} // namespace gproshan

#endif // VIEWER_VIEWER_H

