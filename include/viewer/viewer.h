#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <map>
#include <cstring>

#include "camera.h"
#include "shader.h"
#include "che_viewer.h"

#include "include_opengl.h"

#define N_MESHES 12


// geometry processing and shape analysis framework
namespace gproshan {


typedef void (*function_t) (void);

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

struct process_t
{
	index_t sub_menu;
	std::string name_function;
	function_t function;
};

class viewer
{
	public:
		static void init(const std::vector<che *> & _meshes);

		static che_viewer meshes[N_MESHES];
		static vcorr_t corr_mesh[N_MESHES];
		static size_t n_meshes;
		static index_t current; // current mesh

		static std::vector<index_t> select_vertices;
		static std::vector<vertex> other_vertices;
		static std::vector<vertex> vectors;
		static std::vector<std::string> sub_menus;

		static char * share;

		static const int & window_width();
		static const int & window_height();

		static che_viewer & mesh(); //get current che_viewer mesh
		static color_t & vcolor(const index_t & i);
		static void add_process(const char & key, const std::string & name, function_t function);
		static void add_mesh(const std::vector<che *> & _meshes);

	protected:
		// init
		static void debug_info();
		static void init_glut();
		static void init_menus();
		static void init_glsl();
		static void update_vbo();

		// GLUT callbacks
		static void display();
		static void idle();
		static void keyboard(unsigned char c, int x, int y);
		static void special(int i, int x, int y);
		static void mouse(int button, int state, int x, int y);
		static void motion(int x, int y);

		// menu functions
		static void menu(int value);
		static void menu_process(int value);
		static void menu_meshes(int value);
		static void menu_process(function_t pro);
		static void menu_reset_mesh();
		static void menu_save_mesh();
		static void menu_exit();
		static void menu_zoom_in();
		static void menu_zoom_out();

		// render options
		static void invert_orientation();
		static void set_render_wireframe();
		static void set_render_gradient_field();
		static void set_render_normal_field();
		static void set_render_border();
		static void set_render_lines();
		static void set_render_corr();
		static void set_is_flat();

		// draw routines
		static void set_gl();
		static void set_lighting();
		static void set_mesh_materia();
		static void draw_scene();
		static void draw_corr();
		static void draw_polygons();
		static void draw_wireframe();
		static void draw_gradient_field();
		static void draw_normal_field();
		static void draw_vertices();
		static void draw_border();
		static void draw_selected_vertices();
		static void draw_vectors();
		static void draw_isolated_vertices();
		static void pick_vertex(int x, int y);

		static int window_size[2];
		static double ww, wh;
		static int m_window_size[N_MESHES][2];

		static camera cam;
		// keeps track of view state

		static shader shader_program;
		// shader used to determine appearance of surface

		static bool render_wireframe;
		static bool render_gradient_field;
		static bool render_normal_field;
		static bool render_border;
		static bool render_corr;
		static bool render_lines;
		static bool is_flat;
		static float bgc;

		static std::map<unsigned char, process_t> processes;
};

void draw_str(const char * str, int x, int y, float color[4], void * font);


} // namespace gproshan

#endif // VIEWER_VIEWER_H

