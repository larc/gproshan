#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <GLES3/gl3.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <map>
#include <cstring>

#include "camera.h"
#include "shader.h"
#include "che_viewer.h"

#define N_MESHES 10

typedef void (*function_t) (void);

struct vcorr_t
{
	index_t mesh_i;
	corr_t * corr;

	vcorr_t()
	{
		mesh_i = NIL;
		corr = NULL;
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
		return mesh_i != NIL && corr != NULL;
	}
};

struct process_t
{
	index_t sub_menu;
	string name_function;
	function_t function;

	process_t(index_t _sub_menu = NIL, string _name_function = "", function_t _function = NULL)
	{
		sub_menu = _sub_menu;
		name_function = _name_function;
		function = _function;
	}
};

class viewer
{
	public:
		static void init(const vector<che *> & _meshes);
			
		static che_viewer meshes[N_MESHES];
		static vcorr_t corr_mesh[N_MESHES];
		static size_t n_meshes;
		static index_t current; // current mesh
		static vector<index_t> select_vertices;
		static vector<vertex> other_vertices;
		static vector<vertex> vectors;
		static vector<string> sub_menus;
		
		static che_viewer & mesh(); //get current che_viewer mesh
		static color_t & get_color(const index_t & i);
		static void add_process(const char & key, const string & name, function_t function);
		static void add_mesh(const vector<che *> & _meshes);
	
	protected:
		// init
		static void debug_info();
		static void initGLUT();
		static void init_menus();
		static void initGLSL();
		static void update_vbo();

		// GLUT callbacks
		static void display();
		static void idle();
		static void keyboard( unsigned char c, int x, int y );
		static void special( int i, int x, int y );
		static void mouse( int button, int state, int x, int y );
		static void motion( int x, int y );
		static void menu(int value);
		static void menu_process( int value );
		static void menu_meshes( int value );
		
		// menu functions
		static void mProcess(function_t pro);
		static void mResetMesh();
		static void mWriteMesh();
		static void mExit();
		static void mZoomIn();
		static void mZoomOut();
		static void invert_orientation();
		static void set_render_wireframe();
		static void set_render_gradient_field();
		static void set_render_normal_field();
		static void set_render_border();
		static void set_render_lines();
		static void set_render_corr();
		static void set_is_flat();
		
		// draw routines
		static void setGL();
		static void setLighting();
		static void setMeshMaterial();
		static void callDisplayList();
		static void updateDisplayList();
		static void drawScene();
		static void draw_corr();
		static void drawPolygons();
		static void drawWireframe();
		static void drawGradientField();
		static void drawNormalField();
		static void drawVertices();
		static void drawBorder();
		static void drawSelectedVertices();
		static void drawVectors();
		static void drawIsolatedVertices();
		static void pickVertex(int x, int y);	
		static void storeviewerState();
		static void restoreviewerState();
	
		static int windowSize[2];

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

		static map<unsigned char, process_t> processes;
};

#endif

