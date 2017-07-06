#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <GLES3/gl3.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <map>

#include "Camera.h"
#include "Shader.h"

#include "che_viewer.h"

#define N_MESHES 10

typedef void (*function_t) (void);

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

namespace DDG
{
	class Viewer
	{
	public:
		static void init(const vector<che *> & _meshes);
			
		static che_viewer meshes[N_MESHES];
		static corr_t * corr_mesh[N_MESHES];
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
		static void update_VBO();

		// GLUT callbacks
		static void display();
		static void idle();
		static void keyboard( unsigned char c, int x, int y );
		static void special( int i, int x, int y );
		static void mouse( int button, int state, int x, int y );
		static void motion( int x, int y );
		static void menu( int value );
		static void view( int value );
		static void menu_process( int value );
		static void menu_meshes( int value );
		
		// menu functions
		static void mProcess(function_t pro);
		static void mResetMesh();
		static void mWriteMesh();
		static void mExit();
		static void mWireframe();
		static void mZoomIn();
		static void mZoomOut();
		static void mGradientField();
		static void mNormalField();
		static void mBorder();
		static void mOrientation();
		static void mIsFlat();
		static void mLines();
		
		// unique identifiers for menus
		enum
		{
			menuProcess = 10000,
			menuResetMesh = 10001,
			menuWriteMesh = 10002,
			menuExit = 10003,
			menuWireframe = 10004,
			menuZoomIn = 10005,
			menuZoomOut = 10006,
			menuGradientField = 10007,
			menuNormalField = 1008,
			menuBorder = 1009,
			menuOrientation = 1010,
			menuIsFlat = 1011,
			menuLines = 1012
		};
		
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
		static void storeViewerState();
		static void restoreViewerState();
	
		static int windowSize[2];

		static Camera camera;
		// keeps track of view state
			
		static Shader shader;
		// shader used to determine appearance of surface
		
		static bool renderWireframe;
		static bool renderGradientField;
		static bool renderNormalField;
		static bool renderBorder;
		static bool is_flat;
		static bool lines;
		static float bgc;

		static map<unsigned char, process_t> processes;
	};
}

#endif

