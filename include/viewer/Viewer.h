#ifndef VIEWER_VIEWER_H
#define VIEWER_VIEWER_H

#include <GLES3/gl3.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <map>

#include "Camera.h"
#include "Shader.h"

#include "che_viewer.h"

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
		static void init(che * _mesh);
		
		static che_viewer mesh;
		static vector<index_t> select_vertices;
		static vector<vertex> other_vertices;
		static vector<vertex> vectors;
		static vector<string> sub_menus;

		static color_t & get_color(const index_t & i);
		static void add_process(const char & key, const string & name, function_t function);
	
	protected:
		// init
		static void debug_info();
		static void initGLUT( void );
		static void initGLSL( void );
		static void update_VBO( void );

		// GLUT callbacks
		static void display( void );
		static void idle( void );
		static void keyboard( unsigned char c, int x, int y );
		static void special( int i, int x, int y );
		static void mouse( int button, int state, int x, int y );
		static void motion( int x, int y );
		static void menu( int value );
		static void view( int value );
		static void menu_process( int value );
		
		// menu functions
		static void mProcess(function_t pro);
		static void mResetMesh( void );
		static void mWriteMesh();
		static void mExit( void );
		static void mWireframe( void );
		static void mZoomIn( void );
		static void mZoomOut( void );
		static void mGradientField( void );
		static void mNormalField( void );
		static void mBorder( void );
		static void mOrientation( void );
		static void mIsFlat( void );
		static void mLines( void );
		
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
		static void setGL( void );
		static void setLighting( void );
		static void setMeshMaterial( void );
		static void callDisplayList( void );
		static void updateDisplayList( void );
		static void drawScene( void );
		static void drawPolygons( void );
		static void drawWireframe( void );
		static void drawGradientField( void );
		static void drawNormalField( void );
		static void drawVertices( void );
		static void drawBorder( void );
		static void drawSelectedVertices( void );
		static void drawVectors( void );
		static void drawIsolatedVertices( void );
		static void pickVertex(int x, int y);	
		static void storeViewerState( void );
		static void restoreViewerState( void );
	
		static int windowSize[2];

		static Camera camera;
		// keeps track of view state
			
		static Shader shader;
		// shader used to determine appearance of surface
		
		static bool renderWireframe;
		static bool renderGradientField;
		static bool renderNormalField;
		static bool renderBorder;
		static bool invertOrientation;
		static bool is_flat;
		static bool lines;
		static float bgc;

		static map<unsigned char, process_t> processes;
	};
}

#endif

