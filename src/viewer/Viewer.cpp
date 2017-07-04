#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

using namespace std;

#include "Viewer.h"

void drawText(const char *text, int length, int x, int y)
{
	glColor3f(0, 1, 0);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION); // change the current matrix to PROJECTION
	glPushMatrix(); // push current state of MODELVIEW matrix to stack
	glLoadIdentity(); // reset PROJECTION matrix to identity matrix
	glOrtho(0, 800, 0, 600, -5, 5); // orthographic perspective
	glMatrixMode(GL_MODELVIEW); // change current matrix to MODELVIEW matrix again
	glLoadIdentity(); // reset it to identity matrix
	glPushMatrix(); // push current state of MODELVIEW matrix to stack

	glRasterPos2i(x, y); // raster position in 2D
	for(int i = 0; i < length; i++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, (int)text[i]); // generation of characters in our text with 9 by 15 GLU font
	
	glPopMatrix(); // reset
	glMatrixMode(GL_PROJECTION); // change current matrix mode to PROJECTION
	glPopMatrix(); // reset
	glMatrixMode(GL_MODELVIEW); // change current matrix mode to MODELVIEW
	
	glEnable(GL_LIGHTING);
}

namespace DDG
{
	// declare static member variables
	che_viewer Viewer::mesh;
	vector<index_t> Viewer::select_vertices;
	vector<vertex> Viewer::other_vertices;
	vector<vertex> Viewer::vectors;
	vector<string> Viewer::sub_menus;
	area_t Viewer::factor = 0;

	int Viewer::windowSize[2] = { 1366, 768 };
	Camera Viewer::camera;
	Shader Viewer::shader;
	bool Viewer::renderWireframe = false;
	bool Viewer::renderGradientField = false;
	bool Viewer::renderNormalField = false;
	bool Viewer::renderBorder = false;
	bool Viewer::invertOrientation = false;
	bool Viewer::is_flat = false;
	bool Viewer::lines = false;
	float Viewer::bgc = 0.;

	map<unsigned char, process_t> Viewer::processes;
	
	void Viewer::init(che * _mesh)
	{
		//restoreViewerState();
		initGLUT();
		mesh.init(_mesh);
		glutSetWindowTitle(mesh->filename().c_str());
		
		setGL();
		initGLSL();

		update_VBO();
		
		glutMainLoop();
	}

	void Viewer::debug_mesh_info()
	{
		if(!mesh) return;

		debug(mesh->n_vertices())
		debug(mesh->n_faces())
		debug(mesh->n_half_edges())
		debug(mesh->n_edges())
		debug(mesh->area_surface())
		debug(mesh->is_manifold())
		debug(mesh->n_borders())
		debug(mesh->memory() / 1E6)
	}

	void Viewer::debug_info()
	{
		const GLubyte *renderer = glGetString( GL_RENDERER );
		const GLubyte *vendor = glGetString( GL_VENDOR );
		const GLubyte *version = glGetString( GL_VERSION );
		const GLubyte *glslVersion = glGetString( GL_SHADING_LANGUAGE_VERSION );
		GLint major, minor;
		glGetIntegerv(GL_MAJOR_VERSION, &major);
		glGetIntegerv(GL_MINOR_VERSION, &minor);
		fprintf(stderr, "GL Vendor %s\n", vendor);
		fprintf(stderr, "GL Renderer %s\n", renderer);
		fprintf(stderr, "GL Version (string) %s\n", version);
		fprintf(stderr, "GL Version (integer) %d.%d\n", major, minor);
		fprintf(stderr, "GLSL Version %s\n", glslVersion);
	}

	void Viewer::initGLUT( void )
	{
		int argc = 0;
		vector< vector<char> > argv(1);
		
		// initialize window
		glutInitWindowSize( windowSize[0], windowSize[1] );
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
		glutInit( &argc, (char**)&argv );
		glutCreateWindow( "che_viewer" );
		//glutFullScreen();

		debug_info();

		// specify callbacks
		glutDisplayFunc( Viewer::display );
		glutIdleFunc( Viewer::idle );
		glutKeyboardFunc( Viewer::keyboard );
		glutSpecialFunc( Viewer::special );
		glutMouseFunc( Viewer::mouse );
		glutMotionFunc( Viewer::motion );

		glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION );
		
		// initialize menus
		int viewMenu = glutCreateMenu( Viewer::view );
		glutSetMenu( viewMenu );
		glutAddMenuEntry( "[f] Wireframe", menuWireframe );
		glutAddMenuEntry( "[g] GradientField", menuGradientField );
		glutAddMenuEntry( "[n] NormalField", menuNormalField );
		glutAddMenuEntry( "[b] Border", menuBorder );
		glutAddMenuEntry( "[i] Orientation", menuOrientation );
		glutAddMenuEntry( "[tab] Flat", menuIsFlat );
		glutAddMenuEntry( "[space] Lines", menuLines );

		int * sub_menu = new int[sub_menus.size()];

		for(index_t sm = 0; sm < sub_menus.size(); sm++)
			sub_menu[sm] = glutCreateMenu(Viewer::menu_process);

		for(auto mp: Viewer::processes)
		{
			stringstream ss;
			ss << "[" << mp.first << "] " << mp.second.name_function;
			
			glutSetMenu(sub_menu[mp.second.sub_menu]);
			glutAddMenuEntry( ss.str().c_str(), mp.first );
		}

		int mainMenu = glutCreateMenu( Viewer::menu );
		glutSetMenu( mainMenu );
		glutAddMenuEntry( "[r] Reset Mesh", menuResetMesh );
		glutAddMenuEntry( "[w] Write Mesh", menuWriteMesh );
		glutAddMenuEntry( "[<] Zoom In", menuZoomIn );
		glutAddMenuEntry( "[>] Zoom Out", menuZoomOut );
		glutAddMenuEntry( "[esc] Exit", menuExit );
		glutAddSubMenu( "View", viewMenu );
		
		for(index_t sm = 0; sm < sub_menus.size(); sm++)
			glutAddSubMenu(sub_menus[sm].c_str(), sub_menu[sm]);
		
		glutAttachMenu( GLUT_RIGHT_BUTTON );

		delete [] sub_menu;
	}
	
	void Viewer::initGLSL( void )
	{
		//shader.loadVertex( "shaders/vertex.glsl" );
		//shader.loadFragment( "shaders/fragment.glsl" );
		//shader.loadGeometry( "shaders/geometry.glsl" );
		//shader.loadGeometry( "shaders/new_geometry.glsl" );
		shader.loadVertex( "shaders/new_vertex.glsl" );
		shader.loadFragment( "shaders/new_fragment.glsl" );
	}

	color_t & Viewer::get_color(const index_t & i)
	{
		return mesh.color(i);
	}

	void Viewer::update_VBO( void )
	{
		mesh.update();
	}

	void Viewer::menu( int value )
	{
		switch( value )
		{
			case( menuResetMesh ):
				mResetMesh();
				break;
			case( menuWriteMesh ):
				mWriteMesh();
				break;
			case( menuExit ):
				mExit();
				break;
			default:
				break;
		}
	}
	
	void Viewer::menu_process(int value)
	{
		mProcess(processes[value].function);
	}
	
	void Viewer::add_process(const char & key, const string & name, function_t function)
	{
		if(processes.find(key) == processes.end())
		{
			processes[key] = process_t(sub_menus.size() - 1, name, function);
		}
		else cerr << "Repeat key: " << key << endl;  
	}

	void Viewer::view( int value )
	{
		switch( value )
		{
			case( menuWireframe ):
				mWireframe();
				break;
			case( menuZoomIn ):
				mZoomIn();
				break;
			case( menuZoomOut ):
				mZoomOut();
				break;
			case( menuGradientField ):
				mGradientField();
				break;
			case( menuNormalField ):
				mNormalField();
				break;
			case( menuBorder ):
				mBorder();
				break;
			case( menuOrientation ):
				mOrientation();
				break;
			case( menuIsFlat ):
				mIsFlat();
				break;
			case( menuLines ):
				mLines();
				break;
			default:
				break;
		}
	}
	
	void Viewer::keyboard( unsigned char c, int x, int y )
	{
		if(c >= '0' && c <= '9')
		{
			bgc = (c - '0') / 9.;
			glClearColor( bgc, bgc, bgc, 1. );
		}

		switch( c )
		{
			case 'x':
				mesh.update_colors();
				update_VBO();
				break;
			case 'c':
				select_vertices.clear();
				break;
			case 'f':
				mWireframe();
				break;
			case 'w':
				mWriteMesh();
				break;
			case 'r':
				mResetMesh();
				break;
			case 27:
				mExit();
				break;
			case 'g':
				mGradientField();
				break;
			case 'n':
				mNormalField();
				break;
			case 'b':
				mBorder();
				break;
			case 'i':
				mOrientation();
				break;
			case '\t':
				mIsFlat();
				break;
			case ' ':
				mLines();
				break;
			default:
				mProcess(processes[c].function);
				break;
		}
	}
	
	void Viewer::special( int i, int x, int y )
	{
		switch(i)
		{
			case GLUT_KEY_UP:
				camera.zoomIn();
				break;
			case GLUT_KEY_DOWN:
				camera.zoomOut();
				break;
			case 27:
				mExit();
				break;
			default:
				break;
		}
	}

	void Viewer::mouse( int button, int state, int x, int y )
	{
		if( ( glutGetModifiers() & GLUT_ACTIVE_SHIFT) and state == GLUT_UP )
			pickVertex(x, y);
		else if(button == 6) camera.zoomIn();
		else if(button == 5) camera.zoomOut();
		else camera.mouse( button, state, x, y );
	}
	
	void Viewer::motion( int x, int y )
	{
		camera.motion( x, y );
	}
	
	void Viewer::idle( void )
	{
		camera.idle();
		glutPostRedisplay();
	}
	
	void Viewer::storeViewerState( void )
	{
		ofstream out( ".viewer_state.txt" );
		
		out << camera.rLast[0] << endl;
		out << camera.rLast[1] << endl;
		out << camera.rLast[2] << endl;
		out << camera.rLast[3] << endl;
		
		GLint view[4];
		glGetIntegerv( GL_VIEWPORT, view );
		out << view[2] << endl;
		out << view[3] << endl;
	}
	
	void Viewer::restoreViewerState(void)
	{
		ifstream in( ".viewer_state.txt" );
		if( !in.is_open() ) return;
		
		in >> camera.rLast[0];
		in >> camera.rLast[1];
		in >> camera.rLast[2];
		in >> camera.rLast[3];
		in >> windowSize[0];
		in >> windowSize[1];
	}

	void Viewer::mProcess(function_t pro)
	{
		if(pro) pro();
		update_VBO();
	}
	
	void Viewer::mResetMesh()
	{
		mesh->reload();
		mesh->normalize();
		select_vertices.clear();
		other_vertices.clear();
		vectors.clear();
		factor = mesh->mean_edge();
		debug_mesh_info();
		update_VBO();
	}
	
	void Viewer::mWriteMesh()
	{
		string file = mesh->filename();
		index_t p = file.find_last_of('.');
		file = file.substr(0, p) + "_new.off";
		cout << __FUNCTION__ << " " << file << endl;
		mesh->write_file(file);
	}
	
	void Viewer::mExit( void )
	{
	//	storeViewerState();
		glutLeaveMainLoop();
	}
	
	void Viewer::mWireframe( void )
	{
		renderWireframe = !renderWireframe;
	}
	
	void Viewer::mZoomIn( void )
	{
		camera.zoomIn();
	}
	
	void Viewer::mZoomOut( void )
	{
		camera.zoomOut();
	}
	
	void Viewer::mGradientField( void )
	{
		renderGradientField = !renderGradientField;
	}

	void Viewer::mNormalField( void )
	{
		renderNormalField = !renderNormalField;
	}
		
	void Viewer::mBorder( void )
	{
		renderBorder = !renderBorder;
		if(!renderBorder) select_vertices.clear();
	}
	
	void Viewer::mOrientation( void )
	{
		invertOrientation = !invertOrientation;
		mesh.update_normals();
		update_VBO();
	}

	void Viewer::mIsFlat( void )
	{
		is_flat = !is_flat;
	}
	
	void Viewer::mLines( void )
	{
		lines = !lines;
	}
	
	void Viewer::display( void )
	{
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		shader.enable();

		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		GLint viewport[4];
		glGetIntegerv( GL_VIEWPORT, viewport );
		double aspect = (double) viewport[2] / (double) viewport[3];
		const double fovy = 50.;
		const double clipNear = .01;
		const double clipFar = 1000.;
		gluPerspective( fovy, aspect, clipNear, clipFar );
		
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		
	
		Quaternion eye = vertex( 0., 0., -2.5 * camera.zoom );
		Quaternion center = vertex( 0., 0., 0. );
		Quaternion up = vertex( 0., 1., 0. );		
		gluLookAt(	eye[1],		eye[2],		eye[3],
					center[1],	center[2],	center[3],
					up[1],		up[2],		up[3] );
		
		
		Quaternion r = camera.currentRotation();
		eye = r.conj() * eye * r;
		GLint uniformEye = glGetUniformLocation( shader, "eye" );
		glUniform3f( uniformEye, eye[1], eye[2], eye[3] );
		
		Quaternion light = vertex( -1., 1., -2. );
		light = r.conj() * light * r;
		GLint uniformLight = glGetUniformLocation( shader, "light" );
		glUniform3f( uniformLight, light[1], light[2], light[3] );
	
		camera.setView();
		
		GLint uniformIsFlat = glGetUniformLocation( shader, "is_flat" );
		glUniform1i(uniformIsFlat, is_flat);
		
		GLint uniformLines = glGetUniformLocation( shader, "lines" );
		glUniform1i(uniformLines, lines);

		GLfloat ModelViewMatrix[16]; 
		GLfloat ProjectionMatrix[16];

		glGetFloatv(GL_MODELVIEW_MATRIX, ModelViewMatrix);
		glGetFloatv(GL_PROJECTION_MATRIX, ProjectionMatrix); 

		GLint uniformModelViewMatrix = glGetUniformLocation(shader, "ModelViewMatrix");
		GLint uniformProjectionMatrix = glGetUniformLocation(shader, "ProjectionMatrix");
		
		glUniformMatrix4fv(uniformModelViewMatrix, 1, 0, ModelViewMatrix);
		glUniformMatrix4fv(uniformProjectionMatrix, 1, 0, ProjectionMatrix);

		//glPushAttrib( GL_ALL_ATTRIB_BITS );
		glEnable( GL_DEPTH_TEST );
		glEnable( GL_LIGHTING );
		
		setMeshMaterial();
		drawScene();
		
		//glPopAttrib();
		
		shader.disable();
		glutSwapBuffers();
	}
	
	void Viewer::setGL( void )
	{
		glClearColor( bgc, bgc, bgc, 1. );
		setLighting();
	}
	
	void Viewer::setLighting( void )
	{
		GLfloat position[4] = { 20., 30., 40., 0. };
		glLightfv( GL_LIGHT0, GL_POSITION, position );
		glEnable( GL_LIGHT0 );
		glEnable( GL_NORMALIZE );
	}
	
	void Viewer::setMeshMaterial( void )
	{
		GLfloat diffuse[4] = { .8, .5, .3, 1. };
		GLfloat specular[4] = { .3, .3, .3, 1. };
		GLfloat ambient[4] = { .2, .2, .5, 1. };
		
		glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,	diffuse );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,	ambient );
		glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 16.		);
	}
	
	void Viewer::drawScene( void )
	{
	//	glPushAttrib( GL_ALL_ATTRIB_BITS );

		glEnable( GL_POLYGON_OFFSET_FILL );
		glPolygonOffset( 1., 1. );
		drawPolygons();
		glDisable( GL_POLYGON_OFFSET_FILL );
		
		if( renderWireframe ) drawWireframe();
		if( renderGradientField ) drawGradientField();
		if( renderNormalField ) drawNormalField();
		if( renderBorder ) drawBorder();
		
		drawIsolatedVertices();
		drawVectors();
		drawSelectedVertices();
/*	
		char text[50];
		int n_text;
		n_text = sprintf(text, "%15s: %llu", "Vertices", mesh->n_vertices());
		drawText(text, n_text, 50, 100);
		n_text = sprintf(text, "%15s: %llu", "Faces", mesh->n_faces());
		drawText(text, n_text, 50, 88);
*/		
	//	glPopAttrib();
	}
	
	void Viewer::drawPolygons( void )
	{
		mesh.draw();
	}
	
	void Viewer::drawWireframe( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glDisable( GL_LIGHTING );
		glColor4f( 0., 0., 0., 0.5 );
		glEnable( GL_BLEND );
		glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
		
		glBegin( GL_LINES );

		for(index_t e = 0; e < mesh->n_edges(); e++)
		{
			glVertex3v( &mesh->gt_vt(mesh->et(e)).x );
			glVertex3v( &mesh->gt_vt(next(mesh->et(e))).x );
		}

		glEnd();
		glPopAttrib();
	}
	
	void Viewer::drawVectors( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glDisable( GL_LIGHTING );
		glColor4f( 1., 0., 0., 1 );
		glLineWidth( 3.0 );
		glEnable( GL_BLEND );
		glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
		
		glBegin( GL_LINES );
		
		index_t i = 0;
		for(vertex & v: vectors)
		{
			if(i % 8 == 0) glColor4f( 1., 0., 0., 1 );
			if(i % 8 == 2) glColor4f( 0., 1., 0., 1 );
			if(i % 8 == 4) glColor4f( 0., 0., 1., 1 );
			if(i % 8 == 6) glColor4f( 1., 1., 0., 1 );
			glVertex3v( &v.x );
			i++;
		}

		glEnd();
		glPopAttrib();
	}
	
	void Viewer::drawIsolatedVertices( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glEnable( GL_COLOR_MATERIAL );
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glColor3f( 0.5, 0., 0.5 );
		
		double h = 0.01 * camera.zoom; 
		for(const vertex & v: other_vertices)
		{
			glPushMatrix();
			glTranslated(v.x, v.y, v.z);
			glutSolidSphere(h, 10, 10);
			glPopMatrix();
		}

		glEnd();
		
		glPopAttrib();
		/*shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glPointSize( 5 );
		glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
		glEnable( GL_POINT_SMOOTH );
		glColor3f( 1., 0., 0. );
		
		glBegin( GL_POINTS );

		for(const vertex & v: other_vertices)
			glVertex3v(&v.x);

		glEnd();
		
		glPopAttrib();*/
	}
	
	void Viewer::drawVertices( void )
	{
		for(index_t v = 0; v < mesh->n_vertices(); v++)
		{
			glLoadName(v);
			glBegin(GL_POINTS);
			glVertex3v( &mesh->gt(v).x );
			glEnd();
		}
	}
	
	void Viewer::drawBorder( void )
	{
		select_vertices.clear();
		for(index_t b = 0; b < mesh->n_borders(); b++)
			for_border(he, mesh, mesh->bt(b))
				select_vertices.push_back(mesh->vt(he));
	}
	
	void Viewer::drawSelectedVertices( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glEnable( GL_COLOR_MATERIAL );
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glColor3f( 0., 0.5, 0.5 );
		
		double h = 0.02 * camera.zoom; 
		for(int v: select_vertices)
		{
			glPushMatrix();
			glTranslated(mesh->gt(v).x, mesh->gt(v).y, mesh->gt(v).z);
			glutSolidSphere(h, 10, 10);
			glPopMatrix();
		}

		glEnd();
		
		glPopAttrib();
	}

	void Viewer::drawNormalField( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );

		glDisable( GL_LIGHTING );
		glColor3f( .8, .8, 1. );
		glLineWidth( 1.0 );

		for(index_t v = 0; v < mesh->n_vertices(); v++)
		{
			vertex n = factor * mesh.normal(v);
			vertex a = mesh->get_vertex(v);
			vertex b = a + n;
			
			glBegin( GL_LINES );
			glVertex3v( &a[0] );
			glVertex3v( &b[0] );
			glEnd();
		}

		glPopAttrib();
	}

	void Viewer::drawGradientField( void )
	{
		shader.disable();
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glDisable( GL_LIGHTING );
		glColor3f( .8, 1., .8 );
		glLineWidth( 1.0 );

		double h = 0.3 * factor;

		for(index_t f = 0; f < mesh->n_faces(); f++)
		{
			vertex g = h * mesh->gradient_he(f * P, &mesh.color(0));
			vertex a = mesh->barycenter(f);
			vertex b = a + g;
			vertex n = mesh->normal_he(f * P);

			vertex v = b - a;
			vertex v90 = n * v;
			vertex p0 = b;
			vertex p1 = p0 - 0.25 * v - 0.15 * v90;
			vertex p2 = p0 - 0.25 * v + 0.15 * v90;

			glBegin( GL_LINES );
			glVertex3v( &a[0] );
			glVertex3v( &b[0] );
			glEnd();

			glBegin(GL_TRIANGLES);
			glVertex3v( &p0[0] );
			glVertex3v( &p1[0] );
			glVertex3v( &p2[0] );
			glEnd();
		}
		
		glPopAttrib();
	}

	void Viewer::pickVertex(int x, int y)
	{
		int width = glutGet(GLUT_WINDOW_WIDTH);
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		if( x < 0 || x >= width || y < 0 || y >= height ) return;
		
		int bufSize = mesh->n_vertices();
		GLuint * buf = new GLuint[bufSize];
		glSelectBuffer(bufSize, buf);
		
		GLint viewport[4];
		GLdouble projection[16];
		glGetIntegerv( GL_VIEWPORT, viewport );
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		
		glRenderMode(GL_SELECT);
		glInitNames();
		glPushName(0);
		
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluPickMatrix(x, viewport[3] - y, 10, 10, viewport);
		glMultMatrixd(projection);
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		drawVertices();
		glPopMatrix();
		
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		
		glMatrixMode(GL_MODELVIEW);
		long hits = glRenderMode(GL_RENDER);

		int index = -1;
		double min_z = 1.0e100;
		for( long i = 0; i < hits; ++i )
		{
			double distance = buf[4*i + 1];
			if( distance < min_z )
			{
				index = buf[4*i + 3];
				min_z = distance;
			}
		}
		delete[] buf;

		if (index >= 0 && index < mesh->n_vertices())
		{
			debug(index)
			debug(mesh.color(index))
			select_vertices.push_back(index);
		}
	}
}

