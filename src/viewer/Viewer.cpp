#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
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
	che_viewer Viewer::meshes[N_MESHES];
	vcorr_t Viewer::corr_mesh[N_MESHES]; // zero initialization 
	size_t Viewer::n_meshes = 0;
	index_t Viewer::current = 0;
	vector<index_t> Viewer::select_vertices;
	vector<vertex> Viewer::other_vertices;
	vector<vertex> Viewer::vectors;
	vector<string> Viewer::sub_menus;

	int Viewer::windowSize[2] = { 1366, 768 };
	Camera Viewer::camera;
	Shader Viewer::shader;
	bool Viewer::renderWireframe = false;
	bool Viewer::renderGradientField = false;
	bool Viewer::renderNormalField = false;
	bool Viewer::renderBorder = false;
	bool Viewer::render_corr = false;
	bool Viewer::is_flat = false;
	bool Viewer::lines = false;
	float Viewer::bgc = 0.;

	map<unsigned char, process_t> Viewer::processes;
	
	che_viewer & Viewer::mesh()
	{
		assert(n_meshes > 0);
		return meshes[current];
	}

	void Viewer::init(const vector<che *> & _meshes)
	{
		//restoreViewerState();
		initGLUT();
		add_mesh(_meshes);	

		glutSetWindowTitle(mesh()->filename().c_str());
		init_menus();
		
		debug_info();
		mesh().debug_info();

		setGL();
		initGLSL();

		glutMainLoop();
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

	void Viewer::initGLUT()
	{
		int argc = 0;
		vector< vector<char> > argv(1);
		
		// initialize window
		glutInitWindowSize( windowSize[0], windowSize[1] );
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
		glutInit( &argc, (char**)&argv );
		glutCreateWindow( "che_viewer" );
		//glutFullScreen();

		// specify callbacks
		glutDisplayFunc( Viewer::display );
		glutIdleFunc( Viewer::idle );
		glutKeyboardFunc( Viewer::keyboard );
		glutSpecialFunc( Viewer::special );
		glutMouseFunc( Viewer::mouse );
		glutMotionFunc( Viewer::motion );

		glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION );
	}	

	void Viewer::init_menus()
	{
		// set current mesh menu
		int mesh_menu = glutCreateMenu( Viewer::menu_meshes );
		glutSetMenu(mesh_menu);
		for(index_t i = 0; i < n_meshes; i++)
			glutAddMenuEntry(meshes[i]->filename().c_str(), i);
		// init viewer menu
		sub_menus.push_back("Viewer");
		add_process('f', "Wireframe", mWireframe);
		add_process('g', "Gradient Field", mGradientField);
		add_process('n', "Normal Field", mNormalField);
		add_process('b', "Show Border", mBorder);
		add_process('\t', "Flat", mIsFlat);
		add_process(' ', "Lines", mLines);
		
		// init mesh menu
		sub_menus.push_back("Mesh");
		add_process('r', "Reset Mesh", mResetMesh);
		add_process('w', "Write Mesh", mWriteMesh);
		add_process('<', "Zoom In", mZoomIn);
		add_process('>', "Zoom Out", mZoomOut);
		add_process(27, "Exit", mExit);
		
		
		// init		
		// process sub menus
		int * sub_menu = new int[sub_menus.size()];

		for(index_t sm = 0; sm < sub_menus.size(); sm++)
			sub_menu[sm] = glutCreateMenu(Viewer::menu_process);

		for(auto mp: Viewer::processes)
		{
			char ss[128];
			sprintf(ss, "[%c] %s", mp.first, mp.second.name_function.c_str());
			glutSetMenu(sub_menu[mp.second.sub_menu]);
			glutAddMenuEntry(ss, mp.first);
		}

		int mainMenu = glutCreateMenu(Viewer::menu);
		glutSetMenu(mainMenu);	
		
		glutAddSubMenu("Select Mesh", mesh_menu);
		for(index_t sm = 0; sm < sub_menus.size(); sm++)
			glutAddSubMenu(sub_menus[sm].c_str(), sub_menu[sm]);
		
		glutAttachMenu( GLUT_RIGHT_BUTTON );

		delete [] sub_menu;
	}
	
	void Viewer::initGLSL()
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
		return mesh().color(i);
	}

	void Viewer::update_VBO()
	{
		for(index_t i = 0; i < n_meshes; i++)
			meshes[i].update();
	}
	
	void Viewer::menu(int value)
	{
	
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

	void Viewer::add_mesh(const vector<che *> & _meshes)
	{
		for(che * _mesh: _meshes)
		{
			assert(n_meshes < N_MESHES);
			meshes[n_meshes++].init(_mesh);
		}
		
		angle_t angle = 2 * M_PI / n_meshes;
		vertex_t r = sqrt(n_meshes - 1);
		for(index_t i = 0; i < n_meshes; i++)
		{
			meshes[i].update();
			meshes[i].translate({r * cos(i * angle), r * sin(i * angle), 0});
		}
	}

	void Viewer::keyboard( unsigned char c, int x, int y )
	{
		if(c >= '0' && c <= '9')
		{
			bgc = (c - '0') / 9.;
			glClearColor( bgc, bgc, bgc, 1. );
		}

		mProcess(processes[c].function);
	}

	void Viewer::menu_meshes(int value)
	{
		current = value;
		select_vertices.clear();
		glutSetWindowTitle(mesh()->filename().c_str());	
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
	
	void Viewer::idle()
	{
		camera.idle();
		glutPostRedisplay();
	}
	
	void Viewer::storeViewerState()
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
		mesh()->reload();
		mesh()->normalize();
		select_vertices.clear();
		other_vertices.clear();
		vectors.clear();
		
		mesh().debug_info();

		update_VBO();
	}
	
	void Viewer::mWriteMesh()
	{
		string file = mesh()->filename();
		index_t p = file.find_last_of('.');
		file = file.substr(0, p) + "_new.off";
		cout << __FUNCTION__ << " " << file << endl;
		mesh()->write_file(file);
	}
	
	void Viewer::mExit()
	{
	//	storeViewerState();
		glutLeaveMainLoop();
	}
	
	void Viewer::mWireframe()
	{
		renderWireframe = !renderWireframe;
	}
	
	void Viewer::mZoomIn()
	{
		camera.zoomIn();
	}
	
	void Viewer::mZoomOut()
	{
		camera.zoomOut();
	}
	
	void Viewer::mGradientField()
	{
		renderGradientField = !renderGradientField;
	}

	void Viewer::mNormalField()
	{
		renderNormalField = !renderNormalField;
	}
		
	void Viewer::mBorder()
	{
		renderBorder = !renderBorder;
		if(!renderBorder) select_vertices.clear();
	}
	
	void Viewer::mOrientation()
	{
		mesh().invert_orientation();
		mesh().update_normals();
		update_VBO();
	}

	void Viewer::mIsFlat()
	{
		is_flat = !is_flat;
	}
	
	void Viewer::mLines()
	{
		lines = !lines;
	}
	
	void Viewer::display()
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
	
	void Viewer::setGL()
	{
		glClearColor( bgc, bgc, bgc, 1. );
		setLighting();
	}
	
	void Viewer::setLighting()
	{
		GLfloat position[4] = { 20., 30., 40., 0. };
		glLightfv( GL_LIGHT0, GL_POSITION, position );
		glEnable( GL_LIGHT0 );
		glEnable( GL_NORMALIZE );
	}
	
	void Viewer::setMeshMaterial()
	{
		GLfloat diffuse[4] = { .8, .5, .3, 1. };
		GLfloat specular[4] = { .3, .3, .3, 1. };
		GLfloat ambient[4] = { .2, .2, .5, 1. };
		
		glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,	diffuse );
		glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
		glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,	ambient );
		glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 16.		);
	}
	
	void Viewer::drawScene()
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
		if(render_corr) draw_corr();

		drawIsolatedVertices();
		drawVectors();
		drawSelectedVertices();
/*	
		char text[50];
		int n_text;
		n_text = sprintf(text, "%15s: %llu", "Vertices", mesh()->n_vertices());
		drawText(text, n_text, 50, 100);
		n_text = sprintf(text, "%15s: %llu", "Faces", mesh()->n_faces());
		drawText(text, n_text, 50, 88);
*/		
	//	glPopAttrib();
	}
	
	void Viewer::drawPolygons()
	{
		for(index_t i = 0; i < n_meshes; i++)
			meshes[i].draw();
	}
	
	void Viewer::drawWireframe()
	{
		shader.disable();
		for(index_t i = 0; i < n_meshes; i++)
			meshes[i].draw_wireframe();
	}
	
	void Viewer::drawVectors()
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
	
	void Viewer::drawIsolatedVertices()
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

	void Viewer::draw_corr()
	{
		if(n_meshes < 2) return;
		
		shader.disable();

		glPushAttrib(GL_ALL_ATTRIB_BITS);

		glDisable(GL_LIGHTING);
		glColor3f(0.8, .0, .0);
		glLineWidth(2.0);

		glBegin(GL_LINES);
		for(index_t & v: select_vertices)
		{
			glVertex3v(&meshes[0]->gt(v).x);
			glVertex3v(&meshes[1]->gt(v).x);
		}
		glEnd();

		glPopAttrib();
		
		// spheres corr
		glPushAttrib( GL_ALL_ATTRIB_BITS );
		
		glEnable( GL_COLOR_MATERIAL );
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glColor3f( 0.8, 0.0, 0.0 );
		
		double h = 0.015 * camera.zoom; 
		for(index_t & v: select_vertices)
		{
			glPushMatrix();
			glTranslated(meshes[1]->gt(v).x, meshes[1]->gt(v).y, meshes[1]->gt(v).z);
			glutSolidSphere(h, 10, 10);
			glPopMatrix();
		}

		glEnd();
		
		glPopAttrib();

	}

	void Viewer::drawVertices()
	{
		for(index_t v = 0; v < mesh()->n_vertices(); v++)
		{
			glLoadName(v);
			glBegin(GL_POINTS);
			glVertex3v( &mesh()->gt(v).x );
			glEnd();
		}
	}
	
	void Viewer::drawBorder()
	{
		select_vertices.clear();
		for(index_t b = 0; b < mesh()->n_borders(); b++)
			for_border(he, mesh(), mesh()->bt(b))
				select_vertices.push_back(mesh()->vt(he));
	}
	
	void Viewer::drawSelectedVertices()
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
			glTranslated(mesh()->gt(v).x, mesh()->gt(v).y, mesh()->gt(v).z);
			glutSolidSphere(h, 10, 10);
			glPopMatrix();
		}

		glEnd();
		
		glPopAttrib();
	}

	void Viewer::drawNormalField()
	{
		shader.disable();
		for(index_t i = 0; i < n_meshes; i++)
			meshes[i].draw_normal_field();
	}

	void Viewer::drawGradientField()
	{
		shader.disable();
		for(index_t i = 0; i < n_meshes; i++)
			meshes[i].draw_gradient_field();
	}

	void Viewer::pickVertex(int x, int y)
	{
		int width = glutGet(GLUT_WINDOW_WIDTH);
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		if( x < 0 || x >= width || y < 0 || y >= height ) return;
		
		int bufSize = mesh()->n_vertices();
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

		if (index >= 0 && index < mesh()->n_vertices())
		{
			debug(index)
			debug(mesh().color(index))
			select_vertices.push_back(index);
		}
	}
}

