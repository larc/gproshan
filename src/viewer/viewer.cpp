#include "viewer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>

#include <che_off.h>
#include <che_obj.h>
#include <che_ply.h>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


// declare static member variables
che_viewer viewer::meshes[N_MESHES];
vcorr_t viewer::corr_mesh[N_MESHES]; // zero initialization
size_t viewer::n_meshes = 0;
index_t viewer::current = 0;

vector<index_t> viewer::select_vertices;
vector<vertex> viewer::other_vertices;
vector<vertex> viewer::vectors;
vector<string> viewer::sub_menus;

char * viewer::share = nullptr;

int viewer::window_size[2] = {1366, 768};
int viewer::m_window_size[N_MESHES][2] = {	{1, 1}, {1, 2}, {1, 3}, 
											{2, 2}, {2, 3}, {2, 3},
											{2, 4}, {2, 4}, {2, 5},
											{2, 5}, {3, 4}, {3, 4} };
double viewer::ww = 0;
double viewer::wh = 0;

camera viewer::cam;
shader viewer::shader_program;
bool viewer::render_wireframe = false;
bool viewer::render_gradient_field = false;
bool viewer::render_normal_field = false;
bool viewer::render_border = false;
bool viewer::render_corr = false;
bool viewer::render_lines = false;
bool viewer::is_flat = false;
float viewer::bgc = 0.;

map<unsigned char, process_t> viewer::processes;

const int & viewer::window_width()
{
	return window_size[0];
}

const int & viewer::window_height()
{
	return window_size[1];
}

che_viewer & viewer::mesh()
{
	assert(n_meshes > 0);
	return meshes[current];
}

void viewer::init(const vector<che *> & _meshes)
{
//	restoreviewerState();

	init_glut();
	add_mesh(_meshes);

	glutSetWindowTitle(mesh()->filename().c_str());
	init_menus();

	debug_info();
	mesh().log_info();

	set_gl();
	init_glsl();

	glutMainLoop();
}

void viewer::debug_info()
{
	GLint major, minor;
	glGetIntegerv(GL_MAJOR_VERSION, &major);
	glGetIntegerv(GL_MINOR_VERSION, &minor);
	fprintf(stderr, "GL Vendor %s\n", glGetString(GL_VENDOR));
	fprintf(stderr, "GL Renderer %s\n", glGetString(GL_RENDERER));
	fprintf(stderr, "GL Version (string) %s\n", glGetString(GL_VERSION));
	fprintf(stderr, "GL Version (integer) %d.%d\n", major, minor);
	fprintf(stderr, "GLSL Version %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
}

void viewer::init_glut()
{
	int argc = 0;
	vector< vector<char> > argv(1);

	// initialize window
	glutInitWindowSize(window_size[0], window_size[1]);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInit(&argc, (char**)&argv);
	glutCreateWindow("gproshan");

	glewInit();

	// specify callbacks
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
}

void viewer::init_menus()
{
	// set current mesh menu
	int mesh_menu = glutCreateMenu(viewer::menu_meshes);
	glutSetMenu(mesh_menu);
	for(index_t i = 0; i < n_meshes; i++)
		glutAddMenuEntry(meshes[i]->filename().c_str(), i);

	// init viewer menu
	sub_menus.push_back("viewer");
	add_process('i', "Invert Orientation", invert_orientation);
	add_process('f', "Wireframe", set_render_wireframe);
	add_process('g', "Gradient Field", set_render_gradient_field);
	add_process('n', "Normal Field", set_render_normal_field);
	add_process('b', "Show Border", set_render_border);
	add_process(' ', "Lines", set_render_lines);
	add_process('+', "Show Corr", set_render_corr);
	add_process('\t', "Flat", set_is_flat);

	// init mesh menu
	sub_menus.push_back("Mesh");
	add_process('r', "Reset Mesh", menu_reset_mesh);
	add_process('w', "Write Mesh", menu_save_mesh);
	add_process('<', "Zoom In", menu_zoom_in);
	add_process('>', "Zoom Out", menu_zoom_out);
	add_process(27, "Exit", menu_exit);


	// init
	// process sub menus
	int * sub_menu = new int[sub_menus.size()];

	for(index_t sm = 0; sm < sub_menus.size(); sm++)
		sub_menu[sm] = glutCreateMenu(viewer::menu_process);

	char ss[128];
	for(auto mp: viewer::processes)
	{
		sprintf(ss, "[%c] %s", mp.first, mp.second.name_function.c_str());
		glutSetMenu(sub_menu[mp.second.sub_menu]);
		glutAddMenuEntry(ss, mp.first);
	}

	int mainMenu = glutCreateMenu(viewer::menu);
	glutSetMenu(mainMenu);

	glutAddSubMenu("Select Mesh", mesh_menu);
	for(index_t sm = 0; sm < sub_menus.size(); sm++)
		glutAddSubMenu(sub_menus[sm].c_str(), sub_menu[sm]);

	glutAttachMenu(GLUT_RIGHT_BUTTON);

	delete [] sub_menu;
}

void viewer::init_glsl()
{
	shader_program.load_vertex("../shaders/new_vertex.glsl");
	shader_program.load_fragment("../shaders/new_fragment.glsl");
}

color_t & viewer::vcolor(const index_t & i)
{
	return mesh().color(i);
}

void viewer::update_vbo()
{
	for(index_t i = 0; i < n_meshes; i++)
		meshes[i].update();
}

void viewer::menu(int )
{
}

void viewer::menu_process(int value)
{
	menu_process(processes[value].function);
}

void viewer::add_process(const char & key, const string & name, function_t function)
{
	if(processes.find(key) == processes.end())
	{
		processes[key] = {(index_t) sub_menus.size() - 1, name, function};
	}
	else cerr << "Repeat key: " << key << endl;
}

void viewer::add_mesh(const vector<che *> & _meshes)
{
	for(che * _mesh: _meshes)
	{
		assert(n_meshes < N_MESHES);
		meshes[n_meshes++].init(_mesh);
	}
	
	int * mw = m_window_size[n_meshes - 1];

	index_t m = n_meshes - 1;
	for(int i = mw[1] - 1; i >= 0; i--)
	for(int j = 0; j < mw[0]; j++)
	{
		meshes[m].vx = i;
		meshes[m].vy = j;
		if(!m--) return;
	}
}

void viewer::keyboard(unsigned char c, int , int )
{
	if(c >= '0' && c <= '9')
	{
		bgc = (c - '0') / 9.;
		glClearColor(bgc, bgc, bgc, 1.);
	}

	menu_process(processes[c].function);
}

void viewer::menu_meshes(int value)
{
	current = value;
	select_vertices.clear();
	glutSetWindowTitle(mesh()->filename().c_str());
}

void viewer::special(int i, int , int )
{
	switch(i)
	{
		case GLUT_KEY_UP:
			cam.zoomIn();
			break;
		case GLUT_KEY_DOWN:
			cam.zoomOut();
			break;
		case 27:
			menu_exit();
			break;
		default:
			break;
	}
}

void viewer::mouse(int button, int state, int x, int y)
{
	if((glutGetModifiers() & GLUT_ACTIVE_SHIFT) and state == GLUT_UP)
		pick_vertex(x, y);
	else if(button == 6 || button == 4) cam.zoomIn();
	else if(button == 5 || button == 3) cam.zoomOut();
	else cam.mouse(button, state, x, y);
}

void viewer::motion(int x, int y)
{
	cam.motion(x, y);
}

void viewer::idle()
{
	cam.idle();
	glutPostRedisplay();
}

void viewer::menu_process(function_t pro)
{
	if(pro) pro();
	update_vbo();
}

void viewer::menu_reset_mesh()
{
	select_vertices.clear();
	other_vertices.clear();
	vectors.clear();

	mesh().reload();
//	mesh().debug_info();

	update_vbo();
}

void viewer::menu_save_mesh()
{
	gproshan_log(APP_VIEWER);
	
	gproshan_log(format: [off obj ply]);
	
	string format; cin >> format;
	string file = mesh()->filename() + "_new";
	
	if(format == "off") che_off::write_file(mesh(), file);
	if(format == "obj") che_obj::write_file(mesh(), file);
	if(format == "ply") che_ply::write_file(mesh(), file);

	cerr << "saved: " << file + "." + format << endl;
}

void viewer::menu_exit()
{
	glutLeaveMainLoop();
}

void viewer::menu_zoom_in()
{
	cam.zoomIn();
}

void viewer::menu_zoom_out()
{
	cam.zoomOut();
}

void viewer::invert_orientation()
{
	mesh().invert_orientation();
	mesh().update_normals();
	update_vbo();
}

void viewer::set_render_wireframe()
{
	render_wireframe = !render_wireframe;
}

void viewer::set_render_gradient_field()
{
	render_gradient_field = !render_gradient_field;
}

void viewer::set_render_normal_field()
{
	render_normal_field = !render_normal_field;
}

void viewer::set_render_border()
{
	render_border = !render_border;
	if(!render_border) select_vertices.clear();
}

void viewer::set_render_lines()
{
	render_lines = !render_lines;
}

void viewer::set_render_corr()
{
	render_corr = !render_corr;
}

void viewer::set_is_flat()
{
	is_flat = !is_flat;
}

void viewer::display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	shader_program.enable();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	double aspect = (double) viewport[2] / (double) viewport[3];
	const double fovy = 50.;
	const double clipNear = .01;
	const double clipFar = 1000.;
	gluPerspective(fovy, aspect, clipNear, clipFar);
	//glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
	//glOrtho(-1.0, 1.0, -1.0, 1.0, 100, 1000);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();


	quaternion eye = vertex(0., 0., -2.5 * cam.zoom);
	quaternion center = vertex(0., 0., 0.);
	quaternion up = vertex(0., 1., 0.);
	gluLookAt(	eye[1],		eye[2],		eye[3],
			center[1],	center[2],	center[3],
			up[1],		up[2],		up[3]);


	quaternion r = cam.currentRotation();
	eye = r.conj() * eye * r;
	GLint uniformEye = glGetUniformLocation(shader_program, "eye");
	glUniform3f(uniformEye, eye[1], eye[2], eye[3]);

	quaternion light = vertex(-1., 1., -2.);
	light = r.conj() * light * r;
	GLint uniformLight = glGetUniformLocation(shader_program, "light");
	glUniform3f(uniformLight, light[1], light[2], light[3]);

	cam.setView();

	GLint uniformIsFlat = glGetUniformLocation(shader_program, "is_flat");
	glUniform1i(uniformIsFlat, is_flat);

	GLint uniformLines = glGetUniformLocation(shader_program, "lines");
	glUniform1i(uniformLines, render_lines);

	GLfloat ModelViewMatrix[16];
	GLfloat ProjectionMatrix[16];

	glGetFloatv(GL_MODELVIEW_MATRIX, ModelViewMatrix);
	glGetFloatv(GL_PROJECTION_MATRIX, ProjectionMatrix);

	GLint uniformModelViewMatrix = glGetUniformLocation(shader_program, "ModelViewMatrix");
	GLint uniformProjectionMatrix = glGetUniformLocation(shader_program, "ProjectionMatrix");

	glUniformMatrix4fv(uniformModelViewMatrix, 1, 0, ModelViewMatrix);
	glUniformMatrix4fv(uniformProjectionMatrix, 1, 0, ProjectionMatrix);

	//glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	set_mesh_materia();
	
	ww = (double) glutGet(GLUT_WINDOW_WIDTH) / m_window_size[n_meshes - 1][1];
	wh = (double) glutGet(GLUT_WINDOW_HEIGHT) / m_window_size[n_meshes - 1][0];
	draw_scene();

	//glPopAttrib();

	shader_program.disable();
	glutSwapBuffers();
}

void viewer::set_gl()
{
	glClearColor(bgc, bgc, bgc, 1.);
	set_lighting();
}

void viewer::set_lighting()
{
	GLfloat position[4] = {20, 30, 40, 0};
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
}

void viewer::set_mesh_materia()
{
	GLfloat diffuse[4] = { .8, .5, .3, 1. };
	GLfloat specular[4] = { .3, .3, .3, 1. };
	GLfloat ambient[4] = { .2, .2, .5, 1. };

	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,	diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,	ambient);
	glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 16.);
}

void viewer::draw_scene()
{
	draw_polygons();

	if(render_wireframe) draw_wireframe();
	if(render_gradient_field) draw_gradient_field();
	if(render_normal_field) draw_normal_field();
	if(render_border) draw_border();
	if(render_corr) draw_corr();

	draw_isolated_vertices();
	draw_vectors();
	draw_selected_vertices();

//	shader_program.disable();
//	mesh().draw_mesh_info();
//	shader_program.enable();
}

void viewer::draw_polygons()
{
	shader_program.enable();

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * ww, meshes[i].vy * wh, ww, wh);
		meshes[i].draw();
	}
	
	glDisable(GL_POLYGON_OFFSET_FILL);
}

void viewer::draw_wireframe()
{
	shader_program.disable();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor4f(0., 0., 0., 0.4);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * ww, meshes[i].vy * wh, ww, wh);
		meshes[i].draw();
	}
	
	glPopAttrib();
}

void viewer::draw_vectors()
{
	shader_program.disable();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor4f(1., 0., 0., 1);
	glLineWidth(3.0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glBegin(GL_LINES);

	index_t i = 0;
	for(vertex & v: vectors)
	{
		if(i % 8 == 0) glColor4f(1., 0., 0., 1);
		if(i % 8 == 2) glColor4f(0., 1., 0., 1);
		if(i % 8 == 4) glColor4f(0., 0., 1., 1);
		if(i % 8 == 6) glColor4f(1., 1., 0., 1);
		glVertex3v(&v.x);
		i++;
	}

	glEnd();
	glPopAttrib();
}

void viewer::draw_isolated_vertices()
{
	shader_program.disable();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glEnable(GL_COLOR_MATERIAL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor3f(0.5, 0., 0.5);

	double h = 0.01 * cam.zoom;
	for(const vertex & v: other_vertices)
	{
		glPushMatrix();
		glTranslated(v.x, v.y, v.z);
		glutSolidSphere(h, 10, 10);
		glPopMatrix();
	}

	glEnd();

	glPopAttrib();
	/*shader_program.disable();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glPointSize(5);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POINT_SMOOTH);
	glColor3f(1., 0., 0.);

	glBegin(GL_POINTS);

	for(const vertex & v: other_vertices)
		glVertex3v(&v.x);

	glEnd();

	glPopAttrib();*/
}

void viewer::draw_corr()
{
	if(n_meshes < 2) return;
	if(!corr_mesh[current].is_loaded()) return;

	shader_program.disable();

	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_LIGHTING);
	glColor3f(0.8, .0, .0);
	glLineWidth(2.0);

	vertex corr_v;

	glBegin(GL_LINES);
	for(index_t & v: select_vertices)
	{
		corr_v = meshes[corr_mesh[current].mesh_i]->corr_vertex(corr_mesh[current][v]);
		glVertex3v(&meshes[current]->gt(v).x);
		glVertex3v(&corr_v.x);
	}
	glEnd();

	glPopAttrib();

	// spheres corr
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glEnable(GL_COLOR_MATERIAL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor3f(0.8, 0.0, 0.0);

	double h = 0.008 * cam.zoom;
	for(index_t & v: select_vertices)
	{
		glPushMatrix();
		corr_v = meshes[corr_mesh[current].mesh_i]->corr_vertex(corr_mesh[current][v]);
		glTranslated(corr_v.x, corr_v.y, corr_v.z);
		glutSolidSphere(h, 10, 10);
		glPopMatrix();
	}

	glEnd();

	glPopAttrib();

}

void viewer::draw_vertices()
{
	for(index_t v = 0; v < mesh()->n_vertices(); v++)
	{
		glLoadName(v);
		glBegin(GL_POINTS);
		glVertex3v(&mesh()->gt(v).x);
		glEnd();
	}
}

void viewer::draw_border()
{
	select_vertices.clear();
	for(index_t b = 0; b < mesh()->n_borders(); b++)
		for_border(he, mesh(), mesh()->bt(b))
			select_vertices.push_back(mesh()->vt(he));
}

void viewer::draw_selected_vertices()
{
	shader_program.disable();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glEnable(GL_COLOR_MATERIAL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor3f(0., 0.5, 0.5);

	double h = 0.02 * cam.zoom;
	for(int v: select_vertices)
	{
		glViewport(mesh().vx * ww, mesh().vy * wh, ww, wh);
		
		glPushMatrix();
		glTranslated(mesh()->gt(v).x, mesh()->gt(v).y, mesh()->gt(v).z);
		glutSolidSphere(h, 10, 10);
		glPopMatrix();
	}

	glEnd();

	glPopAttrib();
}

void viewer::draw_normal_field()
{
	shader_program.disable();
	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * ww, meshes[i].vy * wh, ww, wh);
		meshes[i].draw_normal_field();
	}
}

void viewer::draw_gradient_field()
{
	shader_program.disable();
	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * ww, meshes[i].vy * wh, ww, wh);
		meshes[i].draw_gradient_field();
	}
}

void viewer::pick_vertex(int x, int y)
{
	gproshan_debug(VIEWER);

	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	if(x < 0 || x >= width || y < 0 || y >= height) return;

	int bufSize = mesh()->n_vertices();
	GLuint * buf = new GLuint[bufSize];
	glSelectBuffer(bufSize, buf);

	GLint viewport[4];
	GLdouble projection[16];
	glGetIntegerv(GL_VIEWPORT, viewport);
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
	draw_vertices();
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	long hits = glRenderMode(GL_RENDER);

	index_t index = -1;
	double min_z = 1.0e100;
	for(long i = 0; i < hits; ++i)
	{
		double distance = buf[4*i + 1];
		if(distance < min_z)
		{
			index = buf[4*i + 3];
			min_z = distance;
		}
	}
	delete[] buf;

	if(index != NIL && index < mesh()->n_vertices())
	{
		select_vertices.push_back(index);
		
		gproshan_debug_var(index);
		gproshan_debug_var(mesh().color(index));
		gproshan_debug_var(mesh()->evt(index));

		if(corr_mesh[current].is_loaded())
			gproshan_error_var(corr_mesh[current][index].alpha);
	}
}

void draw_str(const char * str, int x, int y, float color[4], void * font)
{
	glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);

	glColor4fv(color);
	glRasterPos2i(x, y);

	while(*str) glutBitmapCharacter(font, *str++);

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glPopAttrib();
}


} // namespace gproshan

