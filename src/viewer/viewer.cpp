#include "viewer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <thread>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "che_off.h"
#include "che_obj.h"
#include "che_ply.h"

#include "CImg.h"

using namespace cimg_library;
using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


const int viewer::m_window_size[N_MESHES][2] = {{1, 1}, {1, 2}, {1, 3}, 
												{2, 2}, {2, 3}, {2, 3},
												{2, 4}, {2, 4}, {2, 5},
												{2, 5}, {3, 4}, {3, 4}};


viewer::viewer()
{
	window = nullptr;
	
	n_meshes = current = 0;

	render_opt = 0;
	r_embree = nullptr;
	frame = nullptr;

	render_wireframe = false;
	render_gradient_field = false;
	render_normal_field = false;
	render_border = false;
	render_lines = false;
	render_corr = false;
	render_flat = false;

	bgc = 0;

	init_gl();
	init_imgui();
	init_menus();
	
	init_glsl();
	
	info_gl();
}

viewer::~viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	if(r_embree) delete r_embree;
	if(frame) delete [] frame;
}

bool viewer::run()
{
	while(!glfwWindowShouldClose(window))
	{
		glfwGetFramebufferSize(window, &viewport_width, &viewport_height);
		viewport_width /= m_window_size[n_meshes - 1][1];
		viewport_height /= m_window_size[n_meshes - 1][0];

		eye		= vertex(0., 0., -2. * cam.zoom);
		center	= vertex(0., 0., 0.);
		up		= vertex(0., 1., 0.);
		
		light = vertex(-1., 1., -2.);

		quaternion r = cam.currentRotation();

		eye = r.conj() * eye * r;
		light = r.conj() * light * r;
		
		proj_mat = glm::perspective(45.f, float(viewport_width) / float(viewport_height), .01f, 1000.f);
		view_mat = glm::lookAt(	glm::vec3(eye[1], eye[2], eye[3]), 
								glm::vec3(center[1], center[2], center[3]), 
								glm::vec3(up[1], up[2], up[3])
								);

		// RENDER 
		switch(render_opt)
		{
			case 1: render_embree(); break;
			case 2: render_optix(); break;
			default: render_gl();
		}
		
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		if(ImGui::BeginMainMenuBar())
    	{
        	if(ImGui::BeginMenu("Select"))
			{
				for(index_t i = 0; i < n_meshes; i++)
					if(ImGui::MenuItem((to_string(i) + ". " + meshes[i]->name()).c_str()))
						current = i;	

            	ImGui::EndMenu();
			}

			for(index_t i = 0; i < sub_menus.size(); i++)
			{
        		if(ImGui::BeginMenu(sub_menus[i].c_str()))
        		{
					for(auto & p: processes)
            			if(p.second.function != nullptr && p.second.sub_menu == i &&
							ImGui::MenuItem(p.second.name.c_str(), ('[' + p.second.key + ']').c_str()))
							p.second.function(this);

            		ImGui::EndMenu();
				}
        	}
        	ImGui::EndMainMenuBar();
		}

		// Rendering
		ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return true;
}

che_viewer & viewer::mesh()
{
	return meshes[current];
}

void viewer::info_gl()
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

void error_callback(int error, const char* description)
{
	fprintf(stderr, "Error %d: %s\n", error, description);
}

void viewer::init_gl()
{
	glfwSetErrorCallback(error_callback);

	glfwInit();

	#ifdef __APPLE__
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
		glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
	#endif
		
	window = glfwCreateWindow(1600, 900, "gproshan", NULL, NULL);

	glfwSetWindowUserPointer(window, this);

	glfwSetKeyCallback(window, keyboard_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwSetScrollCallback(window, scroll_callback);

	glfwMakeContextCurrent(window);
	glfwSwapInterval(0);

	glewInit();
}

void viewer::init_imgui()
{
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void) io;
	
	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");
}

void viewer::init_menus()
{
	// init viewer menu
	sub_menus.push_back("Viewer");
	add_process(GLFW_KEY_F1, {"F1", "Help", menu_help});
	add_process(GLFW_KEY_F2, {"F2", "Invert Orientation", invert_orientation});
	add_process(GLFW_KEY_F3, {"F3", "Render Wireframe", set_render_wireframe});
	add_process(GLFW_KEY_F4, {"F4", "Gradient Field", set_render_gradient_field});
	add_process(GLFW_KEY_F5, {"F5", "Normal Field", set_render_normal_field});
	add_process(GLFW_KEY_F6, {"F6", "Show Borders", set_render_border});
	add_process(GLFW_KEY_F7, {"F7", "Raycasting", raycasting});
	add_process(GLFW_KEY_F8, {"F8", "Render GL", set_render_gl});
	add_process(GLFW_KEY_F9, {"F9", "Render Embree", set_render_embree});
	add_process(GLFW_KEY_F10, {"F10", "Render OptiX", set_render_optix});
	add_process(GLFW_KEY_SPACE, {"SPACE", "Level Curves", set_render_lines});
	add_process(GLFW_KEY_TAB, {"TAB", "Render Flat", set_is_flat});
//	add_process('+', "Show Corr", set_render_corr);

	// init mesh menu
	sub_menus.push_back("Mesh");
	add_process(GLFW_KEY_BACKSPACE, {"BACKSPACE", "Reload/Reset", menu_reset_mesh});
	add_process(GLFW_KEY_W, {"W", "Save Mesh", menu_save_mesh});
	add_process(GLFW_KEY_UP, {"UP", "Zoom in", menu_zoom_in});
	add_process(GLFW_KEY_DOWN, {"DOWN", "Zoom out", menu_zoom_out});
	add_process(GLFW_KEY_RIGHT, {"RIGHT", "Background color inc", menu_bgc_inc});
	add_process(GLFW_KEY_LEFT, {"LEFT", "Background color dec", menu_bgc_dec});
	add_process(GLFW_KEY_1, {"1", "Background color white", menu_bgc_white});
	add_process(GLFW_KEY_0, {"0", "Background color black", menu_bgc_black});
}

void viewer::init_glsl()
{
	shader_program.load_vertex("../shaders/vertex.glsl");
	shader_program.load_fragment("../shaders/fragment.glsl");
}

void viewer::update_vbo()
{
	for(index_t i = 0; i < n_meshes; i++)
		meshes[i].update();
}

void viewer::menu_process(int value)
{
	menu_process(processes[value].function);
}

void viewer::add_process(const int & key, const process_t & process)
{
	if(processes.find(key) == processes.end())
	{
		processes[key] = process;
		processes[key].sub_menu = sub_menus.size() - 1;
	}
	else cerr << "Repeat key: " << key << endl;
}

void viewer::add_mesh(const vector<che *> & _meshes)
{
	for(che * _mesh: _meshes)
	{
		assert(n_meshes < N_MESHES);
		meshes[n_meshes++].init(_mesh);
		meshes[n_meshes - 1].log_info();
	}
	
	if(!n_meshes) return;

	glfwSetWindowTitle(window, mesh()->filename().c_str());

	const int * mw = m_window_size[n_meshes - 1];

	index_t m = n_meshes - 1;
	for(int i = mw[1] - 1; i >= 0; i--)
	for(int j = 0; j < mw[0]; j++)
	{
		meshes[m].vx = i;
		meshes[m].vy = j;
		if(!m--) return;
	}
}

void viewer::keyboard_callback(GLFWwindow * window, int key, int scancode, int action, int mods)
{
	if(action == GLFW_RELEASE) return;
	
	if(key == GLFW_KEY_ESCAPE)
	{
		glfwSetWindowShouldClose(window, GLFW_TRUE);
		return;
	}
	
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	view->menu_process(view->processes[key].function);
	glClearColor(view->bgc, view->bgc, view->bgc, 1.);
}

void viewer::menu_meshes(int value)
{
	current = value;
	select_vertices.clear();
	glfwSetWindowTitle(window, mesh()->filename().c_str());
}

void viewer::mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	
	if(mods == GLFW_MOD_SHIFT && action == GLFW_RELEASE)
		view->pick_vertex(xpos, ypos);
	else view->cam.mouse(button, action, xpos, ypos);
}

void viewer::cursor_callback(GLFWwindow * window, double x, double y)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS)
	{
		viewer * view = (viewer *) glfwGetWindowUserPointer(window);
		view->cam.motion(x, y);
	}
}

void viewer::scroll_callback(GLFWwindow * window, double xoffset, double yoffset)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	
	if(yoffset > 0) view->cam.zoomIn();
	if(yoffset < 0) view->cam.zoomOut();
}

void viewer::idle()
{
	//cam.idle();
	////glutPostRedisplay();
}

void viewer::menu_process(function_t pro)
{
	if(pro) pro(this);

	update_vbo();
}

void viewer::menu_help(viewer * view)
{
	for(auto & p: view->processes)
		if(p.second.function != nullptr)
			fprintf(stderr, "%16s: %s\n", ('[' + p.second.key + ']').c_str(), p.second.name.c_str());
}

void viewer::menu_reset_mesh(viewer * view)
{
	view->select_vertices.clear();
	view->other_vertices.clear();
	view->vectors.clear();

	view->mesh().reload();

	view->update_vbo();
}

void viewer::menu_save_mesh(viewer * view)
{
	gproshan_log(APP_VIEWER);
	
	gproshan_log(format: [off obj ply]);
	
	string format; cin >> format;
	string file = view->mesh()->filename() + "_new";
	
	if(format == "off") che_off::write_file(view->mesh(), file);
	if(format == "obj") che_obj::write_file(view->mesh(), file);
	if(format == "ply") che_ply::write_file(view->mesh(), file);

	cerr << "saved: " << file + "." + format << endl;
}

void viewer::menu_zoom_in(viewer * view)
{
	view->cam.zoomIn();
}

void viewer::menu_zoom_out(viewer * view)
{
	view->cam.zoomOut();
}

void viewer::menu_bgc_inc(viewer * view)
{
	if(view->bgc < 1) view->bgc += 0.05;
}

void viewer::menu_bgc_dec(viewer * view)
{
	if(view->bgc > 0) view->bgc -= 0.05;
}

void viewer::menu_bgc_white(viewer * view)
{
	view->bgc = 1;
}

void viewer::menu_bgc_black(viewer * view)
{
	view->bgc = 0;
}

void viewer::invert_orientation(viewer * view)
{
	view->mesh().invert_orientation();
	view->mesh().update_normals();
	view->update_vbo();
}

void viewer::set_render_gl(viewer * view)
{
	view->render_opt = 0;
}

void viewer::set_render_embree(viewer * view)
{
	view->render_opt = 1;
}

void viewer::set_render_optix(viewer * view)
{
	view->render_opt = 2;
}

void viewer::set_render_wireframe(viewer * view)
{
	view->render_wireframe = !view->render_wireframe;
}

void viewer::set_render_gradient_field(viewer * view)
{
	view->render_gradient_field = !view->render_gradient_field;
}

void viewer::set_render_normal_field(viewer * view)
{
	view->render_normal_field = !view->render_normal_field;
}

void viewer::set_render_border(viewer * view)
{
	view->render_border = !view->render_border;
	if(!view->render_border) view->select_vertices.clear();
}

void viewer::set_render_lines(viewer * view)
{
	view->render_lines = !view->render_lines;
}

void viewer::set_render_corr(viewer * view)
{
	view->render_corr = !view->render_corr;
}

void viewer::set_is_flat(viewer * view)
{
	view->render_flat = !view->render_flat;
}

void viewer::raycasting(viewer * view)
{
	gproshan_log(VIEWER);

	embree rc;
	
	rc.add_mesh(view->mesh());
	rc.build_bvh();
	
	float * frame = rc.raycaster(	glm::uvec2(view->viewport_width, view->viewport_height),
									view->view_mat, view->proj_mat	
									); 

	std::thread([](CImg<float> img)
	{
		img.display();
	},
	CImg<float>(frame, view->viewport_width, view->viewport_height)).detach();

	delete [] frame;
}

void viewer::render_gl()
{
	shader_program.enable();

	GLint uniformEye = glGetUniformLocation(shader_program, "eye");
	glUniform3f(uniformEye, eye[1], eye[2], eye[3]);

	GLint uniformLight = glGetUniformLocation(shader_program, "light");
	glUniform3f(uniformLight, light[1], light[2], light[3]);

	GLint uniform_render_flat = glGetUniformLocation(shader_program, "render_flat");
	glUniform1i(uniform_render_flat, render_flat);

	GLint uniform_render_lines = glGetUniformLocation(shader_program, "render_lines");
	glUniform1i(uniform_render_lines, render_lines);

	GLint uniform_model_view_mat = glGetUniformLocation(shader_program, "model_view_mat");
	glUniformMatrix4fv(uniform_model_view_mat, 1, 0, &view_mat[0][0]);
	
	GLint uniform_proj_mat = glGetUniformLocation(shader_program, "proj_mat");
	glUniformMatrix4fv(uniform_proj_mat, 1, 0, &proj_mat[0][0]);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	
	draw_scene();
	
	shader_program.disable();
}

void viewer::render_embree()
{
	if(r_embree == nullptr)
	{
		double time_build_embree;
		TIC(time_build_embree);

			r_embree = new embree;
			r_embree->add_mesh(mesh());
			r_embree->build_bvh();

			frame_size = viewport_width * viewport_height;
			frame = new float[frame_size];

		TOC(time_build_embree);
		gproshan_log_var(time_build_embree);
	}

	if(frame_size < viewport_width * viewport_height)
	{
		delete [] frame;
		
		frame_size = viewport_width * viewport_height;
		frame = new float[frame_size];
	}
	
	r_embree->raytracing(frame, glm::uvec2(viewport_width, viewport_height), view_mat, proj_mat, 1);	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDrawPixels(viewport_width, viewport_height, GL_LUMINANCE, GL_FLOAT, frame);
}

void viewer::render_optix()
{
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

	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);
		meshes[i].draw();
	}
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
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);
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
		////glutSolidSphere(h, 10, 10);
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
//		//glutSolidSphere(h, 10, 10);
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
		glViewport(mesh().vx * viewport_width, mesh().vy * viewport_height, viewport_width, viewport_height);
		
		glPushMatrix();
		glTranslated(mesh()->gt(v).x, mesh()->gt(v).y, mesh()->gt(v).z);
	//	//glutSolidSphere(h, 10, 10);
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
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);
		meshes[i].draw_normal_field();
	}
}

void viewer::draw_gradient_field()
{
	shader_program.disable();
	
	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);
		meshes[i].draw_gradient_field();
	}
}

void viewer::pick_vertex(int x, int y)
{
	gproshan_log(VIEWER);

	int width, height;
	glfwGetFramebufferSize(window, &width, &height);

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
	//gluPickMatrix(x, viewport[3] - y, 10, 10, viewport);
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
		
		gproshan_log_var(index);
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

	//viewport_heightile(*str) //glutBitmapCharacter(font, *str++);

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glPopAttrib();
}


} // namespace gproshan

