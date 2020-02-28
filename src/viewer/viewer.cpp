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
#include "che_sphere.h"

#include <CImg.h>


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
	render_frame = nullptr;

	#ifdef GPROSHAN_EMBREE
		rt_embree = nullptr;
	#endif // GPROSHAN_EMBREE

	#ifdef GPROSHAN_OPTIX
		rt_optix = nullptr;
	#endif // GPROSHAN_OPTIX

	render_wireframe = false;
	render_wireframe_fill = false;
	render_gradient_field = false;
	render_normal_field = false;
	render_border = false;
	render_lines = false;
	render_flat = false;

	bgc = 0;

	action = false;

	init_gl();
	init_glsl();
	init_imgui();
	init_menus();	
	
	info_gl();
	
	sphere.init(new che_sphere(0.01), false);
}

viewer::~viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	#ifdef GPROSHAN_EMBREE
		if(rt_embree) delete rt_embree;
	#endif // GPROSHAN_EMBREE

	#ifdef GPROSHAN_OPTIX
		if(rt_optix) delete rt_optix;
	#endif // GPROSHAN_OPTIX

	if(render_frame) delete render_frame;
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
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
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
						if(	p.second.function != nullptr &&
							p.second.sub_menu == i &&
							ImGui::MenuItem(p.second.name.c_str(), ('[' + p.second.key + ']').c_str()) )
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
		//glfwWaitEvents();
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
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#endif
		
	window = glfwCreateWindow(1600, 900, "gproshan", NULL, NULL);

	glfwSetWindowUserPointer(window, this);

	glfwSetKeyCallback(window, keyboard_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwSetScrollCallback(window, scroll_callback);

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	glewInit();


	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void viewer::init_imgui()
{
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void) io;
	
	ImGui::StyleColorsDark();

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 410");
}

void viewer::init_menus()
{
	sub_menus.push_back("Viewer");
	add_process(GLFW_KEY_F1, {"F1", "Help", menu_help});
	add_process(GLFW_KEY_UP, {"UP", "Zoom in", menu_zoom_in});
	add_process(GLFW_KEY_DOWN, {"DOWN", "Zoom out", menu_zoom_out});
	add_process(GLFW_KEY_RIGHT, {"RIGHT", "Background color inc", menu_bgc_inc});
	add_process(GLFW_KEY_LEFT, {"LEFT", "Background color dec", menu_bgc_dec});
	add_process(GLFW_KEY_1, {"1", "Background color white", menu_bgc_white});
	add_process(GLFW_KEY_0, {"0", "Background color black", menu_bgc_black});

	sub_menus.push_back("Render");
	add_process(GLFW_KEY_F6, {"F6", "Render Wireframe", set_render_wireframe});
	add_process(GLFW_KEY_F7, {"F7", "Render Wireframe Fill", set_render_wireframe_fill});
	add_process(GLFW_KEY_F8, {"F8", "Render GL", set_render_gl});
	add_process(GLFW_KEY_F9, {"F9", "Render Embree", set_render_embree});
	add_process(GLFW_KEY_F10, {"F10", "Render OptiX", set_render_optix});
	add_process(GLFW_KEY_ENTER, {"ENTER", "Raycasting", raycasting});

	sub_menus.push_back("Mesh");
	add_process(GLFW_KEY_BACKSPACE, {"BACKSPACE", "Reload/Reset", menu_reset_mesh});
	add_process(GLFW_KEY_TAB, {"TAB", "Render Flat", set_render_flat});
	add_process(GLFW_KEY_SPACE, {"SPACE", "Level Curves", set_render_lines});
	add_process(GLFW_KEY_F2, {"F2", "Invert Orientation", invert_orientation});
	add_process(GLFW_KEY_F3, {"F3", "Gradient Field", set_render_gradient_field});
	add_process(GLFW_KEY_F4, {"F4", "Normal Field", set_render_normal_field});
	add_process(GLFW_KEY_F5, {"F5", "Select Border Vertices", set_render_border});
	add_process(GLFW_KEY_W, {"W", "Save Mesh", menu_save_mesh});
}

void viewer::init_glsl()
{
	shader_sphere.load_vertex("../shaders/vertex_sphere.glsl");
	shader_sphere.load_fragment("../shaders/fragment_sphere.glsl");

	shader_program.load_vertex("../shaders/vertex.glsl");
	shader_program.load_geometry("../shaders/geometry.glsl");
	shader_program.load_fragment("../shaders/fragment.glsl");
	
	shader_normals.load_vertex("../shaders/vertex_normals.glsl");
	shader_normals.load_geometry("../shaders/geometry_normals.glsl");
	shader_normals.load_fragment("../shaders/fragment_normals.glsl");
	
	shader_gradient.load_vertex("../shaders/vertex_gradient.glsl");
	shader_gradient.load_geometry("../shaders/geometry_gradient.glsl");
	shader_gradient.load_fragment("../shaders/fragment_gradient.glsl");
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
		view->action = true;
	}
}

void viewer::scroll_callback(GLFWwindow * window, double xoffset, double yoffset)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	
	if(yoffset > 0)
	{
		view->cam.zoomIn();
		view->action = true;
	}

	if(yoffset < 0)
	{
		view->cam.zoomOut();
		view->action = true;
	}
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

void viewer::set_render_wireframe_fill(viewer * view)
{
	view->render_wireframe_fill = !view->render_wireframe_fill;
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

void viewer::set_render_flat(viewer * view)
{
	view->render_flat = !view->render_flat;
}

void viewer::raycasting(viewer * view)
{
#ifdef GPROSHAN_EMBREE

	gproshan_log(VIEWER);

	rt::embree rc({view->mesh()});
	
	float * frame = rc.raycaster(	glm::uvec2(view->viewport_width, view->viewport_height),
									view->view_mat, view->proj_mat	
									); 

	std::thread([](CImg<float> img)
	{
		img.display();
	},
	CImg<float>(frame, view->viewport_width, view->viewport_height)).detach();

	delete [] frame;

#endif // GPROSHAN_EMBREE
}

void viewer::render_gl()
{
	glProgramUniform3f(shader_sphere, shader_sphere("eye"), eye[1], eye[2], eye[3]);
	glProgramUniform3f(shader_sphere, shader_sphere("light"), light[1], light[2], light[3]);
	glProgramUniformMatrix4fv(shader_sphere, shader_sphere("model_view_mat"), 1, 0, &view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_sphere, shader_sphere("proj_mat"), 1, 0, &proj_mat[0][0]);
	glProgramUniform1f(shader_sphere, shader_sphere("scale"), cam.zoom);


	glProgramUniform3f(shader_program, shader_program("eye"), eye[1], eye[2], eye[3]);
	glProgramUniform3f(shader_program, shader_program("light"), light[1], light[2], light[3]);
	glProgramUniform1i(shader_program, shader_program("render_flat"), render_flat);
	glProgramUniform1i(shader_program, shader_program("render_lines"), render_lines);
	glProgramUniform1i(shader_program, shader_program("render_wireframe"), render_wireframe_fill);
	glProgramUniformMatrix4fv(shader_program, shader_program("model_view_mat"), 1, 0, &view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_program, shader_program("proj_mat"), 1, 0, &proj_mat[0][0]);

	
	shader_program.enable();
	draw_scene();


	if(render_normal_field)
	{
		glProgramUniform1f(shader_normals, shader_normals("length"), mesh().factor);
		glProgramUniformMatrix4fv(shader_normals, shader_normals("model_view_mat"), 1, 0, &view_mat[0][0]);
		glProgramUniformMatrix4fv(shader_normals, shader_normals("proj_mat"), 1, 0, &proj_mat[0][0]);
		
		shader_normals.enable();
		draw_scene();
	}
	
	
	if(render_gradient_field)
	{
		glProgramUniform1f(shader_gradient, shader_gradient("length"), mesh().factor);
		glProgramUniformMatrix4fv(shader_gradient, shader_gradient("model_view_mat"), 1, 0, &view_mat[0][0]);	
		glProgramUniformMatrix4fv(shader_gradient, shader_gradient("proj_mat"), 1, 0, &proj_mat[0][0]);
		
		shader_gradient.enable();
		draw_scene();
	}
}

void viewer::render_embree()
{
#ifdef GPROSHAN_EMBREE

	if(!rt_embree)
	{
		double time_build_embree;
		TIC(time_build_embree);

			rt_embree = new rt::embree({mesh()});

		TOC(time_build_embree);
		gproshan_log_var(time_build_embree);
	}

	rt_embree->pathtracing(	glm::uvec2(viewport_width, viewport_height),
							view_mat, proj_mat, {glm::vec3(light[1], light[2], light[3])}, action);
	
	action = false;

	if(!render_frame) render_frame = new frame;	
	render_frame->display(viewport_width, viewport_height, rt_embree->img);

#endif // GPROSHAN_EMBREE
}

void viewer::render_optix()
{
#ifdef GPROSHAN_OPTIX

	if(!rt_optix)
	{
		double time_build_optix;
		TIC(time_build_optix);

			rt_optix = new rt::optix({mesh()});

		TOC(time_build_optix);
		gproshan_log_var(time_build_optix);
	}

	rt_optix->pathtracing(	glm::uvec2(viewport_width, viewport_height),
							view_mat, proj_mat, {glm::vec3(light[1], light[2], light[3])}, action);
	
	action = false;
	
	if(!render_frame) render_frame = new frame;	
	render_frame->display(viewport_width, viewport_height, rt_embree->img);

#endif // GPROSHAN_OPTIX
}

void viewer::draw_scene()
{	
	draw_polygons();

	if(render_border)
		draw_border();

	draw_selected_vertices();
}

void viewer::draw_polygons()
{
	if(render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


	for(index_t i = 0; i < n_meshes; i++)
	{
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);
		meshes[i].draw();
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
	if(sphere_translations.size() != select_vertices.size())
	{
		sphere_translations.resize(select_vertices.size());

		for(index_t i = 0; i < select_vertices.size(); i++)
			sphere_translations[i] = mesh()->gt(select_vertices[i]);

		sphere.update_instances_translations(sphere_translations);
	}
	
	glViewport(mesh().vx * viewport_width, mesh().vy * viewport_height, viewport_width, viewport_height);
	
	shader_sphere.enable();
	if(sphere_translations.size())
		sphere.draw();
	shader_sphere.disable();
}

void viewer::pick_vertex(int x, int y)
{
	gproshan_log(VIEWER);
}


} // namespace gproshan

