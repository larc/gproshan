#include "viewer/viewer.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <vector>
#include <thread>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "mesh/che_off.h"
#include "mesh/che_obj.h"
#include "mesh/che_ply.h"
#include "mesh/che_xyz.h"
#include "mesh/che_sphere.h"

#ifdef GPROSHAN_EMBREE
	#include "raytracing/rt_embree.h"
	#include "raytracing/rt_embree_splat.h"
#endif // GPROSHAN_EMBREE

#ifdef GPROSHAN_OPTIX
	#include "raytracing/rt_optix.h"
#endif // GPROSHAN_OPTIX


#include <CImg.h>


using namespace cimg_library;
using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


const int viewer::m_window_size[N_MESHES + 1][2] = {{1, 1},
													{1, 1}, {1, 2}, {1, 3},
													{2, 2}, {2, 3}, {2, 3},
													{2, 4}, {2, 4}, {2, 5},
													{2, 5}, {3, 4}, {3, 4}
													};


viewer::viewer(int width, int height): window_width(width), window_height(height)
{
	init_gl();
	init_glsl();
	init_imgui();
	init_menus();

	info_gl();
	gproshan_log_var(sizeof(real_t));

	sphere.init(new che_sphere(0.01), false);
}

viewer::~viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	delete rt_embree;
	delete rt_optix;

	delete render_frame;

	delete sphere;
}

bool viewer::run()
{
	while(!glfwWindowShouldClose(window))
	{
		light	= vertex(-1, 1, -2);

		quaternion r = cam.current_rotation();

		light = r.conj() * light * r;

		view_mat = cam.look_at(r);
		proj_mat = glm::perspective(45.f, float(viewport_width) / float(viewport_height), .01f, 1000.f);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		switch(render_opt)
		{
			case R_GL:		render_gl();		break;
		#ifdef GPROSHAN_EMBREE
			case R_EMBREE:	render_embree();	break;
		#endif // GPROSHAN_EMBREE
		#ifdef GPROSHAN_OPTIX
			case R_OPTIX:	render_optix();		break;
		#endif // GPROSHAN_OPTIX
		}

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		if(ImGui::BeginMainMenuBar())
		{
			if(ImGui::BeginMenu("Select"))
			{
				for(index_t i = 0; i < n_meshes; ++i)
					if(ImGui::MenuItem((to_string(i) + ": " + meshes[i]->filename).c_str(), nullptr, i == idx_active_mesh, i != idx_active_mesh))
					{
						idx_active_mesh = i;
						sphere_translations.clear();
						glfwSetWindowTitle(window, active_mesh()->filename.c_str());
					}

				ImGui::EndMenu();
			}

			if(ImGui::BeginMenu("Color"))
			{
				for(index_t i = 0; i < colormap.size(); ++i)
					if(ImGui::MenuItem(colormap[i].c_str(), nullptr, i == idx_colormap, i != idx_colormap))
						idx_colormap = i;

				ImGui::EndMenu();
			}

			for(index_t i = 0; i < sub_menus.size(); ++i)
			{
				if(ImGui::BeginMenu(sub_menus[i].c_str()))
				{
					for(auto & p: processes)
					{
						process_t & pro = p.second;
						if(pro.function != nullptr && pro.sub_menu == i)
							if(ImGui::MenuItem(pro.name.c_str(), ('[' + pro.key + ']').c_str(), &pro.selected))
								sprintf(status_message, "%s", pro.selected ? pro.name.c_str() : "");
					}

					ImGui::EndMenu();
				}
			}

			ImGui::EndMainMenuBar();
		}

		ImGui::SetNextWindowSize(ImVec2(window_width, -1));
		ImGui::SetNextWindowPos(ImVec2(0, window_height - 32));
		ImGui::Begin("status gproshan", nullptr, ImGuiWindowFlags_NoTitleBar);
		ImGui::Text("[] %s", status_message);
		ImGui::SameLine(window_width - 180);
		ImGui::Text("github.com/larc/gproshan");
		ImGui::End();

		ImGui::SetNextWindowSize(ImVec2(320, -1));
		ImGui::SetNextWindowPos(ImVec2(20, 60), ImGuiCond_Once);
		ImGui::Begin("gproshan");

		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5);

		che_viewer & mesh = active_mesh();
		ImGui::Text("%s", mesh->filename.c_str());
		ImGui::Text("%16s: %10lu", "n_vertices", mesh->n_vertices);
		ImGui::Text("%16s: %10lu", "n_faces", mesh->n_faces);

		if(render_pointcloud)
		{
			ImGui::Checkbox("point_normals", &point_normals);
			ImGui::SliderInt("point_size", (int *) &point_size, 1, 32);
		}

		for(auto & p: processes)
		{
			process_t & pro = p.second;
			if(ImGui::CollapsingHeader(("[" + pro.key + "] " + pro.name).c_str(), &pro.selected, ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::PushID(pro.name.c_str());
				ImGui::Indent();
				pro.selected = pro.selected && p.second.function(this);
				ImGui::Unindent();
				ImGui::PopID();
			}
		}

		ImGui::PopItemWidth();
		ImGui::End();

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return true;
}

che_viewer & viewer::active_mesh()
{
	return meshes[idx_active_mesh];
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

	window = glfwCreateWindow(window_width, window_height, "gproshan", NULL, NULL);

	glfwSetWindowUserPointer(window, this);

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);
	glfwSetKeyCallback(window, keyboard_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwSetScrollCallback(window, scroll_callback);

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	glewInit();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glEnable(GL_PROGRAM_POINT_SIZE);
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
	add_process(GLFW_KEY_PERIOD, {"PERIOD", "Save/Load view", menu_save_load_view});
	add_process(GLFW_KEY_UP, {"UP", "Zoom in", menu_zoom_in});
	add_process(GLFW_KEY_DOWN, {"DOWN", "Zoom out", menu_zoom_out});
	add_process(GLFW_KEY_RIGHT, {"RIGHT", "Background color inc", menu_bgc_inc});
	add_process(GLFW_KEY_LEFT, {"LEFT", "Background color dec", menu_bgc_dec});
	add_process(GLFW_KEY_1, {"1", "Background color white", menu_bgc_white});
	add_process(GLFW_KEY_0, {"0", "Background color black", menu_bgc_black});

	sub_menus.push_back("Render");
	add_process(GLFW_KEY_F5, {"F5", "Render Point Cloud", set_render_pointcloud});
	add_process(GLFW_KEY_F6, {"F6", "Render Wireframe", set_render_wireframe});
	add_process(GLFW_KEY_F7, {"F7", "Render Wireframe Fill", set_render_wireframe_fill});
	add_process(GLFW_KEY_F8, {"F8", "Render GL", set_render_gl});
#ifdef GPROSHAN_EMBREE
	add_process(GLFW_KEY_F9, {"F9", "Render Embree", set_render_embree});
	add_process(GLFW_KEY_ENTER, {"ENTER", "Raycasting", raycasting});
#endif // GPROSHAN_EMBREE
#ifdef GPROSHAN_OPTIX
	add_process(GLFW_KEY_F10, {"F10", "Render OptiX", set_render_optix});
#endif // GPROSHAN_OPTIX

	sub_menus.push_back("Mesh");
	add_process(GLFW_KEY_BACKSPACE, {"BACKSPACE", "Reload/Reset", menu_reset_mesh});
	add_process(GLFW_KEY_TAB, {"TAB", "Render Flat", set_render_flat});
	add_process(GLFW_KEY_SPACE, {"SPACE", "Level Curves", set_render_lines});
	add_process(GLFW_KEY_F2, {"F2", "Invert Orientation", invert_orientation});
	add_process(GLFW_KEY_F3, {"F3", "Gradient Field", set_render_gradient_field});
	add_process(GLFW_KEY_F4, {"F4", "Normal Field", set_render_normal_field});
	add_process(GLFW_KEY_APOSTROPHE, {"APOSTROPHE", "Select Border Vertices", set_render_border});
	add_process(GLFW_KEY_W, {"W", "Save Mesh", menu_save_mesh});
}

void viewer::init_glsl()
{
	shader_sphere.load_vertex(shaders_path("vertex_sphere.glsl"));
	shader_sphere.load_fragment(shaders_path("fragment_sphere.glsl"));

	shader_triangles.load_vertex(shaders_path("vertex.glsl"));
	shader_triangles.load_geometry(shaders_path("geometry.glsl"));
	shader_triangles.load_fragment(shaders_path("fragment.glsl"));

	shader_normals.load_vertex(shaders_path("vertex_normals.glsl"));
	shader_normals.load_geometry(shaders_path("geometry_normals.glsl"));
	shader_normals.load_fragment(shaders_path("fragment_normals.glsl"));

	shader_gradient.load_vertex(shaders_path("vertex_gradient.glsl"));
	shader_gradient.load_geometry(shaders_path("geometry_gradient.glsl"));
	shader_gradient.load_fragment(shaders_path("fragment_gradient.glsl"));

	shader_pointcloud.load_vertex(shaders_path("vertex_pointcloud.glsl"));
	shader_pointcloud.load_fragment(shaders_path("fragment_pointcloud.glsl"));
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

void viewer::add_mesh(che * p_mesh)
{
	if(n_meshes == N_MESHES)
	{
		gproshan_log_var(n_meshes == N_MESHES);
		gproshan_log_var(n_meshes);
		return;
	}

	meshes[n_meshes].init(p_mesh);
	meshes[n_meshes].log_info();
	++n_meshes;

	idx_active_mesh = n_meshes - 1;
	glfwSetWindowTitle(window, active_mesh()->filename.c_str());

	const int & rows = m_window_size[n_meshes][0];
	const int & cols = m_window_size[n_meshes][1];
	for(index_t m = 0; m < n_meshes; ++m)
	{
		meshes[m].vx = m % cols;
		meshes[m].vy = rows - (m / cols) - 1;
	}

	glfwGetFramebufferSize(window, &viewport_width, &viewport_height);
	viewport_width /= m_window_size[n_meshes][1];
	viewport_height /= m_window_size[n_meshes][0];
}

void viewer::framebuffer_size_callback(GLFWwindow * window, int width, int height)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	view->viewport_width = width / m_window_size[view->n_meshes][1];
	view->viewport_height = height / m_window_size[view->n_meshes][0];
}

void viewer::window_size_callback(GLFWwindow * window, int width, int height)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	view->window_width = width;
	view->window_height = height;
}

void viewer::keyboard_callback(GLFWwindow * window, int key, int, int action, int)
{
	if(action == GLFW_RELEASE) return;

	if(key == GLFW_KEY_ESCAPE)
	{
		glfwSetWindowShouldClose(window, GLFW_TRUE);
		return;
	}

	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	if(ImGui::GetIO().WantCaptureKeyboard) return;

	process_t & pro = view->processes[key];
	if(pro.function)
	{
		pro.selected = !pro.selected;
		sprintf(view->status_message, "%s", pro.selected ? pro.name.c_str() : "");
	}
}

void viewer::mouse_callback(GLFWwindow * window, int button, int action, int mods)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	if(mods == GLFW_MOD_SHIFT && action == GLFW_RELEASE)
		view->pick_vertex(xpos, ypos);
	else if(button == GLFW_MOUSE_BUTTON_RIGHT)
	{
	}
	else
		view->cam.mouse(action == GLFW_PRESS, xpos, ypos, view->window_width, view->window_height);
}

void viewer::cursor_callback(GLFWwindow * window, double x, double y)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS)
	{
		viewer * view = (viewer *) glfwGetWindowUserPointer(window);
		if(ImGui::GetIO().WantCaptureMouse) return;

		view->cam.motion(x, y, view->window_width, view->window_height);
		view->action = true;
	}
}

void viewer::scroll_callback(GLFWwindow * window, double, double yoffset)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	if(ImGui::GetIO().WantCaptureMouse) return;

	if(yoffset > 0)
	{
		view->cam.zoom_in();
		view->action = true;
	}

	if(yoffset < 0)
	{
		view->cam.zoom_out();
		view->action = true;
	}
}

bool viewer::menu_help(viewer * view)
{
	for(auto & p: view->processes)
		if(p.second.function != nullptr)
			fprintf(stderr, "%16s: %s\n", ('[' + p.second.key + ']').c_str(), p.second.name.c_str());

	return false;
}

bool viewer::menu_save_load_view(viewer * view)
{
	filesystem::create_directory(tmp_file_path("views/"));

	static char file[128] = "new_view";

	ImGui::InputText("##savefile", file, sizeof(file));
	ImGui::SameLine();

	if(ImGui::Button("Save"))
	{
		ofstream os(tmp_file_path("views/" + file));
		os << view->cam;
		os.close();
	}

	static index_t select = 0;
	static vector<string> vfiles;

	vfiles.clear();
	for(auto & p: filesystem::directory_iterator(tmp_file_path("views/")))
		vfiles.push_back(p.path().string());

	if(!vfiles.size()) return true;

	if(ImGui::BeginCombo("##loadfile", vfiles[select].c_str()))
	{
		for(index_t i = 0; i < vfiles.size(); ++i)
		{
			if(ImGui::Selectable(vfiles[i].c_str(), select == i))
				select = i;

			if(select == i)
				ImGui::SetItemDefaultFocus();
		}

		ImGui::EndCombo();
	}

	ImGui::SameLine();

	if(ImGui::Button("Load"))
	{
		ifstream is(vfiles[select]);
		is >> view->cam;
		is.close();
	}

	return true;
}

bool viewer::menu_reset_mesh(viewer * view)
{
	view->active_mesh().selected.clear();
	view->other_vertices.clear();
	view->vectors.clear();

	view->active_mesh().reload();
	view->active_mesh().update_vbo();

	return false;
}

bool viewer::menu_save_mesh(viewer * view)
{
	const che * mesh = view->active_mesh();

	static char file[128] = "copy";
	static int format = 0;
	static int type_off = 0;
	static bool point_cloud = false;
	static bool vertex_color = false;

	ImGui::InputText("file", file, sizeof(file));
	ImGui::Combo("format", &format, ".off\0.obj\0.ply\0.xyz\0\0");

	switch(format)
	{
		case 0:
			ImGui::Combo("type off", &type_off, "OFF\0NOFF\0COFF\0NCOFF\0\0");
			ImGui::Checkbox("point cloud", &point_cloud);
			break;
		case 1:
			ImGui::Checkbox("vertex color", &vertex_color);
			ImGui::Checkbox("point cloud", &point_cloud);
			break;
		case 2:
			ImGui::Checkbox("vertex color", &vertex_color);
			break;
		case 3:
			ImGui::Checkbox("vertex color", &vertex_color);
			break;
	}

	if(ImGui::Button("Save"))
	{
		switch(format)
		{
			case 0: che_off::write_file(mesh, file, che_off::type(type_off), point_cloud);
				break;
			case 1: che_obj::write_file(mesh, file, vertex_color, point_cloud);
				break;
			case 2: che_ply::write_file(mesh, file, vertex_color);
				break;
			case 3: che_xyz::write_file(mesh, file, vertex_color);
				break;
		}

		sprintf(view->status_message, "file '%s' saved.", file);
	}

	return true;
}

bool viewer::menu_zoom_in(viewer * view)
{
	view->cam.zoom_in();

	return false;
}

bool viewer::menu_zoom_out(viewer * view)
{
	view->cam.zoom_out();

	return false;
}

bool viewer::menu_bgc_inc(viewer * view)
{
	if(view->bgc < 1) view->bgc += 0.05;
	else view->bgc = 1;

	glClearColor(view->bgc, view->bgc, view->bgc, 1.);

	return false;
}

bool viewer::menu_bgc_dec(viewer * view)
{
	if(view->bgc > 0) view->bgc -= 0.05;
	else view->bgc = 0;

	glClearColor(view->bgc, view->bgc, view->bgc, 1.);

	return false;
}

bool viewer::menu_bgc_white(viewer * view)
{
	view->bgc = 1;
	glClearColor(view->bgc, view->bgc, view->bgc, 1.);

	return false;
}

bool viewer::menu_bgc_black(viewer * view)
{
	view->bgc = 0;
	glClearColor(view->bgc, view->bgc, view->bgc, 1.);

	return false;
}

bool viewer::invert_orientation(viewer * view)
{
	view->active_mesh().invert_orientation();

	return false;
}

bool viewer::set_render_gl(viewer * view)
{
	view->render_opt = R_GL;

	return false;
}

#ifdef GPROSHAN_EMBREE
bool viewer::set_render_embree(viewer * view)
{
	static int rt_opt = 0;
	static double time = 0;

	ImGui::Combo("rt_opt", &rt_opt, "Mesh\0Splat\0\0");
	ImGui::InputFloat("pc_radius", &rt::embree::pc_radius, 0, 0, "%.3f");

	if(!view->rt_embree)
	{
		if(ImGui::Button("Start"))
		{
			view->render_opt = R_EMBREE;

			TIC(time);
			switch(rt_opt)
			{
				case 0: view->rt_embree = new rt::embree({view->active_mesh()});
						break;
				case 1: view->rt_embree = new rt::embree_splat({view->active_mesh()}, true);
						break;
			}
			TOC(time);

			if(!view->render_frame)
				view->render_frame = new frame;
		}
	}
	else
	{
		view->render_opt = R_EMBREE;

		if(ImGui::Button("Restart"))
		{
			delete view->rt_embree;
			view->rt_embree = nullptr;

			view->render_opt = R_GL;
		}
	}

	return true;
}
#endif // GPROSHAN_EMBREE

#ifdef GPROSHAN_OPTIX
bool viewer::set_render_optix(viewer * view)
{
	view->render_opt = 2;

	return false;
}
#endif // GPROSHAN_OPTIX

bool viewer::set_render_pointcloud(viewer * view)
{
	view->render_pointcloud = !view->render_pointcloud;

	return false;
}

bool viewer::set_render_wireframe(viewer * view)
{
	view->render_wireframe = !view->render_wireframe;

	return false;
}

bool viewer::set_render_wireframe_fill(viewer * view)
{
	view->render_wireframe_fill = !view->render_wireframe_fill;

	return false;
}

bool viewer::set_render_gradient_field(viewer * view)
{
	view->render_gradient_field = !view->render_gradient_field;

	return false;
}

bool viewer::set_render_normal_field(viewer * view)
{
	view->render_normal_field = !view->render_normal_field;

	return false;
}

bool viewer::set_render_border(viewer * view)
{
	view->render_border = !view->render_border;
	if(!view->render_border) view->active_mesh().selected.clear();

	return false;
}

bool viewer::set_render_lines(viewer * view)
{
	view->render_lines = !view->render_lines;

	return false;
}

bool viewer::set_render_flat(viewer * view)
{
	view->render_flat = !view->render_flat;
	view->action = true;

	return false;
}

#ifdef GPROSHAN_EMBREE
bool viewer::raycasting(viewer * view)
{
	rt::embree rc({view->active_mesh()});

	float * frame = rc.raycaster(	glm::uvec2(view->viewport_width, view->viewport_height),
									view->view_mat, view->proj_mat
									);

	std::thread([](CImg<float> img)
	{
		img.display();
	},
	CImg<float>(frame, view->viewport_width, view->viewport_height)).detach();

	delete [] frame;

	return false;
}
#endif // GPROSHAN_EMBREE

void viewer::render_gl()
{
	glProgramUniform3f(shader_sphere, shader_sphere("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);
	glProgramUniform3f(shader_sphere, shader_sphere("light"), light[0], light[1], light[2]);
	glProgramUniformMatrix4fv(shader_sphere, shader_sphere("model_view_mat"), 1, 0, &view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_sphere, shader_sphere("proj_mat"), 1, 0, &proj_mat[0][0]);
	glProgramUniform1f(shader_sphere, shader_sphere("scale"), cam.zoom);

	glProgramUniform1ui(shader_triangles, shader_triangles("idx_colormap"), idx_colormap);
	glProgramUniform3f(shader_triangles, shader_triangles("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);
	glProgramUniform3f(shader_triangles, shader_triangles("light"), light[0], light[1], light[2]);
	glProgramUniform1i(shader_triangles, shader_triangles("render_flat"), render_flat);
	glProgramUniform1i(shader_triangles, shader_triangles("render_lines"), render_lines);
	glProgramUniform1i(shader_triangles, shader_triangles("render_wireframe"), render_wireframe_fill);
	glProgramUniformMatrix4fv(shader_triangles, shader_triangles("model_view_mat"), 1, 0, &view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_triangles, shader_triangles("proj_mat"), 1, 0, &proj_mat[0][0]);

	glProgramUniform1ui(shader_pointcloud, shader_pointcloud("idx_colormap"), idx_colormap);
	glProgramUniform3f(shader_pointcloud, shader_pointcloud("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);
	glProgramUniform3f(shader_pointcloud, shader_pointcloud("light"), light[0], light[1], light[2]);
	glProgramUniformMatrix4fv(shader_pointcloud, shader_pointcloud("model_view_mat"), 1, 0, &view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_pointcloud, shader_pointcloud("proj_mat"), 1, 0, &proj_mat[0][0]);
	glProgramUniform1i(shader_pointcloud, shader_pointcloud("point_normals"), point_normals);
	glProgramUniform1ui(shader_pointcloud, shader_pointcloud("point_size"), point_size);


	draw_meshes(shader_triangles);


	if(render_normal_field)
	{
		glProgramUniform1f(shader_normals, shader_normals("length"), cam.zoom * 0.02);
		glProgramUniformMatrix4fv(shader_normals, shader_normals("model_view_mat"), 1, 0, &view_mat[0][0]);
		glProgramUniformMatrix4fv(shader_normals, shader_normals("proj_mat"), 1, 0, &proj_mat[0][0]);

		draw_meshes(shader_normals, true);
	}


	if(render_gradient_field)
	{
		glProgramUniform1f(shader_gradient, shader_gradient("length"), cam.zoom * 0.02);
		glProgramUniformMatrix4fv(shader_gradient, shader_gradient("model_view_mat"), 1, 0, &view_mat[0][0]);
		glProgramUniformMatrix4fv(shader_gradient, shader_gradient("proj_mat"), 1, 0, &proj_mat[0][0]);

		draw_meshes(shader_gradient);
	}


	if(render_border) select_border_vertices();

	draw_selected_vertices(shader_sphere);
}

#ifdef GPROSHAN_EMBREE
void viewer::render_embree()
{
	rt_embree->pathtracing(	glm::uvec2(viewport_width, viewport_height),
							view_mat, proj_mat, {glm_vec3(light)},
							render_flat, action
							);

	action = false;
	render_frame->display(viewport_width, viewport_height, rt_embree->img);
}
#endif // GPROSHAN_EMBREE

#ifdef GPROSHAN_OPTIX
void viewer::render_optix()
{
	if(!rt_optix)
	{
		double time_build_optix;
		TIC(time_build_optix);

			rt_optix = new rt::optix({active_mesh()});

		TOC(time_build_optix);
		gproshan_log_var(time_build_optix);
	}

	rt_optix->pathtracing(	glm::uvec2(viewport_width, viewport_height),
							view_mat, proj_mat, {glm_vec3(light)}, action);

	action = false;
	render_frame->display(viewport_width, viewport_height, rt_optix->img);
}
#endif // GPROSHAN_OPTIX

void viewer::draw_meshes(shader & program, const bool & normals)
{
	if(render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for(index_t i = 0; i < n_meshes; ++i)
	{
		glViewport(meshes[i].vx * viewport_width, meshes[i].vy * viewport_height, viewport_width, viewport_height);

		if(normals)
			meshes[i].draw_point_cloud(program);
		else if(meshes[i]->is_pointcloud() || render_pointcloud)
			meshes[i].draw_point_cloud(shader_pointcloud);
		else
			meshes[i].draw(program);
	}
}

void viewer::draw_selected_vertices(shader & program)
{
	if(sphere_translations.size() != active_mesh().selected.size())
	{
		sphere_translations.resize(active_mesh().selected.size());

		for(index_t i = 0; i < active_mesh().selected.size(); ++i)
			sphere_translations[i] = active_mesh()->gt(active_mesh().selected[i]);

		sphere.update_instances_translations(sphere_translations);
	}


	if(sphere_translations.size())
	{
		glViewport(active_mesh().vx * viewport_width, active_mesh().vy * viewport_height, viewport_width, viewport_height);
		sphere.draw(program);
	}
}

void viewer::select_border_vertices()
{
	che_viewer & mesh = active_mesh();

	mesh.selected.clear();

	vector<index_t> bounds = mesh->bounds();
	for(const index_t & b: bounds)
		for_boundary(he, mesh, b)
			mesh.selected.push_back(mesh->vt(he));
}

void viewer::pick_vertex(const real_t & x, const real_t & y)
{
	active_mesh().select(x, y, {viewport_width, viewport_height}, view_mat, proj_mat);
}


} // namespace gproshan

