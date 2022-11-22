#include <gproshan/viewer/viewer.h>

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <vector>
#include <thread>

#include <gproshan/mesh/che_off.h>
#include <gproshan/mesh/che_obj.h>
#include <gproshan/mesh/che_ply.h>
#include <gproshan/mesh/che_xyz.h>
#include <gproshan/mesh/che_pts.h>

#include <gproshan/raytracing/rt_embree.h>

#ifdef GPROSHAN_OPTIX
	#include <gproshan/raytracing/rt_optix.h>
#endif // GPROSHAN_OPTIX

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


const std::vector<ivec2> viewer::m_window_split = {	{1, 1},
													{1, 1}, {1, 2}, {1, 3},
													{2, 2}, {2, 3}, {2, 3},
													{2, 4}, {2, 4}, {2, 5},
													{2, 5}, {3, 4}, {3, 4},
													{4, 4}, {4, 4}, {4, 4}, {4, 4},
													{4, 5}, {4, 5}, {4, 5}, {4, 5}
													};

const size_t viewer::max_meshes = m_window_split.size() - 1;

const std::vector<std::string> viewer::colormap = { "vertex color",
													"blue",
													"red",
													"blue/read",
													"set"
													};

che_sphere viewer::sphere_data{0.01};

viewer::viewer(const int & width, const int & height)
{
	window_width = width;
	window_height = height;

	init_gl();
	init_glsl();
	init_imgui();
	init_menus();

	info_gl();
	gproshan_log_var(sizeof(real_t));

	sphere_data.update_normals();
	sphere = new che_viewer(&sphere_data);;

	frames = new frame[max_meshes];

	render_params.add_light({-1, 1, -2});
}

viewer::~viewer()
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();

	delete sphere;

	delete [] frames;

	for(che_viewer * m: meshes)
		delete m;
}

bool viewer::run()
{
	while(!glfwWindowShouldClose(window))
	{
		TIC(render_time)

		const quaternion & r = cam.current_rotation();

		cam_light = render_params.lights[0];
		cam_light = r.conj() * cam_light * r;

		proj_view_mat = proj_mat * cam.look_at(r);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		render_gl();

		TOC(render_time);

		imgui();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return true;
}

void viewer::imgui()
{
	if(hide_imgui) return;

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	che_viewer & mesh = active_mesh();

	if(ImGui::BeginMainMenuBar())
	{
		if(ImGui::BeginMenu("Select"))
		{
			for(index_t i = 0; i < meshes.size(); ++i)
			{
				const che_viewer & m = *meshes[i];
				if(ImGui::MenuItem((std::to_string(i) + ": " + m->filename).c_str(), nullptr, i == idx_active_mesh, i != idx_active_mesh))
				{
					idx_active_mesh = i;
					glfwSetWindowTitle(window, m->filename.c_str());
				}
			}

			ImGui::EndMenu();
		}

		if(ImGui::BeginMenu("Color"))
		{
			for(index_t i = 0; i < colormap.size(); ++i)
			{
				if(ImGui::MenuItem(colormap[i].c_str(), nullptr, i == mesh.idx_colormap, i != mesh.idx_colormap))
					check_apply_all_meshes([&](che_viewer & mesh)
					{
						mesh.idx_colormap = i;
					});
				ImGui::Separator();
			}
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

					//ImGui::Separator();
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

	ImGui::SetNextWindowSize(ImVec2(360, -1));
	ImGui::SetNextWindowPos(ImVec2(20, 60), ImGuiCond_Once);
	ImGui::Begin("gproshan");

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5);

	ImGui::Checkbox("apply options to all meshes\nmenus: [color, render, mesh]", &apply_all_meshes);
	if(ImGui::CollapsingHeader(mesh->filename.c_str(), ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Text("%13lu fps", size_t(1.0 / render_time));
		ImGui::Text("%13lu vertices", mesh->n_vertices);
		ImGui::Text("%13lu faces", mesh->n_faces);

		ImGui::Indent();
		if(ImGui::Combo("fit screen", (int *) &mesh.opt_fit_screen, "none\0box (2x2x2)\0sphere (97.72%)\0\0"))
		{
			mesh.update_model_mat();
		}
		if(mesh.render_pointcloud)
		{
			ImGui::Checkbox("point_normals", &mesh.point_normals);
			ImGui::SliderInt("point_size", (int *) &mesh.point_size, 1, 32);
		}
		ImGui::Unindent();
	}

	static real_t pos_min = -10;
	static real_t pos_max = 10;
	if(ImGui::CollapsingHeader("Camera"))
	{
		ImGui::Indent();
		ImGui::SliderScalarN("position", ImGuiDataType_Real, &cam.pos[0], 3, &pos_min, &pos_max);
		ImGui::Unindent();
	}

	static char slight[32];
	if(ImGui::CollapsingHeader("Scene Lights"))
	{
		ImGui::Indent();

		for(int i = 0; i < render_params.n_lights; ++i)
		{
			sprintf(slight, "light %d", i);
			ImGui::SliderScalarN(slight, ImGuiDataType_Real, &render_params.lights[i], 3, &pos_min, &pos_max);
		}

		if(ImGui::Button("add light"))
			render_params.add_light({0, 0, 0});

		if(render_params.n_lights > 1)
		{
			ImGui::SameLine();
			if(ImGui::Button("remove light"))
				--render_params.n_lights;
		}

		ImGui::SameLine();
		if(ImGui::Button("show lights"))
		{
			sphere_points.clear();
			for(int i = 0; i < render_params.n_lights; ++i)
				sphere_points.push_back(render_params.lights[i]);
		}

		if(ImGui::Button("add selected points as lights"))
		{
			for(const index_t & v: mesh.selected)
				if(!render_params.add_light(mesh.model_mat * vec4(mesh->point(v), 1)))
					break;
		}

		ImGui::Unindent();
	}


	for(auto & p: processes)
	{
		process_t & pro = p.second;
		if(ImGui::CollapsingHeader(("[" + pro.key + "] " + pro.name).c_str(), &pro.selected, ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::PushID(pro.name.c_str());
			ImGui::Indent();
			pro.selected = pro.selected && pro.function(this);
			ImGui::Unindent();
			ImGui::PopID();
		}
	}

	ImGui::PopItemWidth();
	ImGui::End();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

che_viewer & viewer::active_mesh()
{
	return *meshes[idx_active_mesh];
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
	add_process(GLFW_KEY_F1, "F1", "Help", m_help);
	add_process(GLFW_KEY_ESCAPE, "ESCAPE", "Close", m_close);
	add_process(GLFW_KEY_I, "F1", "Hide/Show ImGui", m_hide_show_imgui);
	add_process(GLFW_KEY_PERIOD, "PERIOD", "Save/Load view", m_save_load_view);
	add_process(GLFW_KEY_UP, "UP", "Zoom in", m_zoom_in);
	add_process(GLFW_KEY_DOWN, "DOWN", "Zoom out", m_zoom_out);
	add_process(GLFW_KEY_RIGHT, "RIGHT", "Background color inc", m_bgc_inc);
	add_process(GLFW_KEY_LEFT, "LEFT", "Background color dec", m_bgc_dec);
	add_process(GLFW_KEY_1, "1", "Background color white", m_bgc_white);
	add_process(GLFW_KEY_0, "0", "Background color black", m_bgc_black);

	sub_menus.push_back("Render");
	add_process(GLFW_KEY_F5, "F5", "Render Point Cloud", m_render_pointcloud);
	add_process(GLFW_KEY_F6, "F6", "Render Wireframe", m_render_wireframe);
	add_process(GLFW_KEY_F7, "F7", "Render Triangles", m_render_triangles);
	add_process(GLFW_KEY_F8, "F8", "Render GL", m_render_gl);
	add_process(GLFW_KEY_SPACE, "SPACE", "Level Curves", m_render_lines);
	add_process(GLFW_KEY_TAB, "TAB", "Render Flat", m_render_flat);
	add_process(GLFW_KEY_R, "R", "Setup Raytracing", m_setup_raytracing);
	add_process(GLFW_KEY_F9, "F9", "Render Embree", m_render_embree);
	add_process(GLFW_KEY_ENTER, "ENTER", "Raycasting", m_raycasting);

#ifdef GPROSHAN_OPTIX
	add_process(GLFW_KEY_F10, "F10", "Render OptiX", m_render_optix);
#endif // GPROSHAN_OPTIX

	sub_menus.push_back("Mesh");
	add_process(GLFW_KEY_BACKSPACE, "BACKSPACE", "Reload/Reset", m_reset_mesh);
	add_process(GLFW_KEY_W, "W", "Save Mesh", m_save_mesh);
	add_process(0, "", "Normalize Mesh", m_normalize_mesh);
	add_process(GLFW_KEY_F2, "F2", "Invert Normals", m_invert_normals);
	add_process(GLFW_KEY_F3, "F3", "Gradient Field", m_render_gradients);
	add_process(GLFW_KEY_F4, "F4", "Normal Field", m_render_normals);
	add_process(GLFW_KEY_B, "B", "Select Border Vertices", m_select_border_vertices);
	add_process(GLFW_KEY_C, "C", "Clean Selected Vertices", m_clean_selected_vertices);
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

void viewer::add_process(const int & key, const std::string & skey, const std::string & name, const function_t & f)
{
	if(processes.find(key) == processes.end())
	{
		processes[key] = {skey, name, f};
		processes[key].sub_menu = sub_menus.size() - 1;
	}
	else std::cerr << "Repeat key: " << key << std::endl;
}

bool viewer::add_mesh(che * p_mesh, const bool & reset_normals)
{
	if(meshes.size() == max_meshes)
		return false;

	if(reset_normals)
		p_mesh->update_normals();

	meshes.push_back(new che_viewer(p_mesh));
	che_viewer & mesh = *meshes.back();
	mesh.log_info();

	idx_active_mesh = meshes.size() - 1;
	glfwSetWindowTitle(window, mesh->filename.c_str());

	const int & rows = m_window_split[meshes.size()].x();
	const int & cols = m_window_split[meshes.size()].y();
	for(index_t m = 0; m < meshes.size(); ++m)
	{
		meshes[m]->vx = m % cols;
		meshes[m]->vy = rows - (m / cols) - 1;
	}

	glfwGetFramebufferSize(window, &viewport_width, &viewport_height);
	viewport_width /= cols;
	viewport_height /= rows;
	cam.aspect = real_t(viewport_width) / viewport_height;
	proj_mat = cam.perspective();


	return true;
}

void viewer::framebuffer_size_callback(GLFWwindow * window, int width, int height)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	view->viewport_width = width / m_window_split[view->meshes.size()].y();
	view->viewport_height = height / m_window_split[view->meshes.size()].x();
	view->cam.aspect = real_t(view->viewport_width) / view->viewport_height;
	view->proj_mat = view->cam.perspective();
}

void viewer::window_size_callback(GLFWwindow * window, int width, int height)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	view->window_width = width;
	view->window_height = height;
}

void viewer::keyboard_callback(GLFWwindow * window, int key, int, int action, int)
{
	if(ImGui::GetIO().WantCaptureKeyboard) return;

	if(action == GLFW_RELEASE) return;

	viewer * view = (viewer *) glfwGetWindowUserPointer(window);

	process_t & pro = view->processes[key];
	if(pro.function)
	{
		pro.selected = view->hide_imgui ? pro.function(view) && pro.selected : !pro.selected;
		sprintf(view->status_message, "%s", pro.selected ? pro.name.c_str() : "");
	}

}

void viewer::mouse_callback(GLFWwindow * window, int button, int action, int mods)
{
	if(ImGui::GetIO().WantCaptureMouse) return;

	viewer * view = (viewer *) glfwGetWindowUserPointer(window);

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
	{
		float xscale, yscale;
		glfwGetWindowContentScale(window, &xscale, &yscale);

		const index_t & ix = xpos * xscale;
		const index_t & iy = ypos * yscale;
		const int & cols = m_window_split[view->meshes.size()].y();
		const index_t & idx_mesh = cols * (iy / view->viewport_height) + ix / view->viewport_width;
		if(idx_mesh < view->meshes.size())
			view->idx_active_mesh = idx_mesh;

		if(mods == GLFW_MOD_SHIFT)
			view->pick_vertex(ix % view->viewport_width, iy % view->viewport_height);
	}

	if(button == GLFW_MOUSE_BUTTON_LEFT)
		view->cam.mouse(action == GLFW_PRESS, xpos, ypos, view->window_width, view->window_height);
}

void viewer::cursor_callback(GLFWwindow * window, double x, double y)
{
	if(ImGui::GetIO().WantCaptureMouse) return;

	viewer * view = (viewer *) glfwGetWindowUserPointer(window);

	if(GLFW_PRESS == glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT))
	{
		view->cam.motion(x, y, view->window_width, view->window_height);
		view->render_params.restart = true;
	}

	if(GLFW_PRESS == glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT))
	{
		view->cam.pos.im().x() = 2 * x / view->window_width - 1;
		view->cam.pos.im().y() = 2 * y / view->window_height - 1;
		view->render_params.restart = true;
	}
}

void viewer::scroll_callback(GLFWwindow * window, double, double yoffset)
{
	viewer * view = (viewer *) glfwGetWindowUserPointer(window);
	if(ImGui::GetIO().WantCaptureMouse) return;

	if(yoffset > 0)
	{
		view->cam.zoom_in();
		view->render_params.restart = true;
	}

	if(yoffset < 0)
	{
		view->cam.zoom_out();
		view->render_params.restart = true;
	}
}

bool viewer::m_help(viewer * view)
{
	for(auto & p: view->processes)
		if(p.second.function != nullptr)
			fprintf(stderr, "%16s: %s\n", ('[' + p.second.key + ']').c_str(), p.second.name.c_str());

	return false;
}

bool viewer::m_close(viewer * view)
{
	glfwSetWindowShouldClose(view->window, GLFW_TRUE);
	return false;
}

bool viewer::m_hide_show_imgui(viewer * view)
{
	view->hide_imgui = !view->hide_imgui;
	return false;
}

bool viewer::m_save_load_view(viewer * view)
{
	std::filesystem::create_directory(tmp_file_path("views/"));

	static char file[128] = "new_view";

	ImGui::InputText("##savefile", file, sizeof(file));
	ImGui::SameLine();

	if(ImGui::Button("Save"))
	{
		std::ofstream os(tmp_file_path(std::string("views/") + file));
		os << view->cam;
		os.close();
	}

	static index_t select = 0;
	static std::vector<std::string> vfiles;

	vfiles.clear();
	for(auto & p: std::filesystem::directory_iterator(tmp_file_path("views/")))
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
		std::ifstream is(vfiles[select]);
		is >> view->cam;
		is.close();
	}

	return true;
}

bool viewer::m_reset_mesh(viewer * view)
{
	view->vectors.clear();

	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.selected.clear();
		mesh->reload();
		mesh->update_normals();
		mesh.update();
	});

	return false;
}

bool viewer::m_save_mesh(viewer * view)
{
	const che * mesh = view->active_mesh();

	static char file[128] = "copy";
	static int format = 0;
	static int type_off = 0;
	static bool point_cloud = false;
	static bool vertex_color = false;

	ImGui::InputText("file", file, sizeof(file));
	ImGui::Combo("format", &format, ".off\0.obj\0.ply\0.xyz\0.pts\0");

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
			case 4: che_pts::write_file(mesh, file);
				break;
		}

		sprintf(view->status_message, "file '%s' saved.", file);
	}

	return true;
}

bool viewer::m_normalize_mesh(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh->normalize_sphere();
		mesh.update();
	});

	return false;
}

bool viewer::m_zoom_in(viewer * view)
{
	view->cam.zoom_in();
	return false;
}

bool viewer::m_zoom_out(viewer * view)
{
	view->cam.zoom_out();
	return false;
}

bool viewer::m_bgc_inc(viewer * view)
{
	if(view->bgc < 1) view->bgc += 0.05;
	else view->bgc = 1;

	glClearColor(view->bgc, view->bgc, view->bgc, 1);

	return false;
}

bool viewer::m_bgc_dec(viewer * view)
{
	if(view->bgc > 0) view->bgc -= 0.05;
	else view->bgc = 0;

	glClearColor(view->bgc, view->bgc, view->bgc, 1);

	return false;
}

bool viewer::m_bgc_white(viewer * view)
{
	view->bgc = 1;
	glClearColor(view->bgc, view->bgc, view->bgc, 1);

	return false;
}

bool viewer::m_bgc_black(viewer * view)
{
	view->bgc = 0;
	glClearColor(view->bgc, view->bgc, view->bgc, 1);

	return false;
}

bool viewer::m_setup_raytracing(viewer * view)
{
	che_viewer & mesh = view->active_mesh();

	static int rt = 0;
	static double time = 0;
	static float pc_radius = 0.01;

	ImGui::Combo("rt", &rt, "Select\0Embree\0OptiX\0\0");
	ImGui::InputFloat("pc_radius (if render_pointcloud)", &pc_radius, 0, 0, "%.3f");

	if(ImGui::Button("Build"))
	{
		switch(rt)
		{
			case R_GL: break;

			case R_EMBREE:
				delete mesh.rt_embree;
				TIC(time);
					mesh.rt_embree = new rt::embree({mesh}, {mesh.model_mat}, mesh.render_pointcloud, pc_radius);
				TOC(time);
				sprintf(view->status_message, "build embree in %.3fs", time);
				break;

			case R_OPTIX:
			#ifdef GPROSHAN_OPTIX
				delete mesh.rt_optix;
				TIC(time);
					mesh.rt_optix = new rt::optix({mesh}, {mesh.model_mat});
				TOC(time);
				sprintf(view->status_message, "build optix in %.3fs", time);
			#endif // GPROSHAN_OPTIX
				break;
		}
	}

	return true;
}

bool viewer::m_render_gl(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_opt = R_GL;
	});

	return false;
}

bool viewer::m_render_embree(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_opt = R_EMBREE;
	});
	view->render_params.restart = true;

	return false;
}

bool viewer::m_render_optix(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_opt = R_OPTIX;
	});
	view->render_params.restart = true;

	return false;
}

bool viewer::m_invert_normals(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh->invert_normals();
		mesh.update_vbo_normal();
	});

	return false;
}

bool viewer::m_select_border_vertices(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		for(const index_t & b: mesh->bounds())
			for(const index_t & v: mesh->boundary(b))
				mesh.selected.push_back(v);
	});

	return false;
}

bool viewer::m_clean_selected_vertices(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.selected.clear();
	});

	return false;
}

bool viewer::m_render_pointcloud(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_pointcloud = !mesh.render_pointcloud;
	});

	return false;
}

bool viewer::m_render_wireframe(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_wireframe = !mesh.render_wireframe;
	});

	return false;
}

bool viewer::m_render_triangles(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_triangles = !mesh.render_triangles;
	});

	return false;
}

bool viewer::m_render_gradients(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_gradients = !mesh.render_gradients;
	});

	return false;
}

bool viewer::m_render_normals(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_normals = !mesh.render_normals;
	});

	return false;
}

bool viewer::m_render_lines(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_lines = !mesh.render_lines;
	});

	return false;
}

bool viewer::m_render_flat(viewer * view)
{
	view->check_apply_all_meshes([&](che_viewer & mesh)
	{
		mesh.render_flat = !mesh.render_flat;
		view->render_params.restart = true;
	});

	return false;
}

bool viewer::m_raycasting(viewer * view)
{
	che_viewer & mesh = view->active_mesh();

	rt::embree rc({mesh}, {mesh.model_mat});

	float * frame = rc.raycaster(	{view->viewport_width, view->viewport_height},
									inverse(view->proj_view_mat),
									view->cam.eye
									);

	std::thread([](CImg<float> img)
	{
		img.display();
	},
	CImg<float>(frame, view->viewport_width, view->viewport_height)).detach();

	delete [] frame;

	return false;
}

void viewer::render_gl()
{
	glProgramUniform3f(shader_sphere, shader_sphere("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);
	glProgramUniform3f(shader_triangles, shader_triangles("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);
	glProgramUniform3f(shader_pointcloud, shader_pointcloud("eye"), cam.eye[0], cam.eye[1], cam.eye[2]);

	glProgramUniform3f(shader_sphere, shader_sphere("cam_light"), cam_light[0], cam_light[1], cam_light[2]);
	glProgramUniform3f(shader_triangles, shader_triangles("cam_light"), cam_light[0], cam_light[1], cam_light[2]);
	glProgramUniform3f(shader_pointcloud, shader_pointcloud("cam_light"), cam_light[0], cam_light[1], cam_light[2]);

	glProgramUniformMatrix4fv(shader_sphere, shader_sphere("proj_view_mat"), 1, true, &proj_view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_triangles, shader_triangles("proj_view_mat"), 1, true, &proj_view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_pointcloud, shader_pointcloud("proj_view_mat"), 1, true, &proj_view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_normals, shader_normals("proj_view_mat"), 1, true, &proj_view_mat[0][0]);
	glProgramUniformMatrix4fv(shader_gradient, shader_gradient("proj_view_mat"), 1, true, &proj_view_mat[0][0]);

	glProgramUniform1f(shader_normals, shader_normals("length"), cam.zoom() * 0.02);
	glProgramUniform1f(shader_gradient, shader_gradient("length"), cam.zoom() * 0.02);

	glProgramUniform1f(shader_sphere, shader_sphere("scale"), cam.zoom());


	for(index_t i = 0; i < meshes.size(); ++i)
	{
		che_viewer & mesh = *meshes[i];

		glViewport(mesh.vx * viewport_width, mesh.vy * viewport_height, viewport_width, viewport_height);

		if(mesh.render_opt != R_GL)
		{
			render_rt(mesh, frames[i]);
			continue;
		}

		if(mesh->is_pointcloud() || mesh.render_pointcloud)
			mesh.draw_point_cloud(shader_pointcloud);
		else
			mesh.draw(shader_triangles);

		if(mesh.render_normals)
			mesh.draw_point_cloud(shader_normals);

		if(mesh.render_gradients)
			mesh.draw(shader_gradient);

		mesh.draw_selected_vertices(*sphere, shader_sphere);

		if(sphere_points.size())
		{
			sphere->model_mat = mat4::identity();
			sphere->update_instances_positions(sphere_points);
			sphere->draw(shader_sphere);
		}
	}

	render_params.restart = false;
}

void viewer::render_rt(che_viewer & mesh, frame & rt_frame)
{
	rt::raytracing * rt = nullptr;
	if(mesh.render_opt == R_EMBREE) rt = mesh.rt_embree;
	if(mesh.render_opt == R_OPTIX) rt = mesh.rt_optix;

	if(!rt) return;

	render_params.restart = rt_frame.resize(viewport_width, viewport_height) || render_params.restart;

	//render_params.viewport_x = mesh.vx * viewport_width;
	//render_params.viewport_y = mesh.vy * viewport_height;
	//render_params.viewport_is_window = false;
	render_params.inv_proj_view = inverse(proj_view_mat);
	render_params.cam_pos = cam.eye;

	rt->render(rt_frame.map_pbo(mesh.render_opt == R_OPTIX), render_params, mesh.render_flat);
	rt_frame.unmap_pbo(mesh.render_opt == R_OPTIX);

	rt_frame.display();
}

void viewer::pick_vertex(const int & x, const int & y)
{
	che_viewer & mesh = active_mesh();

	mesh.select({x, y}, {viewport_width, viewport_height}, inverse(proj_view_mat), cam.eye);
}

void viewer::check_apply_all_meshes(const std::function<void(che_viewer &)> & fun)
{
	if(!apply_all_meshes)
	{
		fun(active_mesh());
		return;
	}

	for(auto & m: meshes)
		fun(*m);
}


} // namespace gproshan

