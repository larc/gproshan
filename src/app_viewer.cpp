#include "app_viewer.h"

#include <random>


using namespace std;
using namespace gproshan::mdict;


// geometry processing and shape analysis framework
namespace gproshan {


app_viewer::~app_viewer()
{
	for(che * mesh: meshes)
		delete mesh;
}

che * app_viewer::load_mesh(const string & file_path)
{
	size_t pos = file_path.rfind('.');
	assert(pos != string::npos);

	string extension = file_path.substr(pos + 1);

	if(extension == "off") return new che_off(file_path);
	if(extension == "obj") return new che_obj(file_path);
	if(extension == "ply") return new che_ply(file_path);
	if(extension == "ptx") return new che_ptx(file_path);

	return new che_img(file_path);
}

int app_viewer::main(int nargs, const char ** args)
{
	if(nargs < 2)
	{
		printf("%s [mesh_paths.(off,obj,ply)]\n", args[0]);
		return 0;
	}

	TIC(time)

	for(int i = 1; i < nargs; ++i)
		add_mesh(load_mesh(args[i]));

	TOC(time)

	gproshan_log_var(sizeof(real_t));
	gproshan_log_var(time);

	init();
	run();

	return 0;
}

void app_viewer::init()
{
	sub_menus.push_back("Geometry");
	add_process(GLFW_KEY_H, {"H", "Convex Hull", process_convex_hull});
	add_process(GLFW_KEY_K, {"K", "Gaussian curvature", process_gaussian_curvature});
	add_process(GLFW_KEY_8, {"8", "Edge Collapse", process_edge_collapse});
	add_process(GLFW_KEY_9, {"9", "Multiplicate", process_multiplicate_vertices});
	add_process(GLFW_KEY_DELETE, {"DELETE", "Delete vertices", process_delete_vertices});
	add_process(GLFW_KEY_MINUS, {"MINUS", "Delete non-manifold vertices", process_delete_non_manifold_vertices});

	sub_menus.push_back("Fairing");
	add_process(GLFW_KEY_T, {"T", "Fairing Taubin", process_fairing_taubin});
	add_process(GLFW_KEY_E, {"E", "Fairing Spectral", process_fairing_spectral});

	sub_menus.push_back("Geodesics");
	add_process(GLFW_KEY_F, {"F", "Fast Marching", process_geodesics_fm});
	add_process(GLFW_KEY_C, {"C", "Parallel Toplesets Propagation CPU", process_geodesics_ptp_cpu});
	add_process(GLFW_KEY_M, {"M", "Heat Method", process_geodesics_heat_method});
#ifdef GPROSHAN_CUDA
	add_process(GLFW_KEY_G, {"G", "Parallel Toplesets Propagation GPU", process_geodesics_ptp_gpu});
	add_process(GLFW_KEY_Q, {"Q", "Heat Method GPU", process_geodesics_heat_method_gpu});
#endif // GPROSHAN_CUDA
	add_process(GLFW_KEY_S, {"S", "Geodesic Farthest Point Sampling", process_farthest_point_sampling});
	add_process(GLFW_KEY_R, {"R", "Geodesic Farthest Point Sampling (radio)", process_farthest_point_sampling_radio});
	add_process(GLFW_KEY_V, {"V", "Geodesic Voronoi", process_voronoi});
	add_process(GLFW_KEY_P, {"P", "Toplesets", process_compute_toplesets});

	sub_menus.push_back("Sparse Coding");
	add_process(GLFW_KEY_I, {"I", "Mesh Sparse Coding", process_msparse_coding});
	add_process(GLFW_KEY_J, {"J", "MDICT Patch", process_mdict_patch});
	add_process(GLFW_KEY_D, {"D", "MDICT Mask", process_mask});
	add_process(GLFW_KEY_L, {"L", "PC reconstruction", process_pc_reconstruction});

	sub_menus.push_back("Features");
	add_process(GLFW_KEY_2, {"2", "Functional Maps", process_functional_maps});
	add_process(GLFW_KEY_3, {"3", "GPS", process_gps});
	add_process(GLFW_KEY_4, {"4", "HKS", process_hks});
	add_process(GLFW_KEY_5, {"5", "WKS", process_wks});
	add_process(GLFW_KEY_6, {"6", "Key Points", process_key_points});
	add_process(GLFW_KEY_7, {"7", "Key Components", process_key_components});

	sub_menus.push_back("Hole Filling");
	add_process(GLFW_KEY_X, {"X", "Poisson Membrane surface", process_poisson_laplacian_1});
	add_process(GLFW_KEY_Y, {"Y", "Poisson Thin-plate surface", process_poisson_laplacian_2});
	add_process(GLFW_KEY_Z, {"Z", "Poisson Minimum variation surface", process_poisson_laplacian_3});
	add_process(GLFW_KEY_A, {"A", "Fill hole - planar mesh", process_fill_holes});
	add_process(GLFW_KEY_B, {"B", "Fill hole - biharmonic splines", process_fill_holes_biharmonic_splines});

	sub_menus.push_back("Others");
	add_process(GLFW_KEY_SEMICOLON, {"SEMICOLON", "Select multiple vertices", process_select_multiple});
	add_process(GLFW_KEY_SLASH, {"SLASH", "Threshold", process_threshold});
	add_process(GLFW_KEY_N, {"N", "Noise", process_noise});
	add_process(GLFW_KEY_COMMA, {"COMMA", "Black noise", process_black_noise});
}


// Geometry

bool app_viewer::process_convex_hull(viewer * p_view)
{
	return true;
}

bool app_viewer::process_gaussian_curvature(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	real_t g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

	a_vec gv(mesh->n_vertices);

	#pragma omp parallel for private(g, a, b) reduction(max: g_max) reduction(min: g_min)
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		g = 0;
		for_star(he, mesh, v)
		{
			a = mesh->gt_vt(next(he)) - mesh->gt(v);
			b = mesh->gt_vt(prev(he)) - mesh->gt(v);
			g += acos((a,b) / (*a * *b));
		}
		//gv(v) = (2 * M_PI - g) / mesh->area_vertex(v);
		gv(v) = mesh->mean_curvature(v);

		g_max = max(g_max, gv(v));
		g_min = min(g_min, gv(v));
	}

	g = g_max - g_min;
	gproshan_log_var(g);
	gproshan_log_var(g_min);
	gproshan_log_var(g_max);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		gv(v) = (gv(v) + g_min) / g;

	real_t gm = mean(gv);
	real_t gs = var(gv);

	gproshan_debug_var(gm);
	gproshan_debug_var(gs);

	auto f = [&](real_t x, real_t a = 4) -> real_t
	{
		if(x < gm - a * gs) return 0;
		if(x > gm + a * gs) return 1;
		return (x - gm) / (2 * a * gs) + 0.5;
	};

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		mesh->heatmap(v) = f(gv(v));

	return false;
}

bool app_viewer::process_edge_collapse(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	index_t levels;
	cin >> levels;

	TIC(view->time) simplification sampling(mesh, &mesh->normal(0), levels); TOC(view->time)
	gproshan_debug_var(view->time);

	//if(view->n_meshes < 2)
	//	view->add_mesh(new che(*mesh));

	return false;
}

bool app_viewer::process_multiplicate_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	mesh->multiplicate_vertices();
	mesh.update();

	mesh.log_info();

	return false;
}

bool app_viewer::process_delete_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	if(!mesh.selected.size()) return true;

	gproshan_debug(removing vertex);
	mesh->remove_vertices(mesh.selected);
	mesh.selected.clear();
	gproshan_debug(removing vertex);

	return false;
}

bool app_viewer::process_delete_non_manifold_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	gproshan_debug(removing vertex);
	mesh->remove_non_manifold_vertices();
	gproshan_debug(removing vertex);

	return false;
}


// Fairing

bool app_viewer::process_fairing_spectral(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static int k = 100;
	ImGui::SliderInt("eigenvectors", &k, 1, mesh->n_vertices / 6);

	if(ImGui::Button("Run"))
	{
		fairing_spectral fair(k);
		fair.run(mesh);

		mesh->set_vertices(fair.get_postions());
		mesh->update_normals();
	}

	return true;
}

bool app_viewer::process_fairing_taubin(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static float step = 0.001; //cin >> step;
	ImGui::InputFloat("step", &step, 0.001, 1, "%.3f");

	if(ImGui::Button("Run"))
	{
		fairing_taubin fair(step);
		fair.run(mesh);

		mesh->set_vertices(fair.get_postions());
		mesh->update_normals();
	}

	return true;
}


// Geodesics

bool app_viewer::process_geodesics(viewer * p_view, const geodesics::algorithm & alg)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	if(!mesh.selected.size())
		mesh.selected.push_back(0);


	static vector<real_t> dist;

	if(dist.size() != mesh->n_vertices)
		dist.resize(mesh->n_vertices);

	geodesics::params params;
	params.alg			= alg;
	params.dist_alloc	= dist.data();

	TIC(view->time)
		geodesics G(mesh, mesh.selected, params);
	TOC(view->time)
	gproshan_log_var(view->time);

	mesh->update_heatmap(&G[0]);

	return false;
}

bool app_viewer::process_geodesics_fm(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_geodesics(p_view, geodesics::FM);
}

bool app_viewer::process_geodesics_ptp_cpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_geodesics(p_view, geodesics::PTP_CPU);
}

bool app_viewer::process_geodesics_heat_method(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_geodesics(p_view, geodesics::HEAT_METHOD);
}

#ifdef GPROSHAN_CUDA

bool app_viewer::process_geodesics_ptp_gpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_geodesics(p_view, geodesics::PTP_GPU);
}

bool app_viewer::process_geodesics_heat_method_gpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_geodesics(p_view, geodesics::HEAT_METHOD_GPU);
}

#endif // GPROSHAN_CUDA

bool app_viewer::process_farthest_point_sampling(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static int n = 10;
	static real_t radio;

	ImGui::SliderInt("samples", &n, 1, mesh->n_vertices / 6);
	ImGui::Text("radio: %.3f", radio);

	if(ImGui::Button("Run"))
	{
		TIC(view->time)
		load_sampling(mesh.selected, radio, mesh, n);
		TOC(view->time)
		gproshan_log_var(view->time);
	}

	return true;
}

bool app_viewer::process_farthest_point_sampling_radio(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	gproshan_input(radio);
	real_t radio; cin >> radio;

#ifdef GPROSHAN_CUDA	// IMPLEMENT/REVIEW
	double time_fps;

	TIC(view->time)
	radio = farthest_point_sampling_ptp_gpu(mesh, mesh.selected, time_fps, NIL, radio);
	TOC(view->time)
	gproshan_log_var(time_fps);
#endif // GPROSHAN_CUDA

	gproshan_log_var(radio);
	gproshan_log_var(mesh.selected.size());
	gproshan_log_var(view->time);

	return false;
}

bool app_viewer::process_voronoi(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	geodesics::params params;
	params.cluster = true;

#ifdef GPROSHAN_CUDA
	params.alg = geodesics::PTP_GPU;
#endif

	TIC(view->time)
	geodesics ptp(mesh, mesh.selected, params);
	TOC(view->time)

	gproshan_log_var(view->time);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		mesh->heatmap(i) = ptp.clusters[i];
		mesh->heatmap(i) /= mesh.selected.size() + 1;
	}

	return false;
}

bool app_viewer::process_compute_toplesets(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	if(!mesh.selected.size())
		mesh.selected.push_back(0);

	index_t * toplesets = new index_t[mesh->n_vertices];
	index_t * sorted = new index_t[mesh->n_vertices];
	vector<index_t> limites;
	mesh->compute_toplesets(toplesets, sorted, limites, mesh.selected);

	size_t n_toplesets = limites.size() - 1;

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(toplesets[v] < n_toplesets)
			mesh->heatmap(v) = real_t(toplesets[v]) / (n_toplesets);
	}

	gproshan_debug_var(n_toplesets);

	delete [] toplesets;
	delete [] sorted;

	return false;
}


// Mesh Sparse Coding

bool app_viewer::process_msparse_coding(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static msparse_coding::params params;
	static size_t n = 12;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui_InputReal("nyquist_factor", &patch::nyquist_factor, 0.01, 0.01, "%.2lf");
	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui_InputReal("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui_InputReal("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui_InputReal("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		real_t max_error = msc.execute();
		gproshan_log_var(max_error);

		mesh->update_heatmap(&msc[0]);
		mesh->update_normals();
	}

	return true;
}

bool app_viewer::process_mdict_patch(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	TIC(view->time)
	index_t * toplevel = new index_t[mesh->n_vertices];
	size_t avg_nvp = 0;

	vertex vdir;
	patch p;
	real_t mean_edge = mesh->mean_edge();
	for(auto & v: mesh.selected)
	{
		p.init(mesh, v, msparse_coding::T, msparse_coding::T * mean_edge, toplevel);
		for(auto & u: p.vertices)
			mesh->heatmap(u) = 1;

		vdir.x = p.T(0, 0);
		vdir.y = p.T(0, 1);
		vdir.z = p.T(0, 2);
		view->vectors.push_back(mesh->gt(v));
		view->vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);

		vdir.x = p.T(1, 0);
		vdir.y = p.T(1, 1);
		vdir.z = p.T(1, 2);
		view->vectors.push_back(mesh->gt(v));
		view->vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);

		vdir.x = p.T(2, 0);
		vdir.y = p.T(2, 1);
		vdir.z = p.T(2, 2);
		view->vectors.push_back(mesh->gt(v));
		view->vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);

		view->vectors.push_back(mesh->gt(v));
		view->vectors.push_back(mesh->gt(v) + 3 * mean_edge * mesh->normal(v));

		avg_nvp += p.vertices.size();
	}

	avg_nvp /= mesh.selected.size();
	gproshan_debug_var(avg_nvp);

	delete [] toplevel;
	TOC(view->time)
	gproshan_debug_var(view->time);

	return false;
}

bool app_viewer::process_mask(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static msparse_coding::params params;
	static size_t n = 12;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui_InputReal("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui_InputReal("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui_InputReal("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		msc.init_radial_feature_patches();
		//dict.init_voronoi_patches();
		mesh->update_heatmap(&msc[0]);
		string f_points = tmp_file_path(string(msc) + ".rsampl");

		a_vec points_out;
		gproshan_debug_var(f_points);
		points_out.load(f_points);
		gproshan_debug_var(points_out.size());

		for(index_t i = 0; i < points_out.size(); ++i)
			mesh.selected.push_back(points_out(i));

		mesh->update_normals();
	}

	return true;
}

bool app_viewer::process_pc_reconstruction(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static msparse_coding::params params;
	static size_t n = 12;
	static real_t percentage_size = 100;
	static real_t radio_factor = 1;

	assert(sizeof(ImGuiDataType_U64) != sizeof(size_t));

	ImGui::InputScalar("basis", ImGuiDataType_U64, &n);
	ImGui::InputScalar("atoms", ImGuiDataType_U64, &params.n_atoms);
	ImGui_InputReal("delta", &params.delta, 0.001, 0.1, "%.3lf");
	ImGui_InputReal("proj_thres", &params.sum_thres, 1.001, 0.1, "%.6lf");
	ImGui_InputReal("area_thres", &params.area_thres, 0.001, 0.1, "%6lf");
	ImGui::Checkbox("learn", &params.learn);

	ImGui_InputReal("percentage_size", &percentage_size, 100, 10, "%.3f");
	ImGui_InputReal("radio_factor", &radio_factor, 1, 0.1, "%.3f");

	if(ImGui::Button("Run"))
	{
		basis_dct phi(n);
		msparse_coding msc(mesh, &phi, params);

		view->add_mesh(msc.point_cloud_reconstruction(percentage_size, radio_factor));
	}

	return true;
}


// Features

bool app_viewer::process_functional_maps(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static int K = 20;

	ImGui::InputInt("eigenvectors", &K);

	if(ImGui::Button("Run"))
	{
		a_sp_mat L, A;
		a_vec eigval;
		a_mat eigvec;

		TIC(view->time) K = eigs_laplacian(mesh, eigval, eigvec, L, A, K); TOC(view->time)
		gproshan_log_var(view->time);

		gproshan_log_var(K);

		K = K < N_MESHES ? K : N_MESHES;
		for(index_t k = 0; k < N_MESHES; ++k)
		{
			if(k) view->add_mesh(new che(*mesh));
			view->idx_active_mesh = k;

			eigvec.col(k) -= eigvec.col(k).min();
			eigvec.col(k) /= eigvec.col(k).max();

			#pragma omp parallel for
			for(index_t v = 0; v < mesh->n_vertices; ++v)
				view->active_mesh()->heatmap(v) = eigvec(v, k);

			view->active_mesh().update_vbo();
		}

		view->idx_active_mesh = 0;
	}

	return true;
}

bool app_viewer::process_descriptor_heatmap(viewer * p_view, const descriptor::signature & sig)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static int K = 50;
	static bool status = true;
	ImGui::InputInt("eigenvectors", &K);
	if(!status) ImGui::TextColored({1, 0, 0, 1}, "Error computing features.");

	if(ImGui::Button("Run"))
	{
		descriptor features(sig, mesh, K);

		if(features)
		{
			status = true;

			real_t max_s = 0;
			#pragma omp parallel for reduction(max: max_s)
			for(index_t v = 0; v < mesh->n_vertices; ++v)
			{
				mesh->heatmap(v) = features(v);
				max_s = max(max_s, mesh->heatmap(v));
			}

			#pragma omp parallel for
			for(index_t v = 0; v < mesh->n_vertices; ++v)
				mesh->heatmap(v) /= max_s;
		}
		else status = false;
	}

	return true;
}

bool app_viewer::process_gps(viewer * p_view)
{
	return process_descriptor_heatmap(p_view, descriptor::GPS);
}

bool app_viewer::process_hks(viewer * p_view)
{
	return process_descriptor_heatmap(p_view, descriptor::HKS);
}

bool app_viewer::process_wks(viewer * p_view)
{
	return process_descriptor_heatmap(p_view, descriptor::WKS);
}

bool app_viewer::process_key_points(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	key_points kps(mesh);

	mesh.selected.clear();
	mesh.selected.reserve(kps.size());

	for(index_t i = 0; i < kps.size(); ++i)
		mesh.selected.push_back(kps[i]);

	return false;
}

bool app_viewer::process_key_components(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	key_points kps(mesh);
	key_components kcs(mesh, kps, .25);

	gproshan_debug_var(kcs);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		mesh->heatmap(v) = (real_t) kcs(v) / kcs;

	return false;
}


// Hole Filling

bool paint_holes_vertices(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	size_t nv = mesh->n_vertices;

	mesh.update();

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		if(v >= nv) mesh->heatmap(v) = .25;

	return false;
}

bool app_viewer::process_poisson(viewer * p_view, const index_t & k)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	size_t old_n_vertices = mesh->n_vertices;
	delete [] fill_all_holes(mesh);

	TIC(view->time) poisson(mesh, old_n_vertices, k); TOC(view->time)
	gproshan_log_var(view->time);

//	paint_holes_vertices();
	mesh.update();

	return false;
}

bool app_viewer::process_poisson_laplacian_1(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_poisson(p_view, 1);
}

bool app_viewer::process_poisson_laplacian_2(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_poisson(p_view, 2);
}

bool app_viewer::process_poisson_laplacian_3(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	return process_poisson(p_view, 3);
}

bool app_viewer::process_fill_holes(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	fill_all_holes(mesh);

	paint_holes_vertices(p_view);

	return false;
}

bool app_viewer::process_fill_holes_biharmonic_splines(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	size_t old_n_vertices, n_vertices = mesh->n_vertices;
	size_t n_holes = 0; // FIX_BOUND mesh->n_borders;

	vector<index_t> * border_vertices;
	che ** holes;
	tie(border_vertices, holes) = fill_all_holes_meshes(mesh);
	if(!holes) return true;

	index_t k = 2;

	for(index_t h = 0; h < n_holes; ++h)
		if(holes[h])
		{
			old_n_vertices = n_vertices;
			biharmonic_interp_2(mesh, old_n_vertices, n_vertices += holes[h]->n_vertices - border_vertices[h].size(), border_vertices[h], k);
			delete holes[h];
		}

	delete [] holes;
	delete [] border_vertices;
	paint_holes_vertices(p_view);

	return false;
}


// Others

bool app_viewer::process_select_multiple(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	static char line[128] = "";

	ImGui::InputText("select", line, sizeof(line));

	if(ImGui::Button("Add"))
	{
		stringstream ss(line);
		index_t v;
		while(ss >> v)
			mesh.selected.push_back(v);
	}

	return true;
}

bool app_viewer::process_threshold(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	for(index_t v = 0; v < mesh->n_vertices; ++v)
		mesh->heatmap(v) = mesh->heatmap(v) > 0.5 ? 1 : 0.5;

	return false;
}

bool app_viewer::process_noise(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		real_t r = real_t(d_mod_1000(generator)) / 200000;
		mesh->get_vertex(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh->update_normals();

	return false;
}

bool app_viewer::process_black_noise(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->active_mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		real_t r = real_t(d_mod_1000(generator)) / 200000;
		mesh->get_vertex(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh->update_normals();

	return false;
}


} // namespace gproshan

