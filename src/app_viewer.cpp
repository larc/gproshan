#include "app_viewer.h"

#include <random>

using namespace std;
using namespace gproshan::mdict;

// geometry processing and shape analysis framework
namespace gproshan {


che * load_mesh(const string & file_path)
{
	size_t pos = file_path.rfind('.');
	assert(pos != string::npos);
	
	string extension = file_path.substr(pos + 1);

	if(extension == "off") return new che_off(file_path);
	if(extension == "obj") return new che_obj(file_path);
	if(extension == "ply") return new che_ply(file_path);

	return new che_img(file_path);
}

app_viewer::app_viewer()
{
	dist = nullptr;
	n_dist = 0;
}

app_viewer::~app_viewer()
{
	if(dist) delete [] dist;

	for(che * mesh: meshes)
		delete mesh;
}

int app_viewer::main(int nargs, const char ** args)
{
	if(nargs < 2)
	{
		printf("./gproshan [mesh_paths.(off,obj,ply)]\n");
		return 0;
	}

	TIC(time)

	vector<che *> meshes;
	for(int i = 1; i < nargs; i++)
		meshes.push_back(load_mesh(args[i]));

	TOC(time)

	gproshan_log_var(sizeof(real_t));
	gproshan_log_var(time);
	
	//init mesher
	add_mesh(meshes);

	sub_menus.push_back("Fairing");
	add_process(GLFW_KEY_T, {"T", "Fairing Taubin", process_fairing_taubin});
	add_process(GLFW_KEY_E, {"E", "Fairing Spectral", process_fairing_spectral});

	sub_menus.push_back("Geodesics");
	add_process(GLFW_KEY_F, {"F", "Fast Marching", process_geodesics_fm});
	add_process(GLFW_KEY_U, {"C", "Parallel Toplesets Propagation CPU", process_geodesics_ptp_cpu});
#ifndef SINGLE_P
	add_process(GLFW_KEY_L, {"L", "Heat Method", process_geodesics_heat_flow});
#endif

#ifdef GPROSHAN_CUDA
	add_process(GLFW_KEY_G, {"G", "Parallel Toplesets Propagation GPU", process_geodesics_ptp_gpu});
//	add_process('L', "Geodesics (HEAT_FLOW_GPU)", process_geodesics_heat_flow_gpu);
#endif // GPROSHAN_CUDA

	add_process(GLFW_KEY_S, {"S", "Geodesic Farthest Point Sampling", process_farthest_point_sampling});
	add_process(GLFW_KEY_R, {"R", "Geodesic Farthest Point Sampling (radio)", process_farthest_point_sampling_radio});
	add_process(GLFW_KEY_V, {"V", "Geodesic Voronoi", process_voronoi});
	add_process(GLFW_KEY_P, {"P", "Toplesets", compute_toplesets});

	sub_menus.push_back("Dictionary Learning");
	add_process(GLFW_KEY_J, {"J", "MDICT Patch", process_mdict_patch});
	add_process(GLFW_KEY_D, {"D", "MDICT Denoising", process_denoising});
	add_process(GLFW_KEY_A, {"A", "MDICT Super Resolution", process_super_resolution});
	add_process(GLFW_KEY_I, {"I", "MDICT Inpaiting", process_inpaiting});
//	add_process('A', "IT Inpainting", process_iterative_inpaiting);

	sub_menus.push_back("Signatures");
	add_process(GLFW_KEY_2, {"2", "GPS", process_gps});
	add_process(GLFW_KEY_3, {"3", "HKS", process_hks});
	add_process(GLFW_KEY_4, {"4", "WKS", process_wks});
	add_process(GLFW_KEY_5, {"5", "Functional Maps", process_functional_maps});
	add_process(GLFW_KEY_6, {"6", "Key Points", process_key_points});
	add_process(GLFW_KEY_7, {"7", "Key Components", process_key_components});

	sub_menus.push_back("Repair Holes");
	add_process(GLFW_KEY_X, {"X", "Poisson Membrane surface", process_poisson_laplacian_1});
	add_process(GLFW_KEY_Y, {"Y", "Poisson Thin-plate surface", process_poisson_laplacian_2});
	add_process(GLFW_KEY_Z, {"Z", "Poisson Minimum variation surface", process_poisson_laplacian_3});
	add_process(GLFW_KEY_H, {"H", "Fill hole - planar mesh", process_fill_holes});
	add_process(GLFW_KEY_B, {"B", "Fill hole - biharmonic splines", process_fill_holes_biharmonic_splines});

	sub_menus.push_back("Others");
	add_process(GLFW_KEY_SLASH, {"SLASH", "Threshold", process_threshold});
	add_process(GLFW_KEY_N, {"N", "Noise", process_noise});
	add_process(GLFW_KEY_COMMA, {"COMMA", "Black noise", process_black_noise});
	add_process(GLFW_KEY_M, {"M", "Multiplicate", process_multiplicate_vertices});
	add_process(GLFW_KEY_PERIOD, {"PERIOD", "Delete vertices", process_delete_vertices});
	add_process(GLFW_KEY_MINUS, {"MINUS", "Delete non-manifold vertices", process_delete_non_manifold_vertices});
	add_process(GLFW_KEY_K, {"K", "Gaussian curvature", process_gaussian_curvature});
	add_process(GLFW_KEY_9, {"9", "Edge Collapse", process_edge_collapse});
	add_process(GLFW_KEY_SEMICOLON, {"SEMICOLON", "Select multiple vertices", select_multiple});


	run();


	return 0;
}

void paint_holes_vertices(viewer * p_view)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t nv = mesh->n_vertices();

	mesh.update();

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		if(v >= nv) mesh.color(v) = .25;
}

void app_viewer::process_delete_non_manifold_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	gproshan_debug(removing vertex);
	mesh->remove_non_manifold_vertices();
	gproshan_debug(removing vertex);
}

void app_viewer::process_delete_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size()) return;
	gproshan_debug(removing vertex);
	mesh->remove_vertices(view->select_vertices);
	view->select_vertices.clear();
	gproshan_debug(removing vertex);
}

void app_viewer::process_poisson(viewer * p_view, const index_t & k)
{
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t old_n_vertices = mesh->n_vertices();
	delete [] fill_all_holes(mesh);

	TIC(view->time) poisson(mesh, old_n_vertices, k); TOC(view->time)
	gproshan_log_var(view->time);

//	paint_holes_vertices();
}

void app_viewer::process_poisson_laplacian_1(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	process_poisson(p_view, 1);
}

void app_viewer::process_poisson_laplacian_2(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	process_poisson(p_view, 2);
}

void app_viewer::process_poisson_laplacian_3(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	process_poisson(p_view, 3);
}

void app_viewer::process_fill_holes(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	fill_all_holes(mesh);

	paint_holes_vertices(p_view);
}

void app_viewer::process_noise(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		distance_t r = distance_t(d_mod_1000(generator)) / 200000;
		mesh->get_vertex(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh.update_normals();
}

void app_viewer::process_black_noise(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	std::default_random_engine generator;
	std::uniform_int_distribution<int> d_mod_5(0, 4);
	std::uniform_int_distribution<int> d_mod_1000(0, 999);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		distance_t r = distance_t(d_mod_1000(generator)) / 200000;
		mesh->get_vertex(v) += (!d_mod_5(generator)) * r * mesh->normal(v);
	}

	mesh.update_normals();
}

void app_viewer::process_threshold(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh.color(v) = mesh.color(v) > 0.5 ? 1 : 0.5;
}

void app_viewer::process_functional_maps(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t K = 100;

	a_sp_mat L, A;

	TIC(view->time) laplacian(mesh, L, A); TOC(view->time)
	gproshan_log_var(view->time);

	a_vec eigval;
	a_mat eigvec;

	TIC(view->time) K = eigs_laplacian(eigval, eigvec, mesh, L, A, K); TOC(view->time)
	gproshan_log_var(view->time);

	gproshan_log_var(K);

	K = K < N_MESHES ? K : N_MESHES;
	for(index_t k = 0; k < N_MESHES; k++)
	{
		if(k) view->add_mesh({new che(*mesh)});
		view->current = k;

		eigvec.col(k) -= eigvec.col(k).min();
		eigvec.col(k) /= eigvec.col(k).max();
	
		#pragma omp parallel for
		for(index_t v = 0; v < mesh->n_vertices(); v++)
			view->mesh().color(v) = eigvec(v, k);
	}
	
	view->current = 0;
}

void app_viewer::process_wks(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t K = 50, T = 100;

	a_sp_mat L, A;

	TIC(view->time) laplacian(mesh, L, A); TOC(view->time)
	gproshan_log_var(view->time);

	a_vec eigval;
	a_mat eigvec;

	TIC(view->time) K = eigs_laplacian(eigval, eigvec, mesh, L, A, K); TOC(view->time)
	gproshan_log_var(view->time);

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		a_vec s(T, arma::fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		mesh.color(v) = norm(s);
		max_s = max(max_s, mesh.color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh.color(v) /= max_s;
}

void app_viewer::process_hks(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t K = 100;
	size_t T = 100;

	a_sp_mat L, A;

	TIC(view->time) laplacian(mesh, L, A); TOC(view->time)
	gproshan_log_var(view->time);

	a_vec eigval;
	a_mat eigvec;

	TIC(view->time) K = eigs_laplacian(eigval, eigvec, mesh, L, A, K); TOC(view->time)
	gproshan_log_var(view->time);

	if(!K) return;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		a_vec s(T, arma::fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-abs(eigval(k)) * t) * eigvec(v, k) * eigvec(v, k);

		mesh.color(v) = norm(abs(arma::fft(s, 128)));
		//mesh.color(v) = norm(s);
		max_s = max(max_s, mesh.color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh.color(v) /= max_s;
}

void app_viewer::process_gps(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t K = 50;

	a_sp_mat L, A;

	TIC(view->time) laplacian(mesh, L, A); TOC(view->time)
	gproshan_log_var(view->time);

	a_vec eigval;
	a_mat eigvec;

	TIC(view->time) K = eigs_laplacian(eigval, eigvec, mesh, L, A, K); TOC(view->time)
	gproshan_log_var(view->time);

	eigvec = abs(eigvec);
	eigvec.col(0).zeros();
	for(index_t i = 1; i < K; i++)
		eigvec.col(i) /= sqrt(abs(eigval(i)));

	a_mat data = eigvec.t();
	a_mat means;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		mesh.color(v) = norm(eigvec.row(v));
			max_s = max(max_s, mesh.color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh.color(v) /= max_s;
}

void app_viewer::process_key_points(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();
	
	key_points kps(mesh);

	view->select_vertices.clear();
	view->select_vertices.reserve(kps.size());

	for(index_t i = 0; i < kps.size(); i++)
		view->select_vertices.push_back(kps[i]);
}

void app_viewer::process_key_components(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();
	
	key_points kps(mesh);
	key_components kcs(mesh, kps, .25);
	
	gproshan_debug_var(kcs);
	
	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
		mesh.color(v) = (real_t) kcs(v) / kcs;
}

void app_viewer::process_mdict_patch(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();
	
	TIC(view->time)
	index_t * toplevel = new index_t[mesh->n_vertices()];
	size_t avg_nvp = 0;

	vertex vdir;
	patch p;
	distance_t mean_edge = mesh->mean_edge();
	for(auto & v: view->select_vertices)
	{
		p.init(mesh, v, dictionary::T, dictionary::T * mean_edge, toplevel);
		for(auto & u: p.vertices)
			mesh.color(u) = 1;

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

	avg_nvp /= view->select_vertices.size();
	gproshan_debug_var(avg_nvp);
	
	delete [] toplevel;
	TOC(view->time)
	gproshan_debug_var(view->time);
}

void app_viewer::process_denoising(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_input(n m M f learn);
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	denoising dict(mesh, phi, m, M, f);
	dict.execute();

	delete phi;
	mesh.update_normals();
}

void app_viewer::process_super_resolution(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_log(parameters: (n, m, M, f, learn));
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	super_resolution dict(mesh, phi, m, M, f);
	dict.execute();

	delete phi;
	mesh.update_normals();
}

void app_viewer::process_inpaiting(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_log(parameters: (n, m, M, f, learn));
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	inpainting dict(mesh, phi, m, M, f);
	dict.execute();

	delete phi;
	mesh.update_normals();
}


void app_viewer::process_iterative_inpaiting(viewer * p_view)
{
	gproshan_log(APP_VIEWER);

//	mesh_iterative_inpaiting(mesh, view->select_vertices, freq, rt, m, M, f, learn);
}

void app_viewer::process_multiplicate_vertices(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	mesh->multiplicate_vertices();
	mesh.log_info();
}

void app_viewer::compute_toplesets(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();
	
	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);

	index_t * toplesets = new index_t[mesh->n_vertices()];
	index_t * sorted = new index_t[mesh->n_vertices()];
	vector<index_t> limites;
	mesh->compute_toplesets(toplesets, sorted, limites, view->select_vertices);

	size_t n_toplesets = limites.size() - 1;

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		if(toplesets[v] < n_toplesets) 
			mesh.color(v) = distance_t(toplesets[v]) / (n_toplesets);
	}

	gproshan_debug_var(n_toplesets);

	delete [] toplesets;
	delete [] sorted;
}

void app_viewer::process_voronoi(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	TIC(view->time)
#ifdef GPROSHAN_CUDA
	geodesics ptp(mesh, view->select_vertices, geodesics::PTP_GPU, nullptr, 1);
#else
	geodesics ptp(mesh, view->select_vertices, geodesics::FM, nullptr, 1);
#endif
	TOC(view->time)
	gproshan_log_var(view->time);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		mesh.color(i) = ptp.clusters[i];
		mesh.color(i) /= view->select_vertices.size() + 1;
	}
}

void app_viewer::process_farthest_point_sampling_radio(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;

	gproshan_input(radio);
	distance_t radio; cin >> radio;

#ifdef GPROSHAN_CUDA	// IMPLEMENT/REVIEW
	double time_fps;

	TIC(view->time)
	radio = farthest_point_sampling_ptp_gpu(mesh, view->select_vertices, time_fps, NIL, radio);
	TOC(view->time)
	gproshan_log_var(time_fps);
#endif // GPROSHAN_CUDA

	gproshan_log_var(radio);
	gproshan_log_var(view->select_vertices.size());
	gproshan_log_var(view->time);
}

void app_viewer::process_farthest_point_sampling(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	gproshan_input(samples_number);
	index_t n; cin >> n;

	distance_t radio;
	TIC(view->time)
	load_sampling(view->select_vertices, radio, mesh, n);
	TOC(view->time)
	gproshan_log_var(view->time);
}

void app_viewer::process_fairing_spectral(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();
	
	gproshan_input(k (eigenvectors number));
	size_t k; cin >> k;
	
	fairing * fair = new fairing_spectral(k);
	fair->run(mesh);

	mesh->set_vertices(fair->get_postions());
	delete fair;

	mesh.update_normals();
}

void app_viewer::process_fairing_taubin(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	gproshan_input(step);
	real_t step; cin >> step;
	
	fairing * fair = new fairing_taubin(step);
	fair->run(mesh);

	mesh->set_vertices(fair->get_postions());
	delete fair;

	mesh.update_normals();
}

void app_viewer::process_geodesics_fm(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);

	TIC(view->time)
	geodesics fm(mesh, view->select_vertices);
	TOC(view->time)
	gproshan_log_var(view->time);

	mesh.update_colors(&fm[0]);
}

void app_viewer::process_geodesics_ptp_cpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);
	
	TIC(view->time)
	geodesics ptp(mesh, view->select_vertices, geodesics::PTP_CPU);
	TOC(view->time)
	gproshan_log_var(view->time);

	mesh.update_colors(&ptp[0]);
}

void app_viewer::process_geodesics_heat_flow(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);
	
	TIC(view->time)
	geodesics heat_flow(mesh, view->select_vertices, geodesics::HEAT_FLOW);
	TOC(view->time)
	gproshan_log_var(view->time);

	mesh.update_colors(&heat_flow[0]);
}


#ifdef GPROSHAN_CUDA

void app_viewer::process_geodesics_ptp_gpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);
	
	if(view->dist && view->n_dist != mesh.n_vertices())
	{
		delete [] view->dist;
		view->n_dist = 0;
		view->dist = nullptr;
	}

	if(!view->dist)
	{
		view->n_dist = mesh.n_vertices();
		view->dist = new distance_t[view->n_dist];
	}

	TIC(view->time)
	geodesics ptp(mesh, view->select_vertices, geodesics::PTP_GPU, view->dist);
	TOC(view->time)
	gproshan_log_var(view->time);
	
	mesh.update_colors(&ptp[0]);
}

void app_viewer::process_geodesics_heat_flow_gpu(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	if(!view->select_vertices.size())
		view->select_vertices.push_back(0);
	
	TIC(view->time)
	geodesics heat_flow(mesh, view->select_vertices, geodesics::HEAT_FLOW_GPU);
	TOC(view->time)
	gproshan_log_var(view->time);

	mesh.update_colors(&heat_flow[0]);
}

#endif // GPROSHAN_CUDA


void app_viewer::process_fill_holes_biharmonic_splines(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	size_t old_n_vertices, n_vertices = mesh.n_vertices();
	size_t n_holes = mesh->n_borders();

	vector<index_t> * border_vertices;
	che ** holes;
	tie(border_vertices, holes) = fill_all_holes_meshes(mesh);
	if(!holes) return;

	index_t k = 2;

	for(index_t h = 0; h < n_holes; h++)
		if(holes[h])
		{
			old_n_vertices = n_vertices;
			biharmonic_interp_2(mesh, old_n_vertices, n_vertices += holes[h]->n_vertices() - border_vertices[h].size(), border_vertices[h], k);
			delete holes[h];
		}

	delete [] holes;
	delete [] border_vertices;
	paint_holes_vertices(p_view);
}

void app_viewer::process_gaussian_curvature(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	real_t g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

	a_vec gv(mesh.n_vertices());

	#pragma omp parallel for private(g, a, b) reduction(max: g_max) reduction(min: g_min)
	for(index_t v = 0; v < mesh.n_vertices(); v++)
	{
		g = 0;
		for_star(he, mesh, v)
		{
			a = mesh->gt_vt(next(he)) - mesh->gt(v);
			b = mesh->gt_vt(prev(he)) - mesh->gt(v);
			g += acos((a,b) / (*a * *b));
		}
		gv(v) = (2 * M_PI - g) / mesh->area_vertex(v);
		g_max = max(g_max, gv(v));
		g_min = min(g_min, gv(v));
	}

	g = g_max - g_min;

	#pragma omp parallel for
	for(index_t v = 0; v < mesh.n_vertices(); v++)
		gv(v) = (gv(v) - g_min) / g;

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
	for(index_t v = 0; v < mesh.n_vertices(); v++)
		mesh.color(v) = f(gv(v));
}

void app_viewer::process_edge_collapse(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;
	che_viewer & mesh = view->mesh();

	index_t levels;
	cin >> levels;

	TIC(view->time) decimation sampling(mesh, mesh.normals_ptr(), levels); TOC(view->time)
	gproshan_debug_var(view->time);

	if(view->n_meshes < 2)
		view->add_mesh({new che(*mesh)});

	view->corr_mesh[1].init(view->meshes[1]->n_vertices(), view->current, sampling);
	view->current = 1;
}

void app_viewer::select_multiple(viewer * p_view)
{
	gproshan_log(APP_VIEWER);
	app_viewer * view = (app_viewer *) p_view;

	char line[128];
	if(fgets(line, 128, stdin))
	{
		stringstream ss(line);
		index_t v;
		while(ss >> v)
			view->select_vertices.push_back(v);
	}
}


} // namespace gproshan

