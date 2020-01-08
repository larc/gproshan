#include "app_viewer.h"

using namespace std;
using namespace gproshan::mdict;


// elapsed time in seconds
double load_time;
distance_t * dist;
size_t n_dist;

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

int viewer_main(int nargs, const char ** args)
{
	if(nargs < 2)
	{
		printf("./gproshan [mesh_paths.(off,obj,ply)]\n");
		return 0;
	}

	TIC(load_time)

	vector<che *> meshes;
	for(int i = 1; i < nargs; i++)
		meshes.push_back(load_mesh(args[i]));

	TOC(load_time)

	gproshan_log_var(sizeof(real_t));
	gproshan_log_var(load_time);

	viewer::sub_menus.push_back("Fairing");
	viewer::add_process('T', "Fairing Taubin", viewer_process_fairing_taubin);
	viewer::add_process('E', "Fairing Spectral", viewer_process_fairing_spectral);

	viewer::sub_menus.push_back("Geodesics");
	viewer::add_process('F', "Geodesics (FM)", viewer_process_geodesics_fm);
	viewer::add_process('U', "Geodesics (PTP_CPU)", viewer_process_geodesics_ptp_cpu);
	#ifndef SINGLE_P
		viewer::add_process('l', "Geodesics (HEAT_FLOW)", viewer_process_geodesics_heat_flow);
	#endif

#ifdef GPROSHAN_CUDA
	viewer::add_process('G', "Geodesics (PTP_GPU)", viewer_process_geodesics_ptp_gpu);
	viewer::add_process('L', "Geodesics (HEAT_FLOW_GPU)", viewer_process_geodesics_heat_flow_gpu);
#endif // GPROSHAN_CUDA

	viewer::add_process('S', "Farthest Point Sampling", viewer_process_farthest_point_sampling);
	viewer::add_process('Q', "Farthest Point Sampling radio", viewer_process_farthest_point_sampling_radio);
	viewer::add_process('V', "Voronoi Regions", viewer_process_voronoi);
	viewer::add_process('P', "Toplesets", viewer_compute_toplesets);

	viewer::sub_menus.push_back("Dictionary Learning");
	viewer::add_process('.', "Mark patch", viewer_process_mdict_patch);
	viewer::add_process('D', "Denoising", viewer_process_denoising);
	viewer::add_process('R', "Super Resolution", viewer_process_super_resolution);
	viewer::add_process('I', "Inpainting", viewer_process_inpaiting);
	viewer::add_process('z', "Load mask", viewer_process_mask);

	viewer::add_process('s', "Synthesis", viewer_process_synthesis);
	viewer::add_process('A', "IT Inpainting", viewer_process_iterative_inpaiting);

	viewer::sub_menus.push_back("Signatures");
	viewer::add_process('s', "GPS (norm)", viewer_process_gps);
	viewer::add_process('H', "HKS (norm)", viewer_process_hks);
	viewer::add_process('W', "WKS (norm)", viewer_process_wks);
	viewer::add_process('X', "Functional maps", viewer_process_functional_maps);
	viewer::add_process('*', "Key Points (adaptive mesh)", viewer_process_key_points);
	viewer::add_process('C', "Key Components", viewer_process_key_components);

	viewer::sub_menus.push_back("Repair Holes");
	viewer::add_process('o', "Membrane surface", viewer_process_poisson_laplacian_1);
	viewer::add_process('p', "Thin-plate surface", viewer_process_poisson_laplacian_2);
	viewer::add_process('q', "Minimum variation surface", viewer_process_poisson_laplacian_3);
	viewer::add_process('h', "Fill Holes (mesh only)", viewer_process_fill_holes);
	viewer::add_process('B', "Fill holes (biharmonic splines)", viewer_process_fill_holes_biharmonic_splines);

	viewer::sub_menus.push_back("Others");
	viewer::add_process('t', "Threshold", viewer_process_thresold);
	viewer::add_process('N', "Noise", viewer_process_noise);
	viewer::add_process('M', "Black Noise", viewer_process_black_noise);
	viewer::add_process('m', "Multiplicate Vertices", viewer_process_multiplicate_vertices);
	viewer::add_process('-', "Make holes", viewer_process_delete_vertices);
	viewer::add_process('d', "Delete non manifolds vertices", viewer_process_delete_non_manifold_vertices);
	viewer::add_process('K', "Gaussian curvature", viewer_process_gaussian_curvature);
	viewer::add_process('/', "Decimation", viewer_process_edge_collapse);
	viewer::add_process(':', "Select multiple vertices", viewer_select_multiple);

	dist = nullptr;
	n_dist = 0;

	//init viewer
	viewer::init(meshes);
	
	if(dist) delete [] dist;
	for(che * mesh: meshes)
		delete mesh;

	return 0;
}

void paint_holes_vertices()
{
	size_t nv = viewer::mesh().n_vertices();

	viewer::mesh().update();

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		if(v >= nv) viewer::vcolor(v) = .25;
}

void viewer_process_delete_non_manifold_vertices()
{
	gproshan_log(APP_VIEWER);

	gproshan_debug(removing vertex);
	viewer::mesh()->remove_non_manifold_vertices();
	gproshan_debug(removing vertex);
}

void viewer_process_delete_vertices()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size()) return;
	gproshan_debug(removing vertex);
	viewer::mesh()->remove_vertices(viewer::select_vertices);
	viewer::select_vertices.clear();
	gproshan_debug(removing vertex);
}

void viewer_process_poisson(const index_t & k)
{
	size_t old_n_vertices = viewer::mesh()->n_vertices();
	delete [] fill_all_holes(viewer::mesh());

	TIC(load_time) poisson(viewer::mesh(), old_n_vertices, k); TOC(load_time)
	gproshan_log_var(load_time);

//	paint_holes_vertices();
}

void viewer_process_poisson_laplacian_1()
{
	gproshan_log(APP_VIEWER);
	viewer_process_poisson(1);
}

void viewer_process_poisson_laplacian_2()
{
	gproshan_log(APP_VIEWER);
	viewer_process_poisson(2);
}

void viewer_process_poisson_laplacian_3()
{
	gproshan_log(APP_VIEWER);
	viewer_process_poisson(3);
}

void viewer_process_fill_holes()
{
	gproshan_log(APP_VIEWER);

	fill_all_holes(viewer::mesh());

	paint_holes_vertices();
}

void viewer_process_noise()
{
	viewer::share = (char *) new vertex;
	delete [] viewer::share;
	gproshan_log(APP_VIEWER);

	srand(time(nullptr));

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		distance_t r = distance_t( rand() % 1000 ) / 200000;
		int p = rand() % 5;
		viewer::mesh()->get_vertex(v) += (!p) * r * viewer::mesh()->normal(v);
	}

	viewer::mesh().update_normals();
}

void viewer_process_black_noise()
{
	gproshan_log(APP_VIEWER);

	srand(time(nullptr));

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		distance_t r = distance_t( rand() % 1000 ) / 200000;
		int p = rand() % 5;
		viewer::mesh()->get_vertex(v) += (!p) * r * viewer::mesh()->normal(v);
		if(!p) viewer::vcolor(v) = INFINITY;
	}

	viewer::mesh().update_normals();
}

void viewer_process_thresold()
{
	gproshan_log(APP_VIEWER);

	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) = viewer::vcolor(v) > 0.5 ? 1 : 0.5;
}

void viewer_process_functional_maps()
{
	gproshan_log(APP_VIEWER);

	size_t K = 100;

	a_sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	gproshan_log_var(load_time);

	a_vec eigval;
	a_mat eigvec;

	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, A, K); TOC(load_time)
	gproshan_log_var(load_time);

	gproshan_log_var(K);

	K = K < N_MESHES ? K : N_MESHES;
	for(index_t k = 0; k < N_MESHES; k++)
	{
		if(k) viewer::add_mesh({new che(*viewer::mesh())});
		viewer::current = k;

		eigvec.col(k) -= eigvec.col(k).min();
		eigvec.col(k) /= eigvec.col(k).max();
	
		#pragma omp parallel for
		for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
			viewer::vcolor(v) = eigvec(v, k);
	}
	
	viewer::current = 0;
}

void viewer_process_wks()
{
	gproshan_log(APP_VIEWER);

	size_t K = 50, T = 100;

	a_sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	gproshan_log_var(load_time);

	a_vec eigval;
	a_mat eigvec;

	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, A, K); TOC(load_time)
	gproshan_log_var(load_time);

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		a_vec s(T, arma::fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		viewer::vcolor(v) = norm(s);
		max_s = max(max_s, viewer::vcolor(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) /= max_s;
}

void viewer_process_hks()
{
	gproshan_log(APP_VIEWER);

	size_t K = 100;
	size_t T = 100;

	a_sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	gproshan_log_var(load_time);

	a_vec eigval;
	a_mat eigvec;

	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, A, K); TOC(load_time)
	gproshan_log_var(load_time);

	if(!K) return;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		a_vec s(T, arma::fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-abs(eigval(k)) * t) * eigvec(v, k) * eigvec(v, k);

		viewer::vcolor(v) = norm(abs(arma::fft(s, 128)));
		//viewer::vcolor(v) = norm(s);
		max_s = max(max_s, viewer::vcolor(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) /= max_s;
}

void viewer_process_gps()
{
	gproshan_log(APP_VIEWER);

	size_t K = 50;

	a_sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	gproshan_log_var(load_time);

	a_vec eigval;
	a_mat eigvec;

	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, A, K); TOC(load_time)
	gproshan_log_var(load_time);

	eigvec = abs(eigvec);
	eigvec.col(0).zeros();
	for(index_t i = 1; i < K; i++)
		eigvec.col(i) /= sqrt(abs(eigval(i)));

	a_mat data = eigvec.t();
	a_mat means;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		viewer::vcolor(v) = norm(eigvec.row(v));
			max_s = max(max_s, viewer::vcolor(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) /= max_s;
}

void viewer_process_key_points()
{
	gproshan_log(APP_VIEWER);
	
	key_points kps(viewer::mesh());

	viewer::select_vertices.clear();
	viewer::select_vertices.reserve(kps.size());

	for(index_t i = 0; i < kps.size(); i++)
		viewer::select_vertices.push_back(kps[i]);
}

void viewer_process_key_components()
{
	gproshan_log(APP_VIEWER);
	
	key_points kps(viewer::mesh());
	key_components kcs(viewer::mesh(), kps, .25);
	
	gproshan_debug_var(kcs);
	
	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) = (real_t) kcs(v) / kcs;
}

void viewer_process_mdict_patch()
{
	gproshan_log(APP_VIEWER);
	
	TIC(load_time)
	che * mesh = viewer::mesh();
	index_t * toplevel = new index_t[mesh->n_vertices()];
	size_t avg_nvp = 0;

	vertex vdir;
	patch p;
	distance_t mean_edge = mesh->mean_edge();
	for(auto & v: viewer::select_vertices)
	{
		p.init(mesh, v, dictionary::T, dictionary::T * mean_edge, toplevel);
		for(auto & u: p.vertices)
			viewer::vcolor(u) = 1;

		vdir.x = p.T(0, 0);
		vdir.y = p.T(0, 1);
		vdir.z = p.T(0, 2);
		viewer::vectors.push_back(mesh->gt(v));
		viewer::vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);
		
		vdir.x = p.T(1, 0);
		vdir.y = p.T(1, 1);
		vdir.z = p.T(1, 2);
		viewer::vectors.push_back(mesh->gt(v));
		viewer::vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);
		
		vdir.x = p.T(2, 0);
		vdir.y = p.T(2, 1);
		vdir.z = p.T(2, 2);
		viewer::vectors.push_back(mesh->gt(v));
		viewer::vectors.push_back(mesh->gt(v) + 3 * mean_edge * vdir);
		
		viewer::vectors.push_back(mesh->gt(v));
		viewer::vectors.push_back(mesh->gt(v) + 3 * mean_edge * mesh->normal(v));

		avg_nvp += p.vertices.size();
	}

	avg_nvp /= viewer::select_vertices.size();
	gproshan_debug_var(avg_nvp);
	
	delete [] toplevel;
	TOC(load_time)
	gproshan_debug_var(load_time);
}

void viewer_process_denoising()
{
	gproshan_log(APP_VIEWER);

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_input(n m M f learn);
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	denoising dict(viewer::mesh(), phi, m, M, f, learn);
	dict.execute();
	
	delete phi;
	viewer::mesh().update_colors(&dict[0]);
	
	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::vcolor(v) = 2 * atan(viewer::vcolor(v) * 10) / M_PI;
	viewer::mesh().update_normals();
}

void viewer_process_super_resolution()
{
	gproshan_log(APP_VIEWER);

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_log(parameters: (n, m, M, f, learn));
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	super_resolution dict(viewer::mesh(), phi, m, M, f, learn);
	dict.execute();

	delete phi;
	viewer::mesh().update_normals();
}

void viewer_process_inpaiting()
{
	gproshan_log(APP_VIEWER);

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;
	size_t avg_p; 
	size_t percentage;

	gproshan_input(n m M f learn avg_p percentage);
	cin >> n >> m >> M >> f >> learn >> avg_p >> percentage;

	basis * phi = new basis_dct(n);
	inpainting dict(viewer::mesh(), phi, m, M, f, learn, avg_p, percentage);
	dict.execute();
	delete phi;
	viewer::mesh().update_colors(&dict[0]);
	
	viewer::mesh().update_normals();
}

void viewer_process_mask()
{
	gproshan_log(APP_VIEWER);

	size_t avg_p; 
	size_t percentage;

	size_t n=12; // dct
	size_t m = 144, M = 0;
	distance_t f = 1;
	bool learn = 0;


	gproshan_input(avg_p percentage );
	cin >> avg_p >> percentage;

	basis * phi = new basis_dct(n);
	inpainting dict(viewer::mesh(),  phi, m, M, f, learn, avg_p, percentage);

	dict.init_radial_patches(1.663e-02);
	//dict.init_voronoi_patches();
	delete phi;
	viewer::mesh().update_colors(&dict[0]);
	
	viewer::mesh().update_normals();
}



void viewer_process_synthesis()
{
	gproshan_log(APP_VIEWER);

	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	gproshan_log(parameters: (n, m, M, f, learn));
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	synthesis dict(viewer::mesh(), phi, m, M, f, learn);
	dict.execute();

	delete phi;
	viewer::mesh().update_normals();
}



void viewer_process_iterative_inpaiting()
{
	gproshan_log(APP_VIEWER);

//	mesh_iterative_inpaiting(viewer::mesh(), viewer::select_vertices, freq, rt, m, M, f, learn);
}

void viewer_process_multiplicate_vertices()
{
	gproshan_log(APP_VIEWER);

	viewer::mesh()->multiplicate_vertices();
	viewer::mesh().log_info();
}

void viewer_compute_toplesets()
{
	gproshan_log(APP_VIEWER);
	
	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);

	index_t * toplesets = new index_t[viewer::mesh()->n_vertices()];
	index_t * sorted = new index_t[viewer::mesh()->n_vertices()];
	vector<index_t> limites;
	viewer::mesh()->compute_toplesets(toplesets, sorted, limites, viewer::select_vertices);

	size_t n_toplesets = limites.size() - 1;

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		if(toplesets[v] < n_toplesets) 
			viewer::vcolor(v) = distance_t(toplesets[v]) / (n_toplesets);
	}

	gproshan_debug_var(n_toplesets);

	delete [] toplesets;
	delete [] sorted;
}

void viewer_process_voronoi()
{
	gproshan_log(APP_VIEWER);

	TIC(load_time)
#ifdef GPROSHAN_CUDA
	geodesics ptp(viewer::mesh(), viewer::select_vertices, geodesics::PTP_GPU, nullptr, 1);
#else
	geodesics ptp(viewer::mesh(), viewer::select_vertices, geodesics::FM, nullptr, 1);
#endif
	TOC(load_time)
	gproshan_log_var(load_time);

	#pragma omp parallel for
	for(index_t i = 0; i < viewer::mesh()->n_vertices(); i++)
	{
		viewer::vcolor(i) = ptp.clusters[i];
		viewer::vcolor(i) /= viewer::select_vertices.size() + 1;
	}
}

void viewer_process_farthest_point_sampling_radio()
{
	gproshan_log(APP_VIEWER);

	gproshan_input(radio);
	distance_t radio; cin >> radio;

#ifdef GPROSHAN_CUDA	// IMPLEMENT/REVIEW
	double time_fps;

	TIC(load_time)
	radio = farthest_point_sampling_ptp_gpu(viewer::mesh(), viewer::select_vertices, time_fps, NIL, radio);
	TOC(load_time)
	gproshan_log_var(time_fps);
#endif // GPROSHAN_CUDA

	gproshan_log_var(radio);
	gproshan_log_var(viewer::select_vertices.size());
	gproshan_log_var(load_time);
}

void viewer_process_farthest_point_sampling()
{
	gproshan_log(APP_VIEWER);

	gproshan_input(samples_number);
	index_t n; cin >> n;

	distance_t radio;
	TIC(load_time)
	load_sampling(viewer::select_vertices, radio, viewer::mesh(), n);
	TOC(load_time)
	gproshan_log_var(load_time);
}

void viewer_process_fairing_spectral()
{
	gproshan_log(APP_VIEWER);
	
	gproshan_input(k (eigenvectors number));
	size_t k; cin >> k;
	
	fairing * fair = new fairing_spectral(k);
	fair->run(viewer::mesh());

	viewer::mesh()->set_vertices(fair->get_postions());
	delete fair;

	viewer::mesh().update_normals();
}

void viewer_process_fairing_taubin()
{
	gproshan_log(APP_VIEWER);

	gproshan_input(step);
	real_t step; cin >> step;
	
	fairing * fair = new fairing_taubin(step);
	fair->run(viewer::mesh());

	viewer::mesh()->set_vertices(fair->get_postions());
	delete fair;

	viewer::mesh().update_normals();
}

void viewer_process_geodesics_fm()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);

	TIC(load_time)
	geodesics fm(viewer::mesh(), viewer::select_vertices);
	TOC(load_time)
	gproshan_log_var(load_time);

	viewer::mesh().update_colors(&fm[0]);
}

void viewer_process_geodesics_ptp_cpu()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);
	
	TIC(load_time)
	geodesics ptp(viewer::mesh(), viewer::select_vertices, geodesics::PTP_CPU);
	TOC(load_time)
	gproshan_log_var(load_time);

	viewer::mesh().update_colors(&ptp[0]);
}

void viewer_process_geodesics_heat_flow()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);
	
	TIC(load_time)
	geodesics heat_flow(viewer::mesh(), viewer::select_vertices, geodesics::HEAT_FLOW);
	TOC(load_time)
	gproshan_log_var(load_time);

	viewer::mesh().update_colors(&heat_flow[0]);
}


#ifdef GPROSHAN_CUDA

void viewer_process_geodesics_ptp_gpu()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);
	
	if(dist && n_dist != viewer::mesh().n_vertices())
	{
		delete [] dist;
		n_dist = 0;
		dist = nullptr;
	}

	if(!dist)
	{
		n_dist = viewer::mesh().n_vertices();
		dist = new distance_t[n_dist];
	}

	TIC(load_time)
	geodesics ptp(viewer::mesh(), viewer::select_vertices, geodesics::PTP_GPU, dist);
	TOC(load_time)
	gproshan_log_var(load_time);
	
	viewer::mesh().update_colors(&ptp[0]);
}

void viewer_process_geodesics_heat_flow_gpu()
{
	gproshan_log(APP_VIEWER);

	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);
	
	TIC(load_time)
	geodesics heat_flow(viewer::mesh(), viewer::select_vertices, geodesics::HEAT_FLOW_GPU);
	TOC(load_time)
	gproshan_log_var(load_time);

	viewer::mesh().update_colors(&heat_flow[0]);
}

#endif // GPROSHAN_CUDA


void viewer_process_fill_holes_biharmonic_splines()
{
	gproshan_log(APP_VIEWER);

	size_t old_n_vertices, n_vertices = viewer::mesh().n_vertices();
	size_t n_holes = viewer::mesh()->n_borders();

	vector<index_t> * border_vertices;
	che ** holes;
	tie(border_vertices, holes) = fill_all_holes_meshes(viewer::mesh());
	if(!holes) return;

	index_t k = 2;

	for(index_t h = 0; h < n_holes; h++)
		if(holes[h])
		{
			old_n_vertices = n_vertices;
			biharmonic_interp_2(viewer::mesh(), old_n_vertices, n_vertices += holes[h]->n_vertices() - border_vertices[h].size(), border_vertices[h], k);
			delete holes[h];
		}

	delete [] holes;
	delete [] border_vertices;
	paint_holes_vertices();
}

void viewer_process_gaussian_curvature()
{
	gproshan_log(APP_VIEWER);

	real_t g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

	a_vec gv(viewer::mesh().n_vertices());

	#pragma omp parallel for private(g, a, b) reduction(max: g_max) reduction(min: g_min)
	for(index_t v = 0; v < viewer::mesh().n_vertices(); v++)
	{
		g = 0;
		for_star(he, viewer::mesh(), v)
		{
			a = viewer::mesh()->gt_vt(next(he)) - viewer::mesh()->gt(v);
			b = viewer::mesh()->gt_vt(prev(he)) - viewer::mesh()->gt(v);
			g += acos((a,b) / (*a * *b));
		}
		gv(v) = (2 * M_PI - g) / viewer::mesh()->area_vertex(v);
		g_max = max(g_max, gv(v));
		g_min = min(g_min, gv(v));
	}

	g = g_max - g_min;

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh().n_vertices(); v++)
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
	for(index_t v = 0; v < viewer::mesh().n_vertices(); v++)
		viewer::vcolor(v) = f(gv(v));
}

void viewer_process_edge_collapse()
{
	gproshan_log(APP_VIEWER);

	index_t levels;
	cin >> levels;

	TIC(load_time) decimation sampling(viewer::mesh(), viewer::mesh().normals_ptr(), levels); TOC(load_time)
	gproshan_debug_var(load_time);

	if(viewer::n_meshes < 2)
		viewer::add_mesh({new che(*viewer::mesh())});

	viewer::corr_mesh[1].init(viewer::meshes[1]->n_vertices(), viewer::current, sampling);
	viewer::current = 1;
}

void viewer_select_multiple()
{
	gproshan_log(APP_VIEWER);

	char line[128];
	if(fgets(line, 128, stdin))
	{
		stringstream ss(line);
		index_t v;
		while(ss >> v)
			viewer::select_vertices.push_back(v);
	}
}

