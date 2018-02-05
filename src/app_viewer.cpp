#include "app_viewer.h"

using namespace mdict;

// elapsed time in seconds
float load_time;

int viewer_main(int nargs, char ** args)
{
	if(nargs < 2) return 0;

	TIC(load_time)
	vector<che *> meshes;
	for(int i = 1; i < nargs; i++)
		meshes.push_back(new che_off(args[i]));
	TOC(load_time)
	debug(load_time)	

	viewer::sub_menus.push_back("Fairing");
	viewer::add_process('T', "Fairing Taubin", viewer_process_fairing_taubin);
	viewer::add_process('E', "Fairing Spectral", viewer_process_fairing_spectral);

	viewer::sub_menus.push_back("Fast Marching");
	viewer::add_process('F', "Fast Marching", viewer_process_fastmarching);
	viewer::add_process('G', "Geodesics (fm_che)", viewer_process_geodesics);
	viewer::add_process('S', "Farthest Point Sampling", viewer_process_farthest_point_sampling);
	viewer::add_process('Q', "Farthest Point Sampling radio", viewer_process_farthest_point_sampling_radio);
	viewer::add_process('C', "Geodesics (GPU - fastmarching)", viewer_process_fastmarching_gpu);
	viewer::add_process('U', "Geodesics (CPU - fastmarching)", viewer_process_fastmarching_cpu);
	viewer::add_process('V', "Voronoi Regions", viewer_process_voronoi);
	viewer::add_process('P', "Rings propagations", viewer_sort_by_rings);
	
	viewer::sub_menus.push_back("Dictionary Learning");
	viewer::add_process('D', "Denoising", viewer_process_denoising);
	viewer::add_process('R', "Super Resolution", viewer_process_super_resolution);
	viewer::add_process('I', "Inpaiting", viewer_process_inpaiting);
	viewer::add_process('A', "IT Inpaiting", viewer_process_iterative_inpaiting);
	
	viewer::sub_menus.push_back("Signatures");
	viewer::add_process('s', "GPS (norm)", viewer_process_gps);
	viewer::add_process('H', "HKS (norm)", viewer_process_hks);
	viewer::add_process('W', "WKS (norm)", viewer_process_wks);

	viewer::sub_menus.push_back("Poisson");
	viewer::add_process('o', "Membrane surface", viewer_process_poisson_laplacian_1);
	viewer::add_process('p', "Thin-plate surface", viewer_process_poisson_laplacian_2);
	viewer::add_process('q', "Minimum variation surface", viewer_process_poisson_laplacian_3);
	
	viewer::sub_menus.push_back("Others");
	viewer::add_process('t', "Threshold", viewer_process_thresold);
	viewer::add_process('N', "Noise", viewer_process_noise);
	viewer::add_process('m', "Multiplicate Vertices", viewer_process_multiplicate_vertices);
	viewer::add_process('h', "Fill Holes", viewer_process_fill_holes);
	viewer::add_process('-', "Make holes", viewer_process_delete_vertices);
	viewer::add_process('d', "Delete non manifolds vertices", viewer_process_delete_non_manifold_vertices);
	viewer::add_process('B', "Fill holes (biharmonic splines)", viewer_process_fill_holes_biharmonic_splines);
	viewer::add_process('K', "Gaussian curvature", viewer_process_gaussian_curvature);
	viewer::add_process('/', "Decimation", viewer_process_edge_collapse);
	viewer::add_process(':', "Select multiple vertices", viewer_select_multiple);
	
	//init viewer	
	viewer::init(meshes);
		
	return 0;
}

void paint_holes_vertices()
{
	size_t nv = viewer::mesh().n_vertices();
	
	viewer::mesh().update();

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		if(v >= nv) viewer::get_color(v) = .25;
}

void viewer_process_delete_non_manifold_vertices()
{
	debug_me(APP_VIEWER)

	debug_me(removing vertex);
	viewer::mesh()->remove_non_manifold_vertices();
	debug_me(removing vertex);
}

void viewer_process_delete_vertices()
{
	debug_me(APP_VIEWER)

	if(!viewer::select_vertices.size()) return;
	debug_me(removing vertex);
	viewer::mesh()->remove_vertices(viewer::select_vertices);
	viewer::select_vertices.clear();
	debug_me(removing vertex);
}

void viewer_process_poisson(const index_t & k)
{
	size_t old_n_vertices = viewer::mesh()->n_vertices();
	delete [] fill_all_holes(viewer::mesh());
	
	TIC(load_time) poisson(viewer::mesh(), old_n_vertices, k); TOC(load_time)
	debug(load_time)
	
	paint_holes_vertices();
}

void viewer_process_poisson_laplacian_1()
{
	debug_me(APP_VIEWER)
	viewer_process_poisson(1);
}

void viewer_process_poisson_laplacian_2()
{
	debug_me(APP_VIEWER)
	viewer_process_poisson(2);
}

void viewer_process_poisson_laplacian_3()
{
	debug_me(APP_VIEWER)
	viewer_process_poisson(3);
}

void viewer_process_fill_holes()
{
	debug_me(APP_VIEWER)
	
	fill_all_holes(viewer::mesh());
	
	paint_holes_vertices();
}

void viewer_process_noise()
{
	debug_me(APP_VIEWER)

	srand(time(NULL));

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		distance_t r = distance_t( rand() % 1000 ) / 200000;
		int p = rand() % 5;
		viewer::mesh()->get_vertex(v) += (!p) * r * viewer::mesh()->normal(v);
	}

	viewer::mesh().update_normals();
}

void viewer_process_thresold()
{
	debug_me(APP_VIEWER)

	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::get_color(v) = viewer::get_color(v) > 0.5 ? 1 : 0.5;
}

void viewer_process_wks()
{
	debug_me(APP_VIEWER)
	
	size_t K = 20, T = 10;
	
	sp_mat L, A;
	
	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	debug(load_time)
	
	vec eigval;
	mat eigvec;
	
	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, K); TOC(load_time)
	debug(load_time)

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		vec s(T, fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		viewer::get_color(v) = norm(s);
		max_s = max(max_s, viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::get_color(v) /= max_s;
}

void viewer_process_hks()
{
	debug_me(APP_VIEWER)
	
	size_t K = 20, T = 100;
	
	sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	debug(load_time)

	vec eigval;
	mat eigvec;
	
	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, K); TOC(load_time)
	debug(load_time)

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		vec s(T, fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		viewer::get_color(v) = norm(s);
		max_s = max(max_s, viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::get_color(v) /= max_s;
}

void viewer_process_gps()
{
	debug_me(APP_VIEWER)
	
	size_t K = 20;
	
	sp_mat L, A;

	TIC(load_time) laplacian(viewer::mesh(), L, A); TOC(load_time)
	debug(load_time)

	vec eigval;
	mat eigvec;
	
	TIC(load_time) K = eigs_laplacian(eigval, eigvec, viewer::mesh(), L, K); TOC(load_time)
	debug(load_time)
	
	eigvec.col(0).zeros();
	for(index_t i = 1; i < K; i++)
		eigvec.col(i) /= sqrt(eigval(i));

	mat data = eigvec.t();
	mat means;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		viewer::get_color(v) = norm(eigvec.row(v));
			max_s = max(max_s, viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		viewer::get_color(v) /= max_s;
}

void viewer_process_denoising()
{
	debug_me(APP_VIEWER)
	
	size_t freq, rt; // cosine
	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	d_message(parameters: (n, m, M, f, learn))
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	denoising dict(viewer::mesh(), phi, m, M, f);
	dict.execute();

//	mesh_denoising(viewer::mesh(), viewer::select_vertices, freq, rt, m, M, f, learn);

	delete phi;
	viewer::mesh().update_normals();
}

void viewer_process_super_resolution()
{
	debug_me(APP_VIEWER)
	
	size_t freq, rt; // cosine
	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	d_message(parameters: (n, m, M, f, learn))
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	super_resolution dict(viewer::mesh(), phi, m, M, f);
	dict.execute();

//	mesh_super_resolution(viewer::mesh(), viewer::select_vertices, freq, rt, m, M, f, learn);
	
	delete phi;
	viewer::mesh().update_normals();
}	

void viewer_process_inpaiting()
{
	debug_me(APP_VIEWER)
	
	size_t freq, rt; // cosine
	size_t n; // dct
	size_t m, M;
	distance_t f;
	bool learn;

	d_message(parameters: (n, m, M, f, learn))
	cin >> n >> m >> M >> f >> learn;

	basis * phi = new basis_dct(n);
	inpainting dict(viewer::mesh(), phi, m, M, f);
	dict.execute();

//	mesh_inpainting(viewer::mesh(), viewer::select_vertices, freq, rt, m, M, f, learn);
	
	delete phi;
	viewer::mesh().update_normals();
}


void viewer_process_iterative_inpaiting()
{
	debug_me(APP_VIEWER)
	
//	mesh_iterative_inpaiting(viewer::mesh(), viewer::select_vertices, freq, rt, m, M, f, learn);
}

void viewer_process_multiplicate_vertices()
{
	debug_me(APP_VIEWER)

	viewer::mesh()->multiplicate_vertices();
	
	debug(viewer::mesh()->n_vertices())
}

void viewer_sort_by_rings()
{
	debug_me(APP_VIEWER)
	
	index_t * rings = new index_t[viewer::mesh()->n_vertices()];
	index_t * sorted = new index_t[viewer::mesh()->n_vertices()];
	vector<index_t> limites;
	viewer::mesh()->sort_by_rings(rings, sorted, limites, viewer::select_vertices);
	
	size_t k = limites.size() - 1;

	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
	{
		if(rings[v] < k) //viewer::select_vertices.push_back(v);
		viewer::get_color(v) = distance_t(rings[v]) / (limites.size() - 1);
	}
	debug(k)

	delete [] rings;
	delete [] sorted;
}

void viewer_process_fastmarching_cpu()
{
	if(!viewer::select_vertices.size()) return;
	
	debug_me(APP_VIEWER)
	
	index_t * rings = new index_t[viewer::mesh()->n_vertices()];
	index_t * sorted = new index_t[viewer::mesh()->n_vertices()];
	vector<index_t> limites;
	viewer::mesh()->sort_by_rings(rings, sorted, limites, viewer::select_vertices);

	distance_t * distances = parallel_fastmarching(viewer::mesh(), viewer::select_vertices.data(), viewer::select_vertices.size(), load_time, limites, sorted, true, NULL, false);	
	//distance_t * distances = parallel_fastmarching(viewer::mesh()->filename().c_str(), viewer::select_vertices.data(), viewer::select_vertices.size(), load_time, 9, true, true);

	debug(load_time)

	viewer::mesh().update_colors(distances);

	delete [] distances;
	delete [] rings;
	delete [] sorted;
}
void viewer_process_fastmarching_gpu()
{
	if(!viewer::select_vertices.size()) return;
	
	debug_me(APP_VIEWER)
	
	index_t * rings = new index_t[viewer::mesh()->n_vertices()];
	index_t * sorted = new index_t[viewer::mesh()->n_vertices()];
	vector<index_t> limites;
	viewer::mesh()->sort_by_rings(rings, sorted, limites, viewer::select_vertices);
	
	distance_t * distances = parallel_fastmarching(viewer::mesh(), viewer::select_vertices.data(), viewer::select_vertices.size(), load_time, limites, sorted, true);	
	//distance_t * distances = parallel_fastmarching(viewer::mesh()->filename().c_str(), viewer::select_vertices.data(), viewer::select_vertices.size(), load_time, 9, true, true);
	//distance_t * distances = fast_geodesics(viewer::mesh(), viewer::select_vertices.data(), viewer::select_vertices.size(), limites, sorted);
	debug(load_time)
	
	distance_t max_d = 0;

	#pragma omp parallel for reduction(max: max_d)
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		max_d = max(max_d, distances[v]);
	
	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh()->n_vertices(); v++)
		distances[v] /= max_d;

	viewer::mesh().update_colors(distances);
	
	delete [] distances;
	delete [] rings;
	delete [] sorted;
}

void viewer_process_voronoi()
{
	if(!viewer::select_vertices.size()) return;
	
	debug_me(APP_VIEWER)
	
	index_t * rings = new index_t[viewer::mesh()->n_vertices()];
	index_t * sorted = new index_t[viewer::mesh()->n_vertices()];
	vector<index_t> limites;
	viewer::mesh()->sort_by_rings(rings, sorted, limites, viewer::select_vertices);
	
	index_t * clusters = new index_t[viewer::mesh()->n_vertices()];
	memset(clusters, 255, sizeof(index_t) * viewer::mesh()->n_vertices());
	
	distance_t * distances = parallel_fastmarching(viewer::mesh(), viewer::select_vertices.data(), viewer::select_vertices.size(), load_time, limites, sorted, true, clusters);	
	
	debug(load_time)

	#pragma omp parallel for
	for(index_t i = 0; i < viewer::mesh()->n_vertices(); i++)
	{
		viewer::get_color(i) = clusters[i];
		viewer::get_color(i) /= viewer::select_vertices.size() + 1;
	}

	delete [] distances;
	delete [] rings;
	delete [] sorted;
	delete [] clusters;
}

void viewer_process_farthest_point_sampling_radio()
{
	debug_me(APP_VIEWER)
	
	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);

	distance_t radio;
	cin>>radio;

	float load_time_g;
		
	TIC(load_time)
	radio = farthest_point_sampling_gpu(viewer::select_vertices, load_time_g, viewer::mesh(), viewer::mesh()->n_vertices(), radio);
	TOC(load_time)
	
	debug(radio)
	debug(viewer::select_vertices.size())
	debug(load_time)
	debug(load_time_g)
}

void viewer_process_farthest_point_sampling()
{
	debug_me(APP_VIEWER)
	
	if(!viewer::select_vertices.size())
		viewer::select_vertices.push_back(0);

	index_t n;
	cin>>n;

	distance_t radio;
	TIC(load_time)
	load_sampling(viewer::select_vertices, radio, viewer::mesh(), n);
//	distance_t radio = farthest_point_sampling_gpu(viewer::select_vertices, load_time_g, viewer::mesh(), n);
	TOC(load_time)
}

void viewer_process_fairing_spectral()
{
	debug_me(APP_VIEWER)

	fairing * fair = new fairing_spectral(50);
	fair->run(viewer::mesh());
	viewer::mesh()->set_vertices(fair->get_postions());
	delete fair;
	
	viewer::mesh().update_normals();
}

void viewer_process_fairing_taubin()
{
	debug_me(APP_VIEWER)

	fairing * fair = new fairing_taubin;
	fair->run(viewer::mesh());
	viewer::mesh()->set_vertices(fair->get_postions());
	delete fair;
	
	viewer::mesh().update_normals();
}

void viewer_process_geodesics()
{
	debug_me(APP_VIEWER)
	
	TIC(load_time)
	geodesics geodesic(viewer::mesh(), viewer::select_vertices);
	TOC(load_time)
	debug(load_time)

//	geodesic.path_to(viewer::other_vertices, viewer::mesh(), 0);
	
	geodesic.normalize();
	viewer::mesh().update_colors(geodesic.distances);
/* heat map one color
	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh().n_vertices(); v++)
		viewer::get_color(v) = viewer::get_color(v) * 0.5;	
*/		
		
}

void viewer_process_fastmarching()
{
	debug_me(APP_VIEWER)
	off shape(viewer::mesh()->filename());

	#pragma omp parallel for
	for(index_t v = 0; v < shape.get_nvertices(); v++)
		shape(v) = viewer::mesh()->get_vertex(v);

	TIC(load_time)	
	fastmarching fm(shape, viewer::select_vertices, INFINITY, true, false);
	TOC(load_time)
	debug(load_time)

	viewer::mesh().update_colors(fm.distances);
}

void viewer_process_fill_holes_biharmonic_splines()
{
	debug_me(APP_VIEWER)
	
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
	debug_me(APP_VIEWER)
	
	vertex_t g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

#ifdef SINGLE_P
	fvec gv(viewer::mesh().n_vertices());
#else
	vec gv(viewer::mesh().n_vertices());
#endif

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

	vertex_t gm = mean(gv);
	vertex_t gs = var(gv);
	
	debug(gm)
	debug(gs)

	auto f = [&](vertex_t x, vertex_t a = 4) -> vertex_t
	{
		if(x < gm - a * gs) return 0;
		if(x > gm + a * gs) return 1;
		return (x - gm) / (2 * a * gs) + 0.5;
	};

	#pragma omp parallel for
	for(index_t v = 0; v < viewer::mesh().n_vertices(); v++)
		viewer::get_color(v) = f(gv(v));
}

void viewer_process_edge_collapse()
{
	debug_me(APP_VIEWER)
	
	index_t levels;
	cin >> levels;

	TIC(load_time) decimation sampling(viewer::mesh(), viewer::mesh().normals_ptr(), levels); TOC(load_time)
	debug(load_time)

	if(viewer::n_meshes < 2)
		viewer::add_mesh({new che_off(viewer::mesh()->filename())});
	
	viewer::corr_mesh[1].init(viewer::meshes[1]->n_vertices(), viewer::current, sampling);
	viewer::current = 1;
}

void viewer_select_multiple()
{
	debug_me(APP_VIEWER)
	
	char line[128];
	if(fgets(line, 128, stdin))
	{
		stringstream ss(line);
		index_t v;
		while(ss >> v)
			viewer::select_vertices.push_back(v);
	}
}

