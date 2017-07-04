#include "app_viewer.h"

using namespace DDG;

int viewer_main(int nargs, char ** args)
{
	if(nargs < 2) return 0;

	double time;

	string file = args[1];
	
	TIC(time)
	che * shape_che = new che_off(file);
	TOC(time)

	debug(time)	

	Viewer::sub_menus.push_back("Fairing");
	Viewer::add_process('T', "Fairing Taubin", viewer_process_fairing_taubin);
	Viewer::add_process('E', "Fairing Spectral", viewer_process_fairing_spectral);

	Viewer::sub_menus.push_back("Fast Marching");
	Viewer::add_process('F', "Fast Marching", viewer_process_fastmarching);
	Viewer::add_process('G', "Geodesics (fm_che)", viewer_process_geodesics);
	Viewer::add_process('S', "Farthest Point Sampling", viewer_process_farthest_point_sampling);
	Viewer::add_process('Q', "Farthest Point Sampling radio", viewer_process_farthest_point_sampling_radio);
	Viewer::add_process('C', "Geodesics (GPU - fastmarching)", viewer_process_fastmarching_gpu);
	Viewer::add_process('U', "Geodesics (CPU - fastmarching)", viewer_process_fastmarching_cpu);
	Viewer::add_process('V', "Voronoi Regions", viewer_process_voronoi);
	Viewer::add_process('P', "Rings propagations", viewer_sort_by_rings);
	
	Viewer::sub_menus.push_back("Dictionary Learning");
	Viewer::add_process('D', "Denoising", viewer_process_denoising);
	Viewer::add_process('R', "Super Resolution", viewer_process_super_resolution);
	Viewer::add_process('I', "Inpaiting", viewer_process_inpaiting);
	Viewer::add_process('A', "IT Inpaiting", viewer_process_iterative_inpaiting);
	
	Viewer::sub_menus.push_back("Signatures");
	Viewer::add_process('s', "GPS (norm)", viewer_process_gps);
	Viewer::add_process('H', "HKS (norm)", viewer_process_hks);
	Viewer::add_process('W', "WKS (norm)", viewer_process_wks);

	Viewer::sub_menus.push_back("Poisson");
	Viewer::add_process('o', "Membrane surface", viewer_process_poisson_laplacian_1);
	Viewer::add_process('p', "Thin-plate surface", viewer_process_poisson_laplacian_2);
	Viewer::add_process('q', "Minimum variation surface", viewer_process_poisson_laplacian_3);
	
	Viewer::sub_menus.push_back("Others");
	Viewer::add_process('t', "Threshold", viewer_process_thresold);
	Viewer::add_process('N', "Noise", viewer_process_noise);
	Viewer::add_process('m', "Multiplicate Vertices", viewer_process_multiplicate_vertices);
	Viewer::add_process('h', "Fill Holes", viewer_process_fill_holes);
	Viewer::add_process('-', "Make holes", viewer_process_delete_vertices);
	Viewer::add_process('d', "Delete non manifolds vertices", viewer_process_delete_non_manifold_vertices);
	Viewer::add_process('B', "Fill holes (biharmonic splines)", viewer_process_fill_holes_biharmonic_splines);
	Viewer::add_process('K', "Gaussian curvature", viewer_process_gaussian_curvature);
	Viewer::add_process('/', "Decimation", viewer_process_edge_collapse);
	
	Viewer::factor = shape_che->mean_edge();
	size_t g = shape_che->n_vertices() - shape_che->n_edges() + shape_che->n_faces();
	g = (g - 2)/(-2);
	debug(g)
	debug(Viewer::factor)

	size_t quality = 0;
	#pragma omp parallel for reduction(+: quality)
	for(index_t t = 0; t < shape_che->n_faces(); t++)
		quality += shape_che->pdetriq(t) > 0.6;
	
	debug(quality * 100.0 / shape_che->n_faces())

	Viewer::init(shape_che);
		
	delete shape_che;

	return 0;
}

void paint_holes_vertices()
{
	/*
	size_t nv = Viewer::mesh->n_vertices();
	
	Viewer::mesh.update();

	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
		if(v >= nv) Viewer::get_color(v) = .25;
		*/
}

void viewer_process_delete_non_manifold_vertices()
{
	debug_me(Processing:)

	debug_me(removing vertex);
	Viewer::mesh->remove_non_manifold_vertices();
	debug_me(removing vertex);
}

void viewer_process_delete_vertices()
{
	debug_me(Processing:)

	if(!Viewer::select_vertices.size()) return;
	debug_me(removing vertex);
	Viewer::mesh->remove_vertices(Viewer::select_vertices);
	Viewer::select_vertices.clear();
	debug_me(removing vertex);
}

void viewer_process_poisson(const index_t & k)
{
	size_t old_n_vertices = Viewer::mesh->n_vertices();
	delete [] fill_all_holes(Viewer::mesh);
	
	double time;
	TIC(time) poisson(Viewer::mesh, old_n_vertices, k); TOC(time)
	debug(time)
	
	paint_holes_vertices();
}

void viewer_process_poisson_laplacian_1()
{
	debug_me(Processing:)
	viewer_process_poisson(1);
}

void viewer_process_poisson_laplacian_2()
{
	debug_me(Processing:)
	viewer_process_poisson(2);
}

void viewer_process_poisson_laplacian_3()
{
	debug_me(Processing:)
	viewer_process_poisson(3);
}

void viewer_process_fill_holes()
{
	debug_me(Processing:)
	
	fill_all_holes(Viewer::mesh);
	
	paint_holes_vertices();
}

void viewer_process_noise()
{
	debug_me(Processing:)

	srand(time(NULL));

	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
	{
		distance_t r = distance_t( rand() % 1000 ) / 40000;
		int p = rand() % 5;
		Viewer::mesh->get_vertex(v) += (!p) * r * Viewer::mesh->normal(v);
	}

	Viewer::mesh.update_normals();
}

void viewer_process_thresold()
{
	debug_me(Processing:)

	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
		Viewer::get_color(v) = Viewer::get_color(v) > 0.5 ? 1 : 0.5;
}

void viewer_process_wks()
{
	debug_me(Processing:)
	float time;
	size_t K = 20, T = 10;
	
	sp_mat L, A;
	d_message(init laplacian...)
	TIC(time) laplacian(Viewer::mesh, L, A); TOC(time)
	debug(time)
	
	vec eigval;
	mat eigvec;
	
	d_message(init eigs...)
	TIC(time) eigs_laplacian(eigval, eigvec, Viewer::mesh, L, K); TOC(time)
	debug(time)

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
	{
		vec s(T, fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		Viewer::get_color(v) = norm(s);
		max_s = max(max_s, Viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
		Viewer::get_color(v) /= max_s;
}

void viewer_process_hks()
{
	debug_me(Processing:)
	float time;
	size_t K = 20, T = 100;
	
	sp_mat L, A;

	d_message(init laplacian...)
	TIC(time) laplacian(Viewer::mesh, L, A); TOC(time)
	debug(time)

	vec eigval;
	mat eigvec;
	
	d_message(init eigs...)
	TIC(time) eigs_laplacian(eigval, eigvec, Viewer::mesh, L, K); TOC(time)
	debug(time)

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
	{
		vec s(T, fill::zeros);
		for(index_t t = 0; t < T; t++)
		for(index_t k = 1; k < K; k++)
			s(t) += exp(-eigval(k) * t) * eigvec(v, k) * eigvec(v, k);

		Viewer::get_color(v) = norm(s);
		max_s = max(max_s, Viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
		Viewer::get_color(v) /= max_s;
}

void viewer_process_gps()
{
	debug_me(Processing:)
	float time;
	size_t K = 20;
	
	sp_mat L, A;

	d_message(init laplacian...)
	TIC(time) laplacian(Viewer::mesh, L, A); TOC(time)
	debug(time)

	vec eigval;
	mat eigvec;
	
	d_message(init eigs...)
	TIC(time) eigs_laplacian(eigval, eigvec, Viewer::mesh, L, K); TOC(time)
	debug(time)

	eigvec.col(0).zeros();
	for(index_t i = 1; i < K; i++)
		eigvec.col(i) /= sqrt(eigval(i));

	mat data = eigvec.t();
	mat means;

	distance_t max_s = 0;
	#pragma omp parallel for reduction(max: max_s)
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
	{
		Viewer::get_color(v) = norm(eigvec.row(v));
			max_s = max(max_s, Viewer::get_color(v));
	}

	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
		Viewer::get_color(v) /= max_s;
}

void viewer_process_denoising()
{
	debug_me(Processing:)
	
	size_t K, m, M;
	distance_t f;
	cin >> K >> m >> M >> f;

	mesh_denoising(Viewer::mesh, Viewer::select_vertices, K, m, M, f);

	Viewer::mesh.update_normals();
}

void viewer_process_super_resolution()
{
	debug_me(Processing:)
	
	size_t K, m, M;
	distance_t f;
	cin >> K >> m >> M >> f;

	mesh_super_resolution(Viewer::mesh, Viewer::select_vertices, K, m, M,f);
}	

void viewer_process_inpaiting()
{
	debug_me(Processing:)
	
	size_t K, m, M;
	distance_t f;
	cin >> K >> m >> M >> f;

	mesh_inpaiting(Viewer::mesh, Viewer::select_vertices, K, m, M,f);
	
	paint_holes_vertices();
}


void viewer_process_iterative_inpaiting()
{
	debug_me(Processing:)
	
	size_t K, m, M;
	distance_t f;
	cin >> K >> m >> M >> f;

	size_t n_v = Viewer::mesh->n_vertices();
	
	mesh_iterative_inpaiting(Viewer::mesh, Viewer::select_vertices, K, m, M,f);

	paint_holes_vertices();
}

void viewer_process_multiplicate_vertices()
{
	debug_me(Processing:)

	Viewer::mesh->multiplicate_vertices();
	Viewer::factor = Viewer::mesh->mean_edge();
	
	debug(Viewer::mesh->n_vertices())
}

void viewer_sort_by_rings()
{
	debug_me(Processing:)
	
	index_t * rings = new index_t[Viewer::mesh->n_vertices()];
	index_t * sorted = new index_t[Viewer::mesh->n_vertices()];
	vector<index_t> limites;
	Viewer::mesh->sort_by_rings(rings, sorted, limites, Viewer::select_vertices);
	
	size_t k = limites.size() - 1;

	for(index_t v = 0; v < Viewer::mesh->n_vertices(); v++)
	{
		if(rings[v] < k) //Viewer::select_vertices.push_back(v);
		Viewer::get_color(v) = distance_t(rings[v]) / (limites.size() - 1);
	}
	debug(k)

	delete [] rings;
	delete [] sorted;
}

void viewer_process_fastmarching_cpu()
{
	if(!Viewer::select_vertices.size()) return;
	
	debug_me(Processing:)
	
	index_t * rings = new index_t[Viewer::mesh->n_vertices()];
	index_t * sorted = new index_t[Viewer::mesh->n_vertices()];
	vector<index_t> limites;
	Viewer::mesh->sort_by_rings(rings, sorted, limites, Viewer::select_vertices);
	
	float time;
	distance_t * distances = parallel_fastmarching(Viewer::mesh, Viewer::select_vertices.data(), Viewer::select_vertices.size(), time, limites, sorted, true, NULL, false);	
	//distance_t * distances = parallel_fastmarching(Viewer::mesh->filename().c_str(), Viewer::select_vertices.data(), Viewer::select_vertices.size(), time, 9, true, true);

	debug(time)

	Viewer::mesh.update_colors(distances);

	delete [] distances;
	delete [] rings;
	delete [] sorted;
}
void viewer_process_fastmarching_gpu()
{
	if(!Viewer::select_vertices.size()) return;
	
	debug_me(Processing:)
	
	index_t * rings = new index_t[Viewer::mesh->n_vertices()];
	index_t * sorted = new index_t[Viewer::mesh->n_vertices()];
	vector<index_t> limites;
	Viewer::mesh->sort_by_rings(rings, sorted, limites, Viewer::select_vertices);
	
	float time;
	distance_t * distances = parallel_fastmarching(Viewer::mesh, Viewer::select_vertices.data(), Viewer::select_vertices.size(), time, limites, sorted, true);	
	//distance_t * distances = parallel_fastmarching(Viewer::mesh->filename().c_str(), Viewer::select_vertices.data(), Viewer::select_vertices.size(), time, 9, true, true);

	debug(time)

	Viewer::mesh.update_colors(distances);
	
	delete [] distances;
	delete [] rings;
	delete [] sorted;
}

void viewer_process_voronoi()
{
	if(!Viewer::select_vertices.size()) return;
	
	debug_me(Processing:)
	
	index_t * rings = new index_t[Viewer::mesh->n_vertices()];
	index_t * sorted = new index_t[Viewer::mesh->n_vertices()];
	vector<index_t> limites;
	Viewer::mesh->sort_by_rings(rings, sorted, limites, Viewer::select_vertices);
	
	float time;
	index_t * clusters = new index_t[Viewer::mesh->n_vertices()];
	memset(clusters, 255, sizeof(index_t) * Viewer::mesh->n_vertices());
	
	distance_t * distances = parallel_fastmarching(Viewer::mesh, Viewer::select_vertices.data(), Viewer::select_vertices.size(), time, limites, sorted, true, clusters);	
	
	debug(time)

	#pragma omp parallel for
	for(index_t i = 0; i < Viewer::mesh->n_vertices(); i++)
	{
		Viewer::get_color(i) = clusters[i];
		Viewer::get_color(i) /= Viewer::select_vertices.size() + 1;
	}

	delete [] distances;
	delete [] rings;
	delete [] sorted;
	delete [] clusters;
}

void viewer_process_farthest_point_sampling_radio()
{
	debug_me(Processing:)
	
	if(!Viewer::select_vertices.size())
		Viewer::select_vertices.push_back(0);

	distance_t radio;
	cin>>radio;

	float time, time_g;
		
	TIC(time)
	radio = farthest_point_sampling_gpu(Viewer::select_vertices, time_g, Viewer::mesh, Viewer::mesh->n_vertices(), radio);
	TOC(time)
	
	debug(radio)
	debug(Viewer::select_vertices.size())
	debug(time)
	debug(time_g)
}

void viewer_process_farthest_point_sampling()
{
	debug_me(Processing:)
	
	if(!Viewer::select_vertices.size())
		Viewer::select_vertices.push_back(0);

	index_t n;
	cin>>n;

	float time, time_g;
		
	TIC(time)
	distance_t radio = farthest_point_sampling_gpu(Viewer::select_vertices, time_g, Viewer::mesh, n);
	TOC(time)
	
	debug(radio)
	debug(Viewer::select_vertices.size())
	debug(time)
	debug(time_g)
}

void viewer_process_fairing_spectral()
{
	debug_me(Processing:)

	fairing * fair = new fairing_spectral(50);
	fair->run(Viewer::mesh);
	Viewer::mesh->set_vertices(fair->get_postions());
	Viewer::mesh->normalize();
	delete fair;
	
	Viewer::mesh.update_normals();
}

void viewer_process_fairing_taubin()
{
	debug_me(Processing:)

	fairing * fair = new fairing_taubin;
	fair->run(Viewer::mesh);
	Viewer::mesh->set_vertices(fair->get_postions());
	Viewer::mesh->normalize();
	delete fair;

	Viewer::mesh.update_normals();
}

void viewer_process_geodesics()
{
	debug_me(Processing:)
	
	double start_omp = omp_get_wtime();
	geodesics geodesic(Viewer::mesh, Viewer::select_vertices);
	double time = omp_get_wtime() - start_omp;
	debug(time)

//	geodesic.path_to(Viewer::other_vertices, Viewer::mesh, 0);
	
	geodesic.normalize();
	Viewer::mesh.update_colors(geodesic.distances);
/* heat map one color
	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh.n_vertices(); v++)
		Viewer::get_color(v) = Viewer::get_color(v) * 0.5;	
*/		
		
}

void viewer_process_fastmarching()
{
	debug_me(Processing:)
	off shape(Viewer::mesh->filename());

	#pragma omp parallel for
	for(index_t v = 0; v < shape.get_nvertices(); v++)
		shape(v) = Viewer::mesh->get_vertex(v);

	double time;
	TIC(time)	
	fastmarching fm(shape, Viewer::select_vertices, INFINITY, true, false);
	TOC(time)
	debug(time)

	Viewer::mesh.update_colors(fm.distances);
}

void viewer_process_fill_holes_biharmonic_splines()
{
	debug_me(Processing:)
	
	size_t old_n_vertices, n_vertices = Viewer::mesh.n_vertices();
	size_t n_holes = Viewer::mesh->n_borders();

	vector<index_t> * border_vertices;
	che ** holes;
	tie(border_vertices, holes) = fill_all_holes_meshes(Viewer::mesh);
	if(!holes) return;
	
	index_t k = 2;

	for(index_t h = 0; h < n_holes; h++)
		if(holes[h])
		{
			old_n_vertices = n_vertices;
			biharmonic_interp_2(Viewer::mesh, old_n_vertices, n_vertices += holes[h]->n_vertices() - border_vertices[h].size(), border_vertices[h], k);
			delete holes[h];
		}
	
	delete [] holes;
	delete [] border_vertices;
	paint_holes_vertices();
}

void viewer_process_gaussian_curvature()
{
	debug_me(Processing:)
	
	vertex_t g, g_max = -INFINITY, g_min = INFINITY;
	vertex a, b;

#ifdef SINGLE_P
	fvec gv(Viewer::mesh.n_vertices());
#else
	vec gv(Viewer::mesh.n_vertices());
#endif

	#pragma omp parallel for private(g, a, b) reduction(max: g_max) reduction(min: g_min)
	for(index_t v = 0; v < Viewer::mesh.n_vertices(); v++)
	{
		g = 0;
		for_star(he, Viewer::mesh, v)
		{
			a = Viewer::mesh->gt_vt(next(he)) - Viewer::mesh->gt(v);
			b = Viewer::mesh->gt_vt(prev(he)) - Viewer::mesh->gt(v);
			g += acos((a,b) / (*a * *b));
		}
		gv(v) = (2 * M_PI - g) / Viewer::mesh->area_vertex(v);
		g_max = max(g_max, gv(v));
		g_min = min(g_min, gv(v));
	}

	g = g_max - g_min;
	
	#pragma omp parallel for
	for(index_t v = 0; v < Viewer::mesh.n_vertices(); v++)
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
	for(index_t v = 0; v < Viewer::mesh.n_vertices(); v++)
		Viewer::get_color(v) = f(gv(v));
}

void viewer_process_edge_collapse()
{
	debug_me(Processing:)
	
	double time;

	TIC(time) decimation sampling(Viewer::mesh); TOC(time)
	debug(time)
}
