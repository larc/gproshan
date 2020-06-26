#include "inpainting.h"

#include "che_off.h"

#include <cassert>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <queue>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


inpainting::inpainting(che *const & _mesh, basis *const & _phi_basis, const params & p): dictionary(_mesh, _phi_basis, p.m, p.M, p.f, p.learn, p.plot)
{
	delta = p.delta;
	sum_thres = p.sum_thres;
	avg_p = p.avg_p;
	percent = p.percentage;
	area_thres = p.area_thres;

	M = mesh->n_vertices() / avg_p;
	mask = new bool[mesh->n_vertices()];
	memset(mask, 0, sizeof(bool) * mesh->n_vertices());
}


void inpainting::load_mask()
{
	//string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(avg_p) + '_' + to_string(percent) + '_' + to_string(radio) + ".msk");
	string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(percent) + ".msk");
	arma::uvec V;
	gproshan_log(loading radial mask);
	

	if(V.load(f_mask))
	{
		#pragma omp for 
		for(index_t i = 0; i < mesh->n_vertices(); i++)
		{
			mask[i] = V(i);
		}
	}
	else
	{
		V.zeros(mesh->n_vertices());
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0, mesh->n_vertices()-1);
		size_t percentage = mesh->n_vertices() - ceil(mesh->n_vertices() * (percent/ 100.0)) ;

		size_t k = 0;
		size_t rn = 0;
		while (k < percentage)
		{
		
			rn = distribution(generator);
			if(!mask[rn])
			{
				mask[rn] = 1;
				V(rn) = 1;
				k++;
			}
				
		}
		V.save(f_mask);
	}
	
}

void inpainting::load_mask(const std::vector<index_t> * vertices, const index_t * clusters)
{
	string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(avg_p) + '_' + to_string(percent) + ".msk");
	arma::uvec V;

	if(V.load(f_mask))
	{
		#pragma omp for 
		for(index_t i = 0; i < mesh->n_vertices(); i++)
		{
			mask[i] = V(i);
		}

	}
	else
	{

		V.zeros(mesh->n_vertices());
		size_t * percentages_size = new size_t[M];
		bool cover_cluster[M];

			
		// create initial desired percentage sizes
		#pragma omp for 
		for(index_t s = 0; s < M; s++)
		{
		
			percentages_size[s] = ceil(vertices[s].size() * (percent/ 100.0)) ;
			cover_cluster[s] = 0;
			
		}
		//Generate random mask according to a percentage of patches capacity
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0,n_vertices-1);

		size_t k = 0;
		size_t rn=0;
		

		while( k < M )
		{	
			//gproshan_debug_var(clusters[rn]);

			rn = distribution(generator);
			if(!mask[rn] && percentages_size[clusters[rn]] > 0)
			{

				mask[rn] = 1;
				V(rn) = 1;
				percentages_size[clusters[rn]]--;
							
			}
			if( !cover_cluster[clusters[rn]] && percentages_size[clusters[rn]] == 0)
			{
				cover_cluster[clusters[rn]] = 1; // It is finished
				k++;
			}	
				
		}
	
	
		V.save(f_mask);

	}
	
}

void inpainting::load_sampling(bool save_all)
{
	// saving sampling
	//double delta = PI/6;
	//double sum_thres = 0.008; // best arma 0.0001 the worst with 0.001, now with 0.0005 lets see tomorrow
	//double sum_thres = PI;
	string f_sampl = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".rsampl");

	arma::mat S;

	vector<index_t> all_sorted_features;
	vector<index_t> seeds;
	size_t featsize;
	load_features(all_sorted_features, featsize);

	gproshan_debug_var(all_sorted_features.size());
	string f_points = tmp_file_path(mesh->name_size() + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".points");

//	vector<index_t> features(all_sorted_features.begin(), all_sorted_features.begin() + featsize );
//	gproshan_debug_var(features.size());
//	geodesics geo(mesh, features , geodesics::FM, NULL, false, mesh->n_vertices());
//	index_t * indexes = new index_t[geo.n_sorted_index()];
//	geo.copy_sorted_index(indexes, geo.n_sorted_index());
	size_t count = 0;
//	real_t max_radio = geo[indexes[mesh->n_vertices()-1]] ;
	
	//radio *= 1.1;
	real_t area_mesh = mesh->area_surface();
	real_t max_radio = 0.05 * area_mesh;
	gproshan_debug_var(max_radio);

	if(S.load(f_sampl))
	{
		gproshan_debug(already done);
		size_t n_seeds = S.n_rows;
		real_t euc_radio, geo_radio;
		for(index_t i = 0; i < n_seeds; i++)
		{
			patch p;
		//	gproshan_debug_var(i);
			//p.recover_radial_disjoint( mesh, S(i,1), S(i,0) );
		//	p.init_radial_disjoint(mesh, S(i,1), S(i,0), euc_radio, geo_radio, delta, sum_thres, area_thres);
			patches.push_back(p); 
		}
		
		M = n_seeds;
		string f_points = tmp_file_path(mesh->name_size() + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".points");
		gproshan_debug_var(n_seeds);
		a_vec outlv(n_seeds);
		gproshan_debug(restored);
		for(index_t i = 0; i < n_seeds; i++)
		{
			outlv(i) = S(i,0);
		}
		outlv.save(f_points);
		gproshan_debug(restored);

		return;
	}

	//ELSE IF !S.load(f_sample)
	
	gproshan_debug(SAMPLING);
	
	// compute features will be seeds
	patches_map.resize(mesh->n_vertices());

	//Coverage of the points 
	bool covered[mesh->n_vertices()];

	#pragma omp for 
	for(index_t i = 0; i < mesh->n_vertices(); i++)
		covered[i] = 0;
	
	real_t euc_radio;
	real_t geo_radio;
	vector<real_t> radios;
	vector<real_t> geo_radios;
	size_t count_cov = 0;
	size_t count_cov_patch = 0;

	bool * invalid_seed = new bool[mesh->n_vertices()];
	memset(invalid_seed, 0, mesh->n_vertices() * sizeof(bool));

	for(const index_t & vsf: all_sorted_features)
	{
		if(invalid_seed[vsf]) continue;

		patch p;
		// increasing a bit the radio
		p.init_radial_disjoint(euc_radio, geo_radio, mesh, vsf, max_radio, delta, sum_thres, area_thres, area_mesh);

		//gproshan_debug_var(p.vertices.size());
		count_cov_patch = 0;
		if(p.vertices.size() >= 7 )
		{
			for(const index_t & v: p.vertices)
				if(!covered[v]) count_cov_patch++;

			count_cov += count_cov_patch;
			if(count_cov_patch > 0)
			{
				//gproshan_debug_var(p.vertices.size());
				patches.push_back(p);
				seeds.push_back(vsf);
				radios.push_back(euc_radio);
				geo_radios.push_back(geo_radio);
				count += p.vertices.size();

				for(const index_t & v: p.vertices)
				{
					covered[v] = 1;

					if(!invalid_seed[v])
					{
						const vertex & va = mesh->get_vertex(vsf);
						const vertex & vb = mesh->get_vertex(v);

						invalid_seed[v] = *(va - vb) < 0.4 * euc_radio;
					}
				}
			}
		}
	}						

	delete [] invalid_seed;

	vector<index_t> outliers;
	gproshan_debug_var(sum_thres);
	gproshan_debug_var(count);
	gproshan_debug_var(count_cov);
	gproshan_debug_var(seeds.size());
	M = seeds.size();

	gproshan_debug(SAMPLING);
	///////////////////////////////////////

#ifndef NDEBUG
	gproshan_debug_var(M);
	index_t tmp;
	bool remark[mesh->n_vertices()] = {};
	
	gproshan_debug_var(outliers.size());
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		if(!covered[i] )
		{
			outliers.push_back(i);
		}
	}
	a_vec outlv(outliers.size() );
	//gproshan_debug_var(seeds.size());
	for(index_t i = 0; i < outliers.size(); i++)
	{
		//outlv(i) = seeds[i];
		outlv(i) = outliers[i];
		//gproshan_debug_var(seeds[i]);
	}
		
	gproshan_debug_var(f_points);
	outlv.save(f_points);
	S.resize(seeds.size(),2);
	for(index_t i = 0; i < seeds.size(); i++)
	{
		S(i,0) = seeds[i];
		S(i,1) = geo_radios[i];
	}
	S.save(f_sampl);
#endif // NDEBUG
}

void inpainting::init_radial_feature_patches()
{
	load_sampling(false);

	load_mask();

	//Initializing patches
	gproshan_log(initializing patches);
	n_vertices = mesh->n_vertices();
	patches.resize(M);
	patches_map.resize(n_vertices);
	//initializing patch
	gproshan_debug_var(M);
	#ifndef NDEBUG
		size_t patch_avg_size = 0;
		size_t patch_min_size = NIL;
		size_t patch_max_size = 0;

		#pragma omp parallel for reduction(+: patch_avg_size)
		for(index_t s = 0; s < M; s++)
			patch_avg_size += patches[s].vertices.size();
		#pragma omp parallel for reduction(min: patch_min_size)
		for(index_t s = 0; s < M; s++)
			patch_min_size = min(patches[s].vertices.size(), patch_min_size);
		#pragma omp parallel for reduction(max: patch_max_size)
		for(index_t s = 0; s < M; s++)
			patch_max_size = max(patches[s].vertices.size(), patch_max_size);

		patch_avg_size /= M;
		gproshan_debug_var(patch_avg_size);
		gproshan_debug_var(patch_min_size);
		gproshan_debug_var(patch_max_size);
	#endif

	bool * pmask = mask;
	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz_disjoint(mesh, dist, M, patches_map, s ,[&pmask](const index_t & i) -> bool { return pmask[i]; } );

	gproshan_debug(passed);

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.scale_xyz(phi_basis->get_radio());
		// adding add extra function
		p.add_extra_xyz_disjoint(mesh,patches_map, s);
		p.compute_avg_distance(mesh, patches_map, s);
		p.phi.set_size(p.xyz.n_cols, phi_basis->get_dim());
		phi_basis->discrete(p.phi, p.xyz);
	}	 

	bool save_all = true;
	if(save_all)
	{
		arma::mat AS;
		AS.resize(M,13);
		for(index_t i = 0; i < M; i++)
		{
			patch & p = patches[i];
			
			AS(i,0) = p.x(0);
			AS(i,1) = p.x(1);
			AS(i,2) = p.x(2);
			AS(i,3) = p.radio;
			AS(i,4) = p.T(0,0);
			AS(i,5) = p.T(1,0);
			AS(i,6) = p.T(2,0);
			AS(i,7) = p.T(0,1);
			AS(i,8) = p.T(1,1);
			AS(i,9) = p.T(2,1);
			AS(i,10) = p.T(0,2);
			AS(i,11) = p.T(1,2);
			AS(i,12) = p.T(2,2);

		}
		string f_samplall = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".smp");
		AS.save(f_samplall);
	}
	gproshan_log(radial patches are ready);
}

void inpainting::init_voronoi_patches()
{
	/////
	std::vector<index_t> vertices[M];
	//index_t * clusters_tmp = init_voronoi_sampling(vertices);

	//////
	gproshan_log_var(M);

		//FPS samplif_dictng with desired number of sources
		TIC(d_time) init_sampling(); TOC(d_time)
		gproshan_debug_var(d_time);
	
		geodesics::params params;
		params.cluster = 1;


		// creating disjoint clusters with geodesics aka voronoi
	#ifdef GPROSHAN_CUDA
		params.alg = geodesics::PTP_GPU;
	#endif
		
		geodesics ptp(mesh, sampling, params);
		TOC(d_time)
		gproshan_log_var(d_time);


		//saving first vertex aka seed vertices
		#pragma omp for 
		for(index_t s = 0; s < M; s++)
		{
			vertices[s].push_back(sample(s));
		}

		for(index_t i = 0; i < mesh->n_vertices(); i++)
			{		
				ptp.clusters[i]--;
				if(sample(ptp.clusters[i]) != i)
					vertices[ptp.clusters[i]].push_back(i) ;
			}
	
	
	load_mask(vertices, ptp.clusters);

	//Initializing patches
	gproshan_log(initializing patches);

	patches.resize(M);
	patches_map.resize(n_vertices);
	//initializing patch
	#pragma omp parallel
	{
		index_t * toplevel = new index_t[mesh->n_vertices()];

		#pragma omp for 
		for(index_t s = 0; s < M; s++)
		{
			patches[s].init_disjoint(mesh, sample(s), dictionary::T, vertices[s], toplevel);
			
		}		
		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < M; s++)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < M; s++)
				patch_min_size = min(patches[s].vertices.size(), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < M; s++)
				patch_max_size = max(patches[s].vertices.size(), patch_max_size);

			patch_avg_size /= M;
			//gproshan_debug_var(patch_avg_size);
			//gproshan_debug_var(patch_min_size);
			//gproshan_debug_var(patch_max_size);
		#endif
	}

	bool * pmask = mask;
	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz_disjoint(mesh, dist, M, patches_map, s ,[&pmask](const index_t & i) -> bool { return pmask[i]; } );

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->get_dim());
		phi_basis->discrete(p.phi, p.xyz);
		p.compute_avg_distance(mesh, patches_map, s);

	} 

	gproshan_log(our patches are ready);
}

real_t inpainting::execute()
{

	TIC(d_time) init_radial_feature_patches(); TOC(d_time)
	gproshan_debug_var(d_time);
//	L = 15;

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) learning(); TOC(d_time)
	gproshan_debug_var(d_time);

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);
	string f_alpha = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".alpha");
	save_alpha(f_alpha);

	
	draw_patches(295);
	draw_patches(384);
	draw_patches(319);
	draw_patches(312);
	//patches_map.resize(n_vertices);
	for(index_t i = 0; i < n_vertices; i++)
	{
		patches_map[i].clear();
	}

	
	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz(mesh, patches_map, s, 0);



	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.scale_xyz(phi_basis->get_radio());
		p.phi.set_size(p.xyz.n_cols, phi_basis->get_dim());
		phi_basis->discrete(p.phi, p.xyz);

	}

	bool *pmask;


	/*draw_patches(295);
	draw_patches(384);
	draw_patches(319);
	draw_patches(312);*/
	//draw_patches(76);
	//draw_patches(1486);
	
	//draw_patches(400);
	//draw_patches(500);
	//draw_patches(56);
	//phi_basis->plot_basis();
	//gproshan_debug_var(alpha.col(463));
	

	TIC(d_time) mesh_reconstruction([&pmask](const index_t & i) -> bool { return pmask[i]; }); TOC(d_time)
	gproshan_debug_var(d_time);
	arma::uvec non_zero = find( abs(alpha) > 0.00001);
	gproshan_debug_var(non_zero.size());
	real_t ratio = (M * 13.0 + non_zero.size()) / (3 * mesh->n_vertices());
	gproshan_debug_var(ratio);

}

che * inpainting::point_cloud_reconstruction(real_t per, real_t fr)
{
	arma::mat S;
	arma::mat T(3,3);
	arma::mat alpha;


	string f_smp = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres) + '_' + to_string(area_thres) + ".smp");
	string f_alpha = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres) +'_' + to_string(area_thres) + ".alpha");

	S.load(f_smp);
	alpha.load(f_alpha);
	gproshan_debug_var(S.n_rows);
	M = S.n_rows;
	real_t radio, max_radio = -1;
	patches.resize(M);
	vertex c;

	for(index_t i = 0; i < M; i++)
		if( S(i,3) > max_radio) max_radio = S(i,3); 

	size_t total_points = 0;
	A.eye(phi_basis->get_dim(), m);
	

	for(index_t i = 0; i < M; i++)
	{
		c.x = S(i,0);
		c.y = S(i,1);
		c.z = S(i,2);
		radio = S(i,3);
		T(0,0) = S(i,4);
		T(1,0) = S(i,5);
		T(2,0) = S(i,6);
		T(0,1) = S(i,7);
		T(1,1) = S(i,8);
		T(2,1) = S(i,9);
		T(0,2) = S(i,10);
		T(1,2) = S(i,11);
		T(2,2) = S(i,12);

		patches[i].init_random(c, T, radio, max_radio, per, fr);
		total_points += patches[i].vertices.size();
		patches[i].phi.set_size(patches[i].vertices.size(), phi_basis->get_dim());
		phi_basis->discrete(patches[i].phi, patches[i].xyz);
		

		a_vec x = patches[i].phi * A * alpha.col(i);
		patches[i].xyz.row(2) = x.t();

		for(index_t j = 0; j < patches[i].vertices.size(); j++)
			if (patches[i].xyz(2, j) > 2)
				gproshan_debug_var( i);
		

	}

	gproshan_debug_var(total_points);

	for(index_t i = 0; i < M; i++)
	{

		patches[i].iscale_xyz(phi_basis->get_radio());
		patches[i].itransform();
	}

	vector<vertex> point_cloud;
	point_cloud.reserve(total_points);

	for(index_t i = 0; i < M; i++)
	for(index_t j = 0; j < patches[i].vertices.size(); j++)
		point_cloud.push_back({	patches[i].xyz(0, j), 
								patches[i].xyz(1, j),
								patches[i].xyz(2, j)
								});

	che * nmesh = new che(point_cloud.data(), point_cloud.size(), nullptr, 0);

	gproshan_debug_var(sum_thres);
	string f_pc = tmp_file_path(mesh->name_size() + '_' + to_string(delta) + '_' + to_string(sum_thres)+ '_' + to_string(area_thres) 
	+ '_' + to_string(per) + '_' + to_string(fr) + "_pc");

	che_off::write_file(nmesh, f_pc);
	gproshan_debug(Done!);
	
	return nmesh;
}


real_t inpainting::execute_tmp()
{
	// fill holes
	size_t threshold = mesh->n_vertices();
	delete [] fill_all_holes(mesh);
	TIC(d_time) poisson(mesh, threshold, 2); TOC(d_time)
	gproshan_debug_var(d_time);

	// remove possible non manifold vertices
	mesh->remove_non_manifold_vertices();

	// sampling including new vertices
	TIC(d_time) init_sampling(); TOC(d_time)
	gproshan_debug_var(d_time);

	// initializing patches with threshold
	TIC(d_time) init_patches(1, [&threshold](const index_t & i) -> bool { return i < threshold; }); TOC(d_time)
	gproshan_debug_var(d_time);

	// learning only from valid patches
	TIC(d_time) learning(); TOC(d_time)
	gproshan_debug_var(d_time);

	// including vertices out of threshold
	TIC(d_time) init_patches(0); TOC(d_time)
	gproshan_debug_var(d_time);

	// Update new alphas, propagating the info towards the center
	TIC(d_time) update_alphas(alpha, threshold); TOC(d_time)
	gproshan_debug_var(d_time);

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	gproshan_debug_var(d_time);
}


} // namespace gproshan::mdict

