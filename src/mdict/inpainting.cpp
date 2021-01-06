#include "mdict/inpainting.h"

#include "mesh/che_off.h"
#include "geodesics/geodesics.h"

#include <cassert>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <queue>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


inpainting::inpainting(che *const & _mesh, basis *const & _phi_basis, const params & p): msparse_coding(_mesh, _phi_basis, p)
{
	m_params.n_patches = mesh->n_vertices() / m_params.avg_p;
	
	key_name = mesh->name_size() + '_' + to_string(m_params.delta) + '_' + to_string(m_params.sum_thres) + '_' + to_string(m_params.area_thres);

	mask = new bool[mesh->n_vertices()];
	memset(mask, 0, sizeof(bool) * mesh->n_vertices());
}

inpainting::~inpainting()
{
	delete [] mask;
}

inpainting::operator const std::string & () const
{
	return key_name;
}

void inpainting::load_mask()
{
	//string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(avg_p) + '_' + to_string(percent) + '_' + to_string(radio) + ".msk");
	string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(m_params.percent) + ".msk");
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
		size_t percentage = mesh->n_vertices() - ceil(mesh->n_vertices() * (m_params.percent / 100.0)) ;

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
	string f_mask = tmp_file_path(mesh->name_size() + '_' + to_string(m_params.avg_p) + '_' + to_string(m_params.percent) + ".msk");
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
		size_t * percentages_size = new size_t[m_params.n_patches];
		bool cover_cluster[m_params.n_patches];

			
		// create initial desired percentage sizes
		#pragma omp for 
		for(index_t s = 0; s < m_params.n_patches; s++)
		{
			percentages_size[s] = ceil(vertices[s].size() * (m_params.percent / 100.0)) ;
			cover_cluster[s] = 0;
		}

		//Generate random mask according to a percentage of patches capacity
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0,n_vertices-1);

		size_t k = 0;
		size_t rn = 0;
		
		while(k < m_params.n_patches)
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
	size_t featsize;
	vector<index_t> all_sorted_features;
	vector<index_t> seeds;
	load_features(all_sorted_features, featsize);

	gproshan_debug_var(all_sorted_features.size());

	size_t count = 0;
	real_t area_mesh = mesh->area_surface();

	a_vec S;
	if(S.load(tmp_file_path(key_name + ".rsampl")))
	{
		gproshan_debug(loading sampling);

		size_t n_seeds = S.n_rows;
		real_t euc_radio, geo_radio;
		for(index_t i = 0; i < n_seeds; i++)
		{
			patch p;
			p.init_radial_disjoint(euc_radio, geo_radio, mesh, S(i), m_params.delta, m_params.sum_thres, m_params.area_thres, area_mesh);
			patches.push_back(move(p)); 
		}
		
		m_params.n_patches = n_seeds;
		
		return;
	}

	//ELSE IF !S.load(f_sample)
	
	gproshan_debug(computing sampling);
	
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
		p.init_radial_disjoint(euc_radio, geo_radio, mesh, vsf, m_params.delta, m_params.sum_thres, m_params.area_thres, area_mesh);

		count_cov_patch = 0;
		if(p.vertices.size() >= 7 )
		{
			for(const index_t & v: p.vertices)
				if(!covered[v]) count_cov_patch++;

			count_cov += count_cov_patch;
			if(count_cov_patch > 0)
			{
				patches.push_back(move(p));
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

						invalid_seed[v] = *(va - vb) < 0.8 * euc_radio;
					}
				}
			}
		}
	}						

	delete [] invalid_seed;

	m_params.n_patches = seeds.size();
	
	gproshan_log(saving sampling);

	S.resize(seeds.size());

	#pragma omp parallel for
	for(index_t i = 0; i < seeds.size(); i++)
		S(i) = seeds[i];
	
	S.save(tmp_file_path(key_name + ".rsampl"));

	gproshan_debug_var(m_params.sum_thres);
	gproshan_debug_var(count);
	gproshan_debug_var(count_cov);
	gproshan_debug_var(seeds.size());
	gproshan_debug_var(m_params.n_patches);


#ifndef NDEBUG
	vector<index_t> outliers;
	

	for(index_t i = 0; i < mesh->n_vertices(); i++)
		if(!covered[i])
			outliers.push_back(i);

	a_vec outlv(outliers.size());
	for(index_t i = 0; i < outliers.size(); i++)
		outlv(i) = outliers[i];
		
	outlv.save(tmp_file_path(key_name + ".outlv"));
#endif // NDEBUG
}

void inpainting::init_radial_feature_patches()
{
	load_sampling(false);
	load_mask();

	#ifndef NDEBUG
		size_t patch_avg_size = 0;
		size_t patch_min_size = NIL;
		size_t patch_max_size = 0;

		#pragma omp parallel for reduction(+: patch_avg_size)
		for(index_t s = 0; s < m_params.n_patches; s++)
			patch_avg_size += patches[s].vertices.size();
		#pragma omp parallel for reduction(min: patch_min_size)
		for(index_t s = 0; s < m_params.n_patches; s++)
			patch_min_size = min(patches[s].vertices.size(), patch_min_size);
		#pragma omp parallel for reduction(max: patch_max_size)
		for(index_t s = 0; s < m_params.n_patches; s++)
			patch_max_size = max(patches[s].vertices.size(), patch_max_size);

		patch_avg_size /= m_params.n_patches;
		gproshan_debug_var(patch_avg_size);
		gproshan_debug_var(patch_min_size);
		gproshan_debug_var(patch_max_size);
	#endif

	//Initializing patches
	gproshan_log(initializing patches);
	gproshan_debug_var(m_params.n_patches);
	
	n_vertices = mesh->n_vertices();
	patches.resize(m_params.n_patches);
	patches_map.resize(n_vertices);
	
	bool * pmask = mask;
	
	for(index_t s = 0; s < m_params.n_patches; s++)
		patches[s].reset_xyz_disjoint(mesh, dist, m_params.n_patches, patches_map, s, [&pmask](const index_t & i) -> bool { return pmask[i]; });

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.scale_xyz(phi_basis->radio());
		p.add_extra_xyz_disjoint(mesh,patches_map, s);
		p.compute_avg_distance(mesh, patches_map, s);
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
	}	 

	bool save_all = true;
	if(save_all)
	{
		arma::mat AS;
		AS.resize(m_params.n_patches, 13);
		for(index_t i = 0; i < m_params.n_patches; i++)
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
		
		AS.save(tmp_file_path(key_name + ".smp"));
	}
	
	gproshan_log(radial patches are ready);
}

void inpainting::init_voronoi_patches()
{
	/////
	std::vector<index_t> vertices[m_params.n_patches];
	//index_t * clusters_tmp = init_voronoi_sampling(vertices);

	//////
	gproshan_log_var(m_params.n_patches);

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
		for(index_t s = 0; s < m_params.n_patches; s++)
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

	patches.resize(m_params.n_patches);
	patches_map.resize(n_vertices);
	//initializing patch
	#pragma omp parallel
	{
		index_t * toplevel = new index_t[mesh->n_vertices()];

		#pragma omp for 
		for(index_t s = 0; s < m_params.n_patches; s++)
		{
			patches[s].init_disjoint(mesh, sample(s), msparse_coding::T, vertices[s], toplevel);
			
		}		
		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_min_size = min(patches[s].vertices.size(), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < m_params.n_patches; s++)
				patch_max_size = max(patches[s].vertices.size(), patch_max_size);

			patch_avg_size /= m_params.n_patches;
			//gproshan_debug_var(patch_avg_size);
			//gproshan_debug_var(patch_min_size);
			//gproshan_debug_var(patch_max_size);
		#endif
	}

	bool * pmask = mask;
	for(index_t s = 0; s < m_params.n_patches; s++)
		patches[s].reset_xyz_disjoint(mesh, dist, m_params.n_patches, patches_map, s ,[&pmask](const index_t & i) -> bool { return pmask[i]; } );

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
		p.compute_avg_distance(mesh, patches_map, s);

	} 

	gproshan_log(our patches are ready);
}

real_t inpainting::execute()
{

	TIC(d_time) init_radial_feature_patches(); TOC(d_time)
	gproshan_debug_var(d_time);
	//L = 10;

	// sparse coding and reconstruction with all patches
	//TIC(d_time) sparse_coding(); TOC(d_time)
	//gproshan_debug_var(d_time);

	TIC(d_time) learning(); TOC(d_time)
	gproshan_debug_var(d_time);

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);

	save_alpha(tmp_file_path(key_name + ".alpha"));

	draw_patches(295);
	draw_patches(312);

	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices; i++)
		patches_map[i].clear();
	
	for(index_t s = 0; s < m_params.n_patches; s++)
		patches[s].reset_xyz(mesh, patches_map, s, 0);

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.scale_xyz(phi_basis->radio());
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
	}

	bool * pmask = mask;

	TIC(d_time)
	real_t max_error = mesh_reconstruction([&pmask](const index_t & i) -> bool { return pmask[i]; }); 
	TOC(d_time)
	gproshan_debug_var(d_time);
	
	arma::uvec non_zero = find( abs(alpha) > 0.00001);
	gproshan_debug_var(non_zero.size());
	real_t ratio = (m_params.n_patches * 13.0 + non_zero.size()) / (3 * mesh->n_vertices());
	gproshan_log_var(ratio);

	return max_error;
}

che * inpainting::point_cloud_reconstruction(real_t per, real_t fr)
{
	arma::mat S;
	arma::mat alpha;

	S.load(tmp_file_path(key_name + ".smp"));
	alpha.load(tmp_file_path(key_name + ".alpha"));
	
	m_params.n_patches = S.n_rows;
	patches.resize(m_params.n_patches);

	gproshan_debug_var(m_params.n_patches);

	real_t max_radio = -1;

	#pragma omp parallel for reduction(max: max_radio)
	for(index_t i = 0; i < m_params.n_patches; i++)
		max_radio = max(max_radio, S(i, 3)); 

	A.eye(phi_basis->dim(), m_params.n_atoms);
	
	size_t total_points = 0;
	
	#pragma omp parallel for
	for(index_t i = 0; i < m_params.n_patches; i++)
	{
		arma::mat T(3,3);
		T(0,0) = S(i,4);
		T(1,0) = S(i,5);
		T(2,0) = S(i,6);
		T(0,1) = S(i,7);
		T(1,1) = S(i,8);
		T(2,1) = S(i,9);
		T(0,2) = S(i,10);
		T(1,2) = S(i,11);
		T(2,2) = S(i,12);
		
		patch & p = patches[i];

		p.init_random({S(i, 0), S(i, 1), S(i, 2)}, T, S(i, 3), max_radio, per, fr);
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
		
		a_vec x = p.phi * A * alpha.col(i);
		p.xyz.row(2) = x.t();
		
		p.iscale_xyz(phi_basis->radio());
		p.itransform();

		#pragma omp atomic
		total_points += p.xyz.n_cols;
		
		for(index_t j = 0; j < patches[i].xyz.n_cols; j++)
			if(patches[i].xyz(2, j) > 2)
			{
				gproshan_debug_var(i);
				break;
			}
	}

	gproshan_debug_var(total_points);

	vector<vertex> point_cloud;
	point_cloud.reserve(total_points);

	for(index_t i = 0; i < m_params.n_patches; i++)
	for(index_t j = 0; j < patches[i].xyz.n_cols; j++)
		point_cloud.push_back({	patches[i].xyz(0, j), 
								patches[i].xyz(1, j),
								patches[i].xyz(2, j)
								});

	che * new_mesh = new che(point_cloud.data(), point_cloud.size(), nullptr, 0);
	new_mesh->update_normals();
	
	a_vec n;
	for(index_t v = 0, i = 0; i < m_params.n_patches; i++)
	{
		patch & p = patches[i];

		phi_basis->d_discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t(), 0);
		a_vec dx = p.phi * A * alpha.col(i);
		
		phi_basis->d_discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t(), 1);
		a_vec dy = p.phi * A * alpha.col(i);

		for(index_t j = 0; j < patches[i].xyz.n_cols; j++, v++)
		{
			n = {-dx(j), -dy(j), 1};
			n = normalise(p.T * n);
			new_mesh->normal(v) = {n(0), n(1), n(2)};
		}	
	}

	che_off::write_file(new_mesh, tmp_file_path(key_name + '_' + to_string(per) + '_' + to_string(fr) + "_pc"), che_off::NOFF);
	
	gproshan_debug(saved new point cloud);
	
	return new_mesh;
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

	return 0;
}


} // namespace gproshan::mdict

