#include <gproshan/mdict/msparse_coding.h>

#include <gproshan/mdict/mdict.h>
#include <gproshan/mesh/che_off.h>
#include <gproshan/mesh/che_poisson.h>
#include <gproshan/mesh/che_fill_hole.h>
#include <gproshan/geodesics/sampling.h>

#include <cassert>
#include <fstream>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


size_t msparse_coding::L = 12;
size_t msparse_coding::K = 10;
size_t msparse_coding::T = 5;

msparse_coding::msparse_coding(che *const & _mesh, basis *const & _phi_basis, const params & p): mesh(_mesh), phi_basis(_phi_basis), m_params(p)
{
	A.eye(phi_basis->dim(), m_params.n_atoms);
	dist = new real_t[mesh->n_vertices];

	m_params.n_patches = mesh->n_vertices / m_params.avg_p;

	key_name = mesh->name_size() + '_' + std::to_string(m_params.delta) + '_' + std::to_string(m_params.sum_thres) + '_' + std::to_string(m_params.area_thres);

	mask = new bool[mesh->n_vertices];
	memset(mask, 0, sizeof(bool) * mesh->n_vertices);
}

msparse_coding::~msparse_coding()
{
	delete [] mask;
}

msparse_coding::operator const std::string & () const
{
	return key_name;
}

void msparse_coding::load_mask()
{
	//std::string f_mask = tmp_file_path(mesh->name_size() + '_' + std::to_string(avg_p) + '_' + std::to_string(percent) + '_' + std::to_string(radio) + ".msk");
	std::string f_mask = tmp_file_path(mesh->name_size() + '_' + std::to_string(m_params.percent) + ".msk");
	arma::uvec V;
	gproshan_log(loading radial mask);


	if(V.load(f_mask))
	{
		#pragma omp for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
		{
			mask[i] = V(i);
		}
	}
	else
	{
		V.zeros(mesh->n_vertices);
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0, mesh->n_vertices-1);
		size_t percentage = mesh->n_vertices - ceil(mesh->n_vertices * (m_params.percent / 100.0)) ;

		size_t k = 0;
		size_t rn = 0;
		while (k < percentage)
		{
			rn = distribution(generator);
			if(!mask[rn])
			{
				mask[rn] = 1;
				V(rn) = 1;
				++k;
			}
		}
		V.save(f_mask);
	}
}

void msparse_coding::load_mask(const std::vector<index_t> * vertices, const index_t * clusters)
{
	std::string f_mask = tmp_file_path(mesh->name_size() + '_' + std::to_string(m_params.avg_p) + '_' + std::to_string(m_params.percent) + ".msk");
	arma::uvec V;

	if(V.load(f_mask))
	{
		#pragma omp for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
		{
			mask[i] = V(i);
		}

	}
	else
	{

		V.zeros(mesh->n_vertices);
		size_t * percentages_size = new size_t[m_params.n_patches];
		bool cover_cluster[m_params.n_patches];


		// create initial desired percentage sizes
		#pragma omp for
		for(index_t s = 0; s < m_params.n_patches; ++s)
		{
			percentages_size[s] = ceil(size(vertices[s]) * (m_params.percent / 100.0)) ;
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
				++k;
			}

		}


		V.save(f_mask);

	}

}

void msparse_coding::load_sampling()
{
	size_t featsize;
	std::vector<index_t> all_sorted_features;
	std::vector<index_t> seeds;
	load_features(all_sorted_features, featsize);

	gproshan_debug_var(size(all_sorted_features));

	size_t count = 0;
	real_t area_mesh = mesh->area_surface();

	a_vec S;
	if(S.load(tmp_file_path(key_name + ".rsampl")))
	{
		gproshan_debug(loading sampling);

		size_t n_seeds = S.n_rows;
		real_t euc_radio, geo_radio;
		for(index_t i = 0; i < n_seeds; ++i)
		{
			patch p;
			p.init_radial_disjoint(euc_radio, geo_radio, mesh, S(i), m_params.delta, m_params.sum_thres, m_params.area_thres, area_mesh);
			patches.push_back(std::move(p));
		}

		m_params.n_patches = n_seeds;

		return;
	}

	//ELSE IF !S.load(f_sample)

	gproshan_debug(computing sampling);

	// compute features will be seeds
	patches_map.resize(mesh->n_vertices);

	//Coverage of the points
	bool covered[mesh->n_vertices];

	#pragma omp for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
		covered[i] = 0;

	real_t euc_radio;
	real_t geo_radio;
	std::vector<real_t> radios;
	std::vector<real_t> geo_radios;
	size_t count_cov = 0;
	size_t count_cov_patch = 0;

	bool * invalid_seed = new bool[mesh->n_vertices];
	memset(invalid_seed, 0, mesh->n_vertices * sizeof(bool));

	for(const index_t & vsf: all_sorted_features)
	{
		if(invalid_seed[vsf]) continue;

		patch p;
		p.init_radial_disjoint(euc_radio, geo_radio, mesh, vsf, m_params.delta, m_params.sum_thres, m_params.area_thres, area_mesh);

		count_cov_patch = 0;
		if(size(p.vertices) >= 7 )
		{
			for(const index_t & v: p.vertices)
				if(!covered[v]) ++count_cov_patch;

			count_cov += count_cov_patch;
			if(count_cov_patch > 0)
			{
				patches.push_back(std::move(p));
				seeds.push_back(vsf);
				radios.push_back(euc_radio);
				geo_radios.push_back(geo_radio);
				count += p.vertices.size();

				for(const index_t & v: p.vertices)
				{
					covered[v] = 1;

					if(!invalid_seed[v])
					{
						const vertex & va = mesh->point(vsf);
						const vertex & vb = mesh->point(v);

						invalid_seed[v] = norm(va - vb) < 0.8 * euc_radio;
					}
				}
			}
		}
	}

	delete [] invalid_seed;

	m_params.n_patches = seeds.size();

	gproshan_log(saving sampling);

	S.resize(size(seeds));

	#pragma omp parallel for
	for(index_t i = 0; i < seeds.size(); ++i)
		S(i) = seeds[i];

	S.save(tmp_file_path(key_name + ".rsampl"));

	gproshan_debug_var(m_params.sum_thres);
	gproshan_log_var(count);
	gproshan_log_var(count_cov);
	gproshan_debug_var(size(seeds));
	gproshan_debug_var(m_params.n_patches);


#ifndef NDEBUG
	std::vector<index_t> outliers;


	for(index_t i = 0; i < mesh->n_vertices; ++i)
		if(!covered[i])
			outliers.push_back(i);

	a_vec outlv(size(outliers));
	for(index_t i = 0; i < outliers.size(); ++i)
		outlv(i) = outliers[i];

	outlv.save(tmp_file_path(key_name + ".outlv"));
#endif // NDEBUG
}

void msparse_coding::init_radial_feature_patches()
{
	load_sampling();
	load_mask();

	#ifndef NDEBUG
		size_t patch_avg_size = 0;
		size_t patch_min_size = NIL;
		size_t patch_max_size = 0;

		#pragma omp parallel for reduction(+: patch_avg_size)
		for(index_t s = 0; s < m_params.n_patches; ++s)
			patch_avg_size += patches[s].vertices.size();
		#pragma omp parallel for reduction(min: patch_min_size)
		for(index_t s = 0; s < m_params.n_patches; ++s)
			patch_min_size = std::min(size(patches[s].vertices), patch_min_size);
		#pragma omp parallel for reduction(max: patch_max_size)
		for(index_t s = 0; s < m_params.n_patches; ++s)
			patch_max_size = std::max(size(patches[s].vertices), patch_max_size);

		patch_avg_size /= m_params.n_patches;
		gproshan_debug_var(patch_avg_size);
		gproshan_debug_var(patch_min_size);
		gproshan_debug_var(patch_max_size);
	#endif

	//Initializing patches
	gproshan_log(initializing patches);
	gproshan_debug_var(m_params.n_patches);

	n_vertices = mesh->n_vertices;
	patches.resize(m_params.n_patches);
	patches_map.resize(n_vertices);

	bool * pmask = mask;

	for(index_t s = 0; s < m_params.n_patches; ++s)
		patches[s].reset_xyz_disjoint(mesh, dist, m_params.n_patches, patches_map, s, [&pmask](const index_t & i) -> bool { return pmask[i]; });

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; ++s)
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
		a_mat AS;
		AS.resize(m_params.n_patches, 13);
		for(index_t i = 0; i < m_params.n_patches; ++i)
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

void msparse_coding::init_voronoi_patches()
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
		for(index_t s = 0; s < m_params.n_patches; ++s)
		{
			vertices[s].push_back(sample(s));
		}

		for(index_t i = 0; i < mesh->n_vertices; ++i)
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
		index_t * toplevel = new index_t[mesh->n_vertices];

		#pragma omp for
		for(index_t s = 0; s < m_params.n_patches; ++s)
		{
			patches[s].init_disjoint(mesh, sample(s), msparse_coding::T, vertices[s], toplevel);

		}
		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_min_size = std::min(size(patches[s].vertices), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_max_size = std::max(size(patches[s].vertices), patch_max_size);

			patch_avg_size /= m_params.n_patches;
			//gproshan_debug_var(patch_avg_size);
			//gproshan_debug_var(patch_min_size);
			//gproshan_debug_var(patch_max_size);
		#endif
	}

	bool * pmask = mask;
	for(index_t s = 0; s < m_params.n_patches; ++s)
		patches[s].reset_xyz_disjoint(mesh, dist, m_params.n_patches, patches_map, s, [&pmask](const index_t & i) -> bool { return pmask[i]; });

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; ++s)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
		p.compute_avg_distance(mesh, patches_map, s);

	}

	gproshan_log(our patches are ready);
}

real_t msparse_coding::execute()
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
	for(index_t i = 0; i < n_vertices; ++i)
		patches_map[i].clear();

	for(index_t s = 0; s < m_params.n_patches; ++s)
		patches[s].reset_xyz(mesh, patches_map, s, 0);

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; ++s)
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
	gproshan_debug_var(size(non_zero));
	real_t ratio = (m_params.n_patches * 13.0 + non_zero.size()) / (3 * mesh->n_vertices);
	gproshan_log_var(ratio);

	return max_error;
}

che * msparse_coding::point_cloud_reconstruction(real_t per, real_t fr)
{
	a_mat S;
	a_mat alpha;

	S.load(tmp_file_path(key_name + ".smp"));
	alpha.load(tmp_file_path(key_name + ".alpha"));

	m_params.n_patches = S.n_rows;
	patches.resize(m_params.n_patches);

	gproshan_debug_var(m_params.n_patches);

	real_t max_radio = -1;

	#pragma omp parallel for reduction(max: max_radio)
	for(index_t i = 0; i < m_params.n_patches; ++i)
		max_radio = std::max(max_radio, S(i, 3));

	A.eye(phi_basis->dim(), m_params.n_atoms);

	size_t total_points = 0;

	#pragma omp parallel for
	for(index_t i = 0; i < m_params.n_patches; ++i)
	{
		a_mat T(3,3);
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

		for(index_t j = 0; j < patches[i].xyz.n_cols; ++j)
			if(patches[i].xyz(2, j) > 2)
			{
				gproshan_debug_var(i);
				break;
			}
	}

	gproshan_debug_var(total_points);

	std::vector<vertex> point_cloud;
	point_cloud.reserve(total_points);

	for(index_t i = 0; i < m_params.n_patches; ++i)
	for(index_t j = 0; j < patches[i].xyz.n_cols; ++j)
		point_cloud.push_back({	patches[i].xyz(0, j),
								patches[i].xyz(1, j),
								patches[i].xyz(2, j)
								});

	che * new_mesh = new che(point_cloud.data(), point_cloud.size(), nullptr, 0);
	new_mesh->update_normals();

	a_vec n;
	for(index_t v = 0, i = 0; i < m_params.n_patches; ++i)
	{
		patch & p = patches[i];

		phi_basis->d_discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t(), 0);
		a_vec dx = p.phi * A * alpha.col(i);

		phi_basis->d_discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t(), 1);
		a_vec dy = p.phi * A * alpha.col(i);

		for(index_t j = 0; j < patches[i].xyz.n_cols; ++j, ++v)
		{
			n = {-dx(j), -dy(j), 1};
			n = normalise(p.T * n);
			new_mesh->normal(v) = {n(0), n(1), n(2)};
		}
	}

	che_off::write_file(new_mesh, tmp_file_path(key_name + '_' + std::to_string(per) + '_' + std::to_string(fr) + "_pc"), che_off::NOFF);

	gproshan_debug(saved new point cloud);

	return new_mesh;
}


real_t msparse_coding::execute_tmp()
{
	// fill holes
	size_t threshold = mesh->n_vertices;
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

void msparse_coding::learning()
{
	gproshan_log(MDICT);

	std::string f_dict = tmp_file_path(mesh->name_size() + '_' + std::to_string(phi_basis->dim()) + '_' + std::to_string(m_params.n_atoms) + '_' + std::to_string(m_params.f) + '_' + std::to_string(L) + ".dict");

	if(m_params.learn)
	{
		gproshan_log_var(f_dict);

		if(!A.load(f_dict))
		{
			//A.eye(phi_basis->dim(), m);
			//initialize with some alpha
			//random
			//arma::uvec r_ind = arma::randi<arma::uvec>(m, arma::distr_param(0, M));
			//A = alpha.cols(r_ind);
			a_mat R, E, U, V;
			a_vec s;
			svd(U, s, V, alpha);
			gproshan_debug(svd done!);
			A = U.cols(0, m_params.n_atoms);
			gproshan_debug(svd done!);
			A = normalise(A);
			gproshan_debug_var(phi_basis->radio());
			gproshan_debug_var(m_params.n_atoms);
			//

			phi_basis->plot_atoms(A);
			KSVD(A, patches, L, K);
			phi_basis->plot_atoms(A);
			A.save(f_dict);
		}
	}
	else A.eye(phi_basis->dim(), m_params.n_atoms);
	gproshan_debug_var(phi_basis->radio());
	assert(A.n_rows == phi_basis->dim());
	assert(A.n_cols == m_params.n_atoms);
	if(m_params.plot)
	{
		phi_basis->plot_basis();
		phi_basis->plot_atoms(A);
	}
}

void msparse_coding::sparse_coding()
{
	gproshan_log(MDICT);

	std::vector<locval_t> locval;
	alpha = OMP_all(patches, phi_basis, A, L);
}

void msparse_coding::init_sampling()
{
	gproshan_log(MDICT);

	n_vertices = mesh->n_vertices;

	// load sampling
	if(m_params.n_patches == 0)
	{
		m_params.n_patches = mesh->n_vertices;
		phi_basis->radio() = mesh->mean_edge();
	}
	else
	{
		sampling.reserve(m_params.n_patches);
		if(!gproshan::load_sampling(sampling, phi_basis->radio(), mesh, m_params.n_patches))
			std::cerr << "Failed to load sampling" << std::endl;
	}

	s_radio = phi_basis->radio();
	phi_basis->radio() *= m_params.f;

	gproshan_debug_var(s_radio);
	gproshan_debug_var(phi_basis->radio());
}

void msparse_coding::load_features(std::vector<index_t> & v_feat, size_t & featsize)
{
	std::string f_feat = tmp_file_path(mesh->name() + ".int");
	std::ifstream inp;
	inp.open(f_feat.c_str(), std::ifstream::in);

	size_t tam;
	index_t tmp;

	gproshan_debug_var(f_feat);
	if(inp.fail())
	{
		inp.clear(std::ios::failbit);
		// call the function using system
		//g++ -O3 *.cpp -lgsl -lCGAL -o harris3d
		//cmake -DCMAKE_BUILD_TYPE=Debug ..

		std::string command = "../../harris3d/harris3d " + mesh->filename + " " + f_feat + " " + tmp_file_path("example.prop");
		gproshan_debug_var(command);
		system(command.c_str());
		gproshan_debug(created);
		inp.close();
		inp.open(f_feat.c_str(), std::ifstream::in);
	}

	gproshan_debug(exists);
	inp >> featsize;
	for(index_t i = 0; i < featsize; ++i)
	{
		inp>>tmp;
		v_feat.push_back(tmp);
	}

	inp >> tam;
	for(index_t i = 0; i < tam; ++i)
	{
		inp >> tmp;
		v_feat.push_back(tmp);
	}

	inp.close();
}

void msparse_coding::init_patches(const bool & reset, const fmask_t & mask)
{
	gproshan_log(MDICT);

	if(reset)
	{
		patches.resize(m_params.n_patches);
		patches_map.resize(n_vertices);

		#pragma omp parallel
		{
			index_t * toplevel = new index_t[n_vertices];

			#pragma omp for
			for(index_t s = 0; s < m_params.n_patches; ++s)
			{
				index_t v = sample(s);
				patches[s].init(mesh, v, msparse_coding::T, phi_basis->radio(), toplevel);
			}


			delete [] toplevel;
		}

		#ifndef NDEBUG
			size_t patch_avg_size = 0;
			size_t patch_min_size = NIL;
			size_t patch_max_size = 0;

			#pragma omp parallel for reduction(+: patch_avg_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_avg_size += patches[s].vertices.size();
			#pragma omp parallel for reduction(min: patch_min_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_min_size = std::min(size(patches[s].vertices), patch_min_size);
			#pragma omp parallel for reduction(max: patch_max_size)
			for(index_t s = 0; s < m_params.n_patches; ++s)
				patch_max_size = std::max(size(patches[s].vertices), patch_max_size);

			patch_avg_size /= m_params.n_patches;
			gproshan_debug_var(patch_avg_size);
			gproshan_debug_var(patch_min_size);
			gproshan_debug_var(patch_max_size);
		#endif
	}

	for(index_t s = 0; s < m_params.n_patches; ++s)
		patches[s].reset_xyz(mesh, patches_map, s, mask);

	#pragma omp parallel for
	for(index_t s = 0; s < m_params.n_patches; ++s)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim());
		phi_basis->discrete(p.phi, p.xyz.row(0).t(), p.xyz.row(1).t());
		p.phi = normalise(p.phi);
	}

/*
#ifndef NDEBUG
	CImgList<real_t> imlist;
	for(index_t s = 0; s < m_params.n_patches; ++s)
		patches[s].save(phi_basis->radio(), 16, imlist);
	imlist.save_ffmpeg_external("tmp/patches.mpg", 5);
#endif

*/

	/*Saving Patches*/
/*
	std::ofstream os(tmp_file_path("patch-mat"));
	for(index_t s = 0; s < m_params.n_patches; ++s)
	{
		patch & p = patches[s];
		p.save_z(os);
	}
	os.close();
	// DRAW NORMALS DEBUG
	for(index_t s = 0; s < m_params.n_patches; ++s)
	{
		viewer::vectors.push_back({patches[s].x(0), patches[s].x(1), patches[s].x(2)});
		a_vec r = patches[s].x() + 0.02 * patches[s].normal();
		viewer::vectors.push_back({r(0), r(1), r(2)});
	}
	*/
}

real_t msparse_coding::mesh_reconstruction(const fmask_t & mask)
{
	gproshan_log(MDICT);

	assert(n_vertices == mesh->n_vertices);

	a_mat V(3, mesh->n_vertices, arma::fill::zeros);

	patches_error.resize(m_params.n_patches);

	#pragma omp parallel for
	for(index_t p = 0; p < m_params.n_patches; ++p)
	{
		patch & rp = patches[p];

		a_vec x = rp.phi * A * alpha.col(p);

		patches_error[p] = { accu(abs(x - rp.xyz.row(2).t())) / rp.vertices.size(), p };

		rp.xyz.row(2) = x.t();
	}

	std::sort(patches_error.begin(), patches_error.end());

	fprintf(stderr, "error %16s%16s\n", "best", "worst");
	for(index_t i = 0; i < 10; ++i)
	{
		const index_t & best = patches_error[i].second;
		const index_t & worst = patches_error[m_params.n_patches - i - 1].second;

		fprintf(stderr, "%5d:%8u>%8u%8u>%8u\n", i, best, draw_patches(best), worst, draw_patches(worst));
	}

	#pragma omp parallel for
	for(index_t p = 0; p < m_params.n_patches; ++p)
	{
		patch & rp = patches[p];
		rp.iscale_xyz(phi_basis->radio());
		rp.itransform();
	}

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		// simple means vertex
		if(size(patches_map[v]) && (!mask || mask(v)))
		{
			a_vec mv = arma::zeros<a_vec>(3);
			for(auto p: patches_map[v])
				mv += patches[p.first].xyz.col(p.second);

			V.col(v) = mv / patches_map[v].size();
		}
		else
		{
			V(0, v) = mesh->point(v).x();
			V(1, v) = mesh->point(v).y();
			V(2, v) = mesh->point(v).z();
		}
	}

	// ------------------------------------------------------------------------

	vertex * new_vertices = (vertex *) V.memptr();

	real_t error = 0;
	real_t max_error = -1;

	#pragma omp parallel for reduction(+: error) reduction(max: max_error)
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		dist[v] = length(new_vertices[v] - mesh->point(v));
		error += dist[v];
		max_error = std::max(max_error, dist[v]);
	}

	error /= mesh->n_vertices;

	gproshan_debug_var(mesh->n_vertices);
	gproshan_debug_var(error);
	gproshan_debug_var(max_error);

	mesh->update_vertices(new_vertices, mesh->n_vertices);
	che_off::write_file(mesh, "../tmp/recon_mesh");

	return max_error;
}

void msparse_coding::update_alphas(a_mat & alpha, size_t threshold)
{
	size_t np_new = m_params.n_patches - threshold;
	bool patches_covered[np_new];
	memset(patches_covered, 0, sizeof(patches_covered));
	size_t count = 0;

	// Choose the border patches using the threshold
	while(count < threshold)
	{
		#pragma omp parallel for
		for(index_t s = threshold; s < m_params.n_patches; ++s)
		{

			if(!patches_covered[s-threshold])
			{
				a_vec sum;
				sum.zeros();
				size_t c = 0;
				// Here updating alphas, we need a structure between patches and neighboor patches
				//We can simulate that structure by using patches map
				for(auto p: patches_map[s])
				{
					if(p.first < threshold || patches_covered[p.first-threshold])
					{
						sum += alpha.col(p.first);
					}
					sum /= c;

				}
				alpha.col(s) = sum;
				patches_covered[s-threshold] = 1;
				++count;
			}

		}
	}

	// update alphas of choosed patches
	// update the threshold
	// repeat until threshold reachs all patches
}

index_t msparse_coding::sample(const index_t & s)
{
	assert(s < m_params.n_patches);
	if(size(sampling)) return sampling[s];
	return s;
}

const real_t & msparse_coding::operator[](const index_t & i) const
{
	assert(i < mesh->n_vertices);
	return dist[i];
}

const index_t & msparse_coding::draw_patches(const index_t & p)
{
	phi_basis->plot_patch(A * alpha.col(p), patches[p].xyz, patches[p].vertices[0]);
	return patches[p].vertices[0];
}

void msparse_coding::save_alpha(std::string file)
{
	alpha.save(file);
}


} // namespace gproshan::mdict

