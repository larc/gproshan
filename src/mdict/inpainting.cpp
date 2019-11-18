#include "inpainting.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


inpainting::inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, const bool & _plot): dictionary(_mesh, _phi_basis, _m, _M, _f, _learn, _plot)
{
}

distance_t inpainting::execute()
{
	// M: the number of seeds
	M = mesh->n_vertices()/36;
	gproshan_log_var(M);
	gproshan_debug_var(M);

	//FPS sampling with desired number of sources
	TIC(d_time) init_sampling(); TOC(d_time)
	gproshan_debug_var(d_time);

	// creating disjoint clusters with geodesics aka voronoi
#ifdef GPROSHAN_CUDA
	geodesics ptp( mesh, sampling, geodesics::PTP_GPU, nullptr, 1);
#else
	geodesics ptp( mesh, sampling, geodesics::FM, nullptr, 1);
#endif
	TOC(d_time)
	gproshan_log_var(d_time);

	//mapping the vertices to disjoint patches 1 vertex belongs only to one patch
	// creating new patches
	std::vector<index_t> vertices[M];

	//saving first vertex aka seed vertices
	#pragma omp for 
	for(index_t s = 0; s < M; s++)
	{
		vertices[s].push_back(sample(s));
	}

	bool * mask = new bool[mesh->n_vertices()];
	size_t * percentages_size = new size_t[M];
	
	// vertices contains the mapping.
	//#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		mask[i] = 0;
		ptp.clusters[i]--;
		if(sample(ptp.clusters[i]) != i)
			vertices[ ptp.clusters[i] ].push_back(i) ;
	}
	//Randomly remove a percentage of points for each patch
	// create initial desired percentage sizes
	size_t percent = 50;

	#pragma omp for 
	for(index_t s = 0; s < M; s++)
	{
		percentages_size[s] = ceil(vertices[s].size() * percent/ 100) ;
	}

	std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution(0,M);

	int k = 0;
	size_t rn=0;
	while( k < M )
	{
		rn = distribution(generator);
		if(!mask[rn] && percentages_size[ptp.clusters[rn] ] > 0)
		{
			mask[rn] = 1;
			percentages_size[ ptp.clusters[rn] ]--;			
		}
		if(percentages_size[ptp.clusters[rn] ] == 0)	
			k++;
	}
	//Initializing patches
	gproshan_log(our mask is ready);

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

	
	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz_disjoint(mesh, dist, patches_map, s, 0 ,[&mask](const index_t & i) -> bool { return mask[i]; } );

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->get_dim());
		phi_basis->discrete(p.phi, p.xyz);

	} 

	gproshan_log(our patches are ready); 

	// sparse coding and reconstruction with all patches
	TIC(d_time) sparse_coding(); TOC(d_time)
	gproshan_debug_var(d_time);

	//patches_map.resize(n_vertices);
	for(index_t  i = 0; i < n_vertices; i++)
	{
		patches_map[i].clear();
	}

	for(index_t s = 0; s < M; s++)
//	patches[s].reset_xyz_disjoint(mesh, dist, patches_map, s, 1 ,[&mask](const index_t & i) -> bool { return !mask[i]; } );
		patches[s].reset_xyz(mesh, patches_map, s, 0);

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->get_dim());
		phi_basis->discrete(p.phi, p.xyz);
	}

	TIC(d_time) mesh_reconstruction(); TOC(d_time)
	gproshan_debug_var(d_time);

}

distance_t inpainting::execute_tmp()
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

