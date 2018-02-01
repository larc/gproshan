#include "dictionary.h"

#include "sampling.h"
#include "d_dict_learning.h"

#include "che_poisson.h"
#include "che_fill_hole.h"
#include <cassert>
#include <set>

// mesh dictionary learning and sparse coding namespace
namespace mdict {

const size_t dictionary::min_nvp = 36;
const size_t dictionary::L = 10;

dictionary::dictionary(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & f, const bool & _d_plot):
					mesh(_mesh), phi_basis(_phi_basis), m(_m), M(_M), d_plot(_d_plot)
{
	n_vertices = mesh->n_vertices();

	// load sampling
	if(M == 0)
	{
		M = mesh->n_vertices();
		phi_basis->radio = 3 * mesh->mean_edge();
	}
	else
	{
		sampling.reserve(M);
		assert(load_sampling(sampling, phi_basis->radio, mesh, M));
	}

	s_radio = phi_basis->radio;
	phi_basis->radio *= f;

	patches.resize(M);
	patches_map.resize(n_vertices);

	patch_t::del_index = false;
	init_patches();

	A.eye(phi_basis->dim, m);
	alpha.zeros(m, M);

	if(d_plot) phi_basis->plot_basis();
}

dictionary::~dictionary()
{
	patch_t::del_index = true;
}

void dictionary::learning()
{
	string f_dict = "tmp/" + mesh->name() + '_' + to_string(phi_basis->dim) + '_' + to_string(m) + ".dict";
	debug(f_dict)

	if(!A.load(f_dict))
	{
		A.eye(phi_basis->dim, m);
		// A.random(phi_basis->dim, m);

		d_message(dictionary learning...)
		TIC(d_time) KSVDT(A, patches, M, L); TOC(d_time)
		debug(d_time)

		A.save(f_dict);
	}

	assert(A.n_rows == phi_basis->dim);
	assert(A.n_cols == m);

	if(d_plot) phi_basis->plot_atoms(A);
}
/*
void dictionary::inpaiting()
{
	d_message(iterative inpainting)
	// Computing hole borders
	size_t n_borders = mesh->n_borders();
	if(!n_borders) return;

	vector<index_t> * border_vertices = fill_all_holes(mesh);

	//with poisson
	poisson(mesh, n_vertices, 2);
	debug(mesh->n_vertices())

	//Contains index of patches which are borders
	set<index_t> border_patches;	
		
	for(index_t nb = 0; nb < n_borders; nb++)
		for(index_t b: border_vertices[nb])
			for(auto p_aux: patches_map[b])
				border_patches.insert(p_aux.first);

	delete [] border_vertices;
		
	debug(border_patches.size())
	
	// First iteration
	index_t v;
	size_t size;
	vec tmp_alpha(m);
	tmp_alpha.zeros();

	patches_map.resize(mesh->n_vertices());

	for(index_t bp: border_patches)
	{
		//updating patch
		patch_t & p = patches[bp];
			
		v = p[0];
		geodesics fm(mesh, { v }, NIL, phi_basis->radio);

		delete [] p.indexes;
		p.indexes = new index_t[fm.get_n_radio()];
	
		fm.get_sort_indexes(p.indexes, fm.get_n_radio());
		size = fm.get_n_radio();

		p.n = size;
		p.reset_xyz(mesh, patches_map, bp);

	//	principal_curvatures(p, mesh);
		if(p.n > min_nvp)
		{
			jet_fit_directions(p);
		}
		else
			principal_curvatures(p, mesh);

		p.transform();

		p.phi.set_size(p.n, phi_basis->dim);
		phi_basis->discrete(p.phi, p.xyz);
	}
		
	// Second iteration
	size_t M_;
	if(sampling.size()) M_ = mesh->n_vertices();
	else
	{	
		float time_g;
		farthest_point_sampling_gpu(sampling, time_g, mesh, n_vertices, s_radio);
		debug(time_g)

		M_ = sampling.size();
	}

	patches.resize(M_);
	alpha.set_size(m, M_);
	patches_map.resize(mesh->n_vertices());

	debug(M_)	
	debug(M)

	index_t * levels = new index_t[M_];
	memset(levels, 0, sizeof(index_t) * M);
	memset(levels + M, 255, sizeof(index_t) * (M_ - M));


//
	// Including all patches to the mapping...
	debug_me(it inpaiting)

	index_t count = M;
	index_t level = 0;
	auto is_level = [&](const index_t & v, const index_t & level) -> bool
	{
		for(auto & m: patches_map[v])
			if(levels[m.first] < level)
				return true;

		return false;
	};

	while(count < M_)
	{
		level++;
		debug(level)
		for(index_t s = M; s < M_; s++)
		{
			v = sample(s);
			patch_t & p = patches[s];

			if(levels[s] == NIL && is_level(v, level))
			{	
				levels[s] = level;
				count++;
				
				geodesics fm(mesh, { v }, NIL, phi_basis->radio );
				p.n = fm.get_n_radio();
				p.indexes = new index_t[p.n];
				fm.get_sort_indexes(p.indexes, p.n);
					
				p.reset_xyz(mesh, patches_map, s);

		//		jet_fit_directions(p);
				if(p.n > min_nvp)
				{
					jet_fit_directions(p);
				}
				else
					principal_curvatures(p, mesh);
			//	principal_curvatures(p, mesh);
					
				p.transform();

				p.phi.set_size(p.n, phi_basis->dim);
				phi_basis->discrete(p.phi, p.xyz);
			
			//	OMP_patch(alpha, A, s, p, L);	
				
				map<index_t, index_t> patches_alphas;

				for(index_t i = 0; i < p.n; i++)
				for(auto & pi: patches_map[p[i]])
			//		if(pi.first != s) patches_alphas[pi.first]++;
					if(levels[pi.first] < level) patches_alphas[pi.first]++;
				
				alpha.col(s).zeros();
				distance_t sum = 0;
				for(auto & pa: patches_alphas)
				{
					alpha.col(s) += pa.second * alpha.col(pa.first);
					sum += pa.second;
				}

				assert(sum > 0);
				alpha.col(s) /= sum;
			}
		}
		debug(count)
	}
		
	debug(count)

	debug(patches.size() - M)
		
	M = M_;

	n_vertices = mesh->n_vertices();

	delete [] levels;
	
}*/
/*
void dictionary::denoising()
{
	d_message(sparse coding...)
	TIC(d_time)
	OMP_all_patches_ksvt(alpha, A, patches, M, L);
	TOC(d_time)

	d_message(mesh reconstruction...)
	assert(n_vertices == mesh->n_vertices());

	TIC(d_time)
	mesh_reconstruction(mesh, M, patches, patches_map, A, alpha);
	TOC(d_time)
	debug(d_time)
}
*/
void dictionary::init_patches()
{
	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		index_t v = sample(s);
		patch_t & p = patches[s];

		geodesics fm(mesh, {v}, NIL, phi_basis->radio);

		p.n = fm.get_n_radio();
		p.indexes = new index_t[p.n];
		fm.get_sort_indexes(p.indexes, p.n);
	}
	
	for(index_t s = 0; s < M; s++)
	{
		patch_t & p = patches[s];
		p.reset_xyz(mesh, patches_map, s);
	}

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch_t & p = patches[s];
		
		assert(p.n > min_nvp); // old code change to principal_curvatures
		jet_fit_directions(p);

		p.transform();
		p.phi.set_size(p.n, phi_basis->dim);
		phi_basis->discrete(p.phi, p.xyz);
	}

	#ifndef NDEBUG
		size_t patch_mean_size = 0;
	
		#pragma omp parallel for reduction(+: patch_mean_size)
		for(index_t s = 0; s < M; s++)
			patch_mean_size += patches[s].n;
		
		patch_mean_size /= M;
		debug(patch_mean_size)
	#endif
}

index_t dictionary::sample(const index_t & s)
{
	assert(s < M);
	if(sampling.size()) return sampling[s];
	return s;
}

} // mdict

