#include "dictionary.h"

#include "sampling.h"
#include "mdict.h"
#include "che_poisson.h"
#include "che_fill_hole.h"

#include <cassert>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


size_t dictionary::L = 10;
size_t dictionary::T = 5;

dictionary::dictionary(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _d_plot):
					mesh(_mesh), phi_basis(_phi_basis), m(_m), M(_M), f(_f), d_plot(_d_plot)
{
	A.eye(phi_basis->dim, m);
}

dictionary::~dictionary()
{
	patch_t::del_index = true;
}

void dictionary::learning()
{
	gproshan_debug(MDICT);

	string f_dict = tmp_file_path(mesh->name_size() + '_' + to_string(phi_basis->dim) + '_' + to_string(m) + ".dict");
	gproshan_debug_var(f_dict);

	if(!A.load(f_dict))
	{
		A.eye(phi_basis->dim, m);
		// A.random(phi_basis->dim, m);

		KSVDT(A, patches, M, L);
		A.save(f_dict);
	}

	assert(A.n_rows == phi_basis->dim);
	assert(A.n_cols == m);

	if(d_plot)
	{
		phi_basis->plot_basis();
		phi_basis->plot_atoms(A);
	}
}

void dictionary::sparse_coding()
{
	gproshan_debug(MDICT);

	alpha.zeros(m, M);
	OMP_all_patches_ksvt(alpha, A, patches, M, L);
}

void dictionary::init_sampling()
{
	gproshan_debug(MDICT);

	n_vertices = mesh->n_vertices();

	// load sampling
	if(M == 0)
	{
		M = mesh->n_vertices();
		phi_basis->radio = T * mesh->mean_edge();
	}
	else
	{
		sampling.reserve(M);
		if(!load_sampling(sampling, phi_basis->radio, mesh, M))
			cerr << "Failed to load sampling" << endl;
	}

	s_radio = phi_basis->radio;
	//phi_basis->radio *= f;
}

void dictionary::init_patches(const bool & reset, const fmask_t & mask)
{
	gproshan_debug(MDICT);

	if(reset)
	{
		patches.resize(M);
		patches_map.resize(n_vertices);

		#pragma omp parallel
		{
			index_t * toplevel = new index_t[n_vertices];

			#pragma omp for 
			for(index_t s = 0; s < M; s++)
			{
				index_t v = sample(s);
				patches[s].init(mesh, v, dictionary::T, phi_basis->radio);
			}

			delete [] toplevel;
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
			gproshan_debug_var(patch_avg_size);
			gproshan_debug_var(patch_min_size);
			gproshan_debug_var(patch_max_size);
		#endif
	}

	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz(mesh, patches_map, s, mask);

	#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];

		p.transform();
		p.phi.set_size(p.xyz.n_cols, phi_basis->dim);
		phi_basis->discrete(p.phi, p.xyz);
	}
}

void dictionary::mesh_reconstruction()
{
	gproshan_debug(MDICT);

	assert(n_vertices == mesh->n_vertices());
	mdict::mesh_reconstruction(mesh, M, patches, patches_map, A, alpha);
}

index_t dictionary::sample(const index_t & s)
{
	assert(s < M);
	if(sampling.size()) return sampling[s];
	return s;
}


} // namespace gproshan::mdict

