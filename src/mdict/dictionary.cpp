#include "dictionary.h"

#include "sampling.h"
#include "mdict.h"
#include "che_poisson.h"
#include "che_fill_hole.h"

#include "viewer/viewer.h"

#include <cassert>
#include <CImg.h>


using namespace cimg_library;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


size_t dictionary::L = 12;
size_t dictionary::K = 10;
size_t dictionary::T = 5;

dictionary::dictionary(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, const bool & _d_plot):
					mesh(_mesh), phi_basis(_phi_basis), m(_m), M(_M), f(_f), learn(_learn), d_plot(_d_plot)
{
	A.eye(phi_basis->dim, m);
	dist = new distance_t[mesh->n_vertices()]; 
}

dictionary::~dictionary()
{
	patch_t::del_index = true;
}

void dictionary::learning()
{
	gproshan_log(MDICT);

	string f_dict = tmp_file_path(mesh->name_size() + '_' + to_string(phi_basis->dim) + '_' + to_string(m) + '_' + to_string(f) + '_' + to_string(L) + ".dict");

	if(learn)
	{
		gproshan_log_var(f_dict);

		if(!A.load(f_dict))
		{
			A.eye(phi_basis->dim, m);
			A = normalise(A);
			gproshan_debug_var(phi_basis->radio);
			gproshan_debug_var(m);
			phi_basis->plot_atoms(A);
			KSVD(A, patches, L, K);
			A.save(f_dict);
		}
	}
	else A.eye(phi_basis->dim, m);
	gproshan_debug_var(phi_basis->radio);
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
	gproshan_log(MDICT);
	
	vector<locval_t> locval;
	alpha = OMP_all(patches, phi_basis, A, L);
}

void dictionary::init_sampling()
{
	gproshan_log(MDICT);

	n_vertices = mesh->n_vertices();

	// load sampling
	if(M == 0)
	{
		M = mesh->n_vertices();
		phi_basis->radio = mesh->mean_edge();
	}
	else
	{
		sampling.reserve(M);
		if(!load_sampling(sampling, phi_basis->radio, mesh, M))
			cerr << "Failed to load sampling" << endl;
	}

	s_radio = phi_basis->radio;
	phi_basis->radio *= f;

	gproshan_debug_var(s_radio);
	gproshan_debug_var(phi_basis->radio);
}

void dictionary::init_patches(const bool & reset, const fmask_t & mask)
{
	gproshan_log(MDICT);

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
				patches[s].init(mesh, v, dictionary::T, phi_basis->radio, toplevel);
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
		p.phi = normalise(p.phi);
	}

/*	
#ifndef NDEBUG
	CImgList<real_t> imlist;
	for(index_t s = 0; s < M; s++)
		patches[s].save(phi_basis->radio, 16, imlist);
	imlist.save_ffmpeg_external("tmp/patches.mpg", 5);
#endif	

*/

	/*Saving Patches*/
/*
	ofstream os(tmp_file_path("patch-mat"));
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];
		p.save_z(os);
	}
	os.close();
	// DRAW NORMALS DEBUG
	for(index_t s = 0; s < M; s++)
	{
		viewer::vectors.push_back({patches[s].x(0), patches[s].x(1), patches[s].x(2)});
		a_vec r = patches[s].x + 0.02 * patches[s].normal();
		viewer::vectors.push_back({r(0), r(1), r(2)});
	}
	*/
}

distance_t dictionary::mesh_reconstruction(const fmask_t & mask)
{
	gproshan_log(MDICT);

	assert(n_vertices == mesh->n_vertices());
	return mdict::mesh_reconstruction(mesh, M, patches, patches_map, A, alpha, dist);
}
void dictionary::update_alphas(a_mat & alpha, size_t threshold)
{
	size_t np_new = M - threshold;
	bool patches_covered[np_new];
	memset(patches_covered, 0, sizeof(patches_covered));
	size_t count = 0;

	// Choose the border patches using the threshold
	while(count < threshold)
	{
		#pragma omp parallel for
		for(index_t s = threshold; s < M; s++)
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
				count++;
			}	

		}
	}
	
	// update alphas of choosed patches
	// update the threshold
	// repeat until threshold reachs all patches
}

index_t dictionary::sample(const index_t & s)
{
	assert(s < M);
	if(sampling.size()) return sampling[s];
	return s;
}

const distance_t & dictionary::operator[](const index_t & i) const
{
	assert(i < mesh->n_vertices());
	return dist[i];
}

void dictionary::draw_patches(index_t i)
{
	phi_basis->plot_patch(A*alpha.col(i),patches[i].xyz, i);
}


} // namespace gproshan::mdict

