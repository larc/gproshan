#ifndef D_MESH_H
#define D_MESH_H

#include "include.h"
#include "che.h"
#include "patch.h"
#include "geodesics.h"

#include "include_arma.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


typedef void * params_t[];
typedef void (* phi_function_t) (a_mat &, a_mat &, params_t);

typedef map< index_t, index_t > patches_map_t;

struct patch_t;

void jet_fit_directions(patch_t & rp);
void PCA(patch_t & rp);
void principal_curvatures(patch_t & rp, che * mesh);

struct patch_t
{
	static bool del_index;
	static size_t min_nvp;

	size_t n;
	index_t * indexes;
	a_mat xyz;
	a_vec avg;
	a_mat E;
	a_mat phi;

	patch_t()
	{
		indexes = nullptr;
	}

	~patch_t()
	{
		if(del_index)
		if(indexes) delete [] indexes;
	}

	index_t operator[](index_t i)
	{
		return indexes[i];
	}

	// xyz = E.t * (xyz - avg)
	void transform()
	{
		xyz.each_col() -= avg;
		xyz = E.t() * xyz;
	}

	void itransform()
	{
		xyz = E * xyz;
		xyz.each_col() += avg;
	}

	bool valid_xyz()
	{
		return xyz.n_cols > min_nvp;
	}

	void reset_xyz(che * mesh, std::vector<patches_map_t> & patches_map, const index_t & p, const index_t & threshold = NIL)
	{
		size_t m = n;
		if(threshold != NIL)
		{
			m = 0;
			for(index_t i = 0; i < n; i++)
				if(indexes[i] < threshold) m++;
		}

		xyz.set_size(3, m);
		for(index_t j = 0, i = 0; i < n; i++)
		{
			if(indexes[i] < threshold)
			{
				const vertex & v = mesh->gt(indexes[i]);
				xyz(0, j) = v.x;
				xyz(1, j) = v.y;
				xyz(2, j) = v.z;

				patches_map[indexes[i]][p] = j++;
			}
		}
	}
};

a_vec gaussian(a_mat & xy, real_t sigma, real_t cx, real_t cy);

a_vec cossine(a_mat & xy, distance_t radio, size_t K);

void phi_gaussian(a_mat & phi, a_mat & xy, void ** params);

void get_centers_gaussian(a_vec & cx, a_vec & cy, real_t radio, size_t K);

void save_patches_coordinates(std::vector<patch_t> & patches, std::vector<std::pair<index_t,index_t> > * lpatches, size_t NV);

void save_patches(std::vector<patch_t> patches, size_t M);

void partial_mesh_reconstruction(size_t old_n_vertices, che * mesh, size_t M, std::vector<patch_t> & patches, std::vector<patches_map_t> & patches_map, a_mat & A, a_mat & alpha);

distance_t mesh_reconstruction(che * mesh, size_t M, std::vector<patch> & patches, std::vector<vpatches_t> & patches_map, a_mat & A, a_mat & alpha, distance_t * dist,  const fmask_t & mask = nullptr, const index_t & v_i = 0);

a_vec non_local_means_vertex(a_mat & alpha, const index_t & v, std::vector<patch> & patches, std::vector<vpatches_t> & patches_map, const distance_t & h);

[[deprecated]]
void mesh_reconstruction(che * mesh, size_t M, std::vector<patch_t> & patches, std::vector<patches_map_t> & patches_map, a_mat & A, a_mat & alpha, const index_t & v_i = 0);

[[deprecated]]
a_vec non_local_means_vertex(a_mat & alpha, const index_t & v, std::vector<patch_t> & patches, std::vector<patches_map_t> & patches_map, const distance_t & h);

a_vec simple_means_vertex( const index_t & v, std::vector<patch> & patches, std::vector<vpatches_t> & patches_map);


} // namespace gproshan::mdict

#endif // D_MESH_H

