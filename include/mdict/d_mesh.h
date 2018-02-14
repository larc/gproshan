#ifndef D_MESH_H
#define D_MESH_H

#include "include.h"
#include "che_off.h"
#include "geodesics.h"

#include <armadillo>

#ifndef CGAL_PATCH_DEFS
	#define CGAL_PATCH_DEFS
	#define CGAL_EIGEN3_ENABLED
	#define CGAL_USE_BOOST_PROGRAM_OPTIONS
	#define CGAL_USE_GMP
	#define DCGAL_USE_MPFR
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>

using namespace arma;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

typedef void * params_t[];
typedef void (* phi_function_t) (mat &, mat &, params_t);

typedef map< index_t, index_t > patches_map_t;


typedef vertex_t DFT;
typedef CGAL::Simple_cartesian<DFT> Data_Kernel;
typedef Data_Kernel::Point_3 DPoint;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;

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
	mat xyz;
	vec avg;
	mat E;
	mat phi;

	patch_t()
	{
		indexes = NULL;
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

	void reset_xyz(che * mesh, vector<patches_map_t> & patches_map, const index_t & p, const index_t & threshold = NIL)
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

vec gaussian(mat & xy, vertex_t sigma, vertex_t cx, vertex_t cy);

vec cossine(mat & xy, distance_t radio, size_t K);

void phi_gaussian(mat & phi, mat & xy, void ** params);

void get_centers_gaussian(vec & cx, vec & cy, vertex_t radio, size_t K);

void save_patches_coordinates( vector<patch_t> & patches, vector< pair<index_t,index_t> > * lpatches, size_t NV);

void save_patches(vector<patch_t> patches, size_t M);

void partial_mesh_reconstruction(size_t old_n_vertices, che * mesh, size_t M, vector<patch_t> & patches, vector<patches_map_t> & patches_map, mat & A, mat & alpha);

void mesh_reconstruction(che * mesh, size_t M, vector<patch_t> & patches, vector<patches_map_t> & patches_map, mat & A, mat & alpha, const index_t & v_i = 0);

vec non_local_means_vertex(mat & alpha, const index_t & v, vector<patch_t> & patches, vector<patches_map_t> & patches_map, const distance_t & h);

vec simple_means_vertex(mat & alpha, const index_t & v, vector<patch_t> & patches, vector<patches_map_t> & patches_map, const distance_t & h);

} // mdict

#endif // D_MESH_H

