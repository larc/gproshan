 #include "d_mesh.h"
	#include "che_off.h"

#ifndef CGAL_PATCH_DEFS
	#define CGAL_PATCH_DEFS
	#define CGAL_EIGEN3_ENABLED
	#define CGAL_USE_BOOST_PROGRAM_OPTIONS
	#define CGAL_USE_GMP
	#define DCGAL_USE_MPFR
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>



// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


typedef real_t DFT;
typedef CGAL::Simple_cartesian<DFT> Data_Kernel;
typedef Data_Kernel::Point_3 DPoint;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;


size_t patch_t::min_nvp = 36;
bool patch_t::del_index = false;

a_vec gaussian(a_mat & xy, real_t sigma, real_t cx, real_t cy)
{
	a_vec x = xy.row(0).t() - cx;
	a_vec y = xy.row(1).t() - cy;

	x.for_each( [] (a_mat::elem_type & val) { val *= val; } );
	y.for_each( [] (a_mat::elem_type & val) { val *= val; } );

	return exp( - ( x + y ) / ( 2 * sigma * sigma ) );
}

a_vec cossine(a_mat & xy, real_t radio, size_t K)
{
	a_vec x = xy.row(0).t() + 0.5;
	a_vec y = xy.row(1).t() + 0.5;


	size_t k = sqrt(K);
	a_vec sum(x.n_elem, arma::fill::zeros);
	a_vec tmp;

	for(index_t nx = 0; nx < k; nx++)
	for(index_t ny = 0; ny < k; ny++)
	{
		tmp = cos( (M_PI*x*(nx-1)+nx)/radio ) % cos((M_PI*y*(ny-1)+ny)/radio );
	}	sum += tmp;

	return sum;
}

void phi_gaussian(a_mat & phi, a_mat & xy, params_t params)
{
	a_vec & cx = *( (a_vec * ) params[0] );
	a_vec & cy = *( (a_vec * ) params[1] );
	real_t sigma = *( (real_t * ) params[2] );

	size_t K = phi.n_cols;

	for(index_t k = 0 ; k < K; k++)
		phi.col(k) = gaussian(xy, sigma, cx(k), cy(k));
}

void get_centers_gaussian(a_vec & cx, a_vec & cy, real_t radio, size_t K)
{
	if(K == 1)
	{
		cx(0) = cy(0) = 0;
		return;
	}

	size_t k = sqrt(K);
	real_t d = 2 * radio / (k - 1);

	for(index_t c = 0, i = 0; i < k; i++)
	for(index_t j = 0; j < k; j++, c++)
	{
		cx(c) = -radio + d * i;
		cy(c) = -radio + d * j;
	}
}

void jet_fit_directions(patch_t & rp)
{
	vector<DPoint> in_points;
	in_points.reserve(rp.n);
	for(index_t i = 0; i < rp.xyz.n_cols; i++)
		in_points.push_back(DPoint(rp.xyz(0, i), rp.xyz(1, i), rp.xyz(2, i)));

	size_t d_fitting = 4;
	size_t d_monge = 4;

	My_Monge_form monge_form;
	My_Monge_via_jet_fitting monge_fit;
	monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

	rp.avg.set_size(3);
	rp.avg(0) = monge_form.origin()[0];
	rp.avg(1) = monge_form.origin()[1];
	rp.avg(2) = monge_form.origin()[2];
	rp.E.set_size(3,3);
	rp.E(0, 0) = monge_form.maximal_principal_direction()[0];
	rp.E(1, 0) = monge_form.maximal_principal_direction()[1];
	rp.E(2, 0) = monge_form.maximal_principal_direction()[2];
	rp.E(0, 1) = monge_form.minimal_principal_direction()[0];
	rp.E(1, 1) = monge_form.minimal_principal_direction()[1];
	rp.E(2, 1) = monge_form.minimal_principal_direction()[2];
	rp.E(0, 2) = monge_form.normal_direction()[0];
	rp.E(1, 2) = monge_form.normal_direction()[1];
	rp.E(2, 2) = monge_form.normal_direction()[2];
}

void PCA(patch_t & rp)
{
	rp.avg = mean(rp.xyz, 1);
	rp.xyz.each_col() -= rp.avg;

	a_mat C = rp.xyz * rp.xyz.t();
	a_vec eigval;
	eig_sym(eigval, rp.E, C);

	rp.E.swap_cols(0, 2);
}

void principal_curvatures( patch_t & rp, che * mesh)
{
	rp.avg = rp.xyz.col(0);

	vertex N = mesh->normal(rp[0]);

	vertex max;
	real_t K = -INFINITY;
	real_t k;

	for_star(he, mesh, rp[0])
	{
		vertex d = mesh->gt_vt(next(he)) - mesh->gt_vt(he);
		d /= *d;
		k = (N,d);

		if(K < k)
		{
			max = d;
			K = k;
		}
	}

	rp.E.set_size(3, 3);
	rp.E(0, 2) = N.x;
	rp.E(1, 2) = N.y;
	rp.E(2, 2) = N.z;

	max -= K * N;
	rp.E(0, 0) = max.x;
	rp.E(1, 0) = max.y;
	rp.E(2, 0) = max.z;

	rp.E.col(1) = cross(rp.E.col(0), rp.E.col(2));
	rp.E.col(1) /= norm(rp.E.col(1));
}

void save_patches_coordinates(vector<patch_t> & patches, vector< pair<index_t,index_t> > * lpatches, size_t NV)
{
	ofstream os(tmp_file_path("test-patches_coordinates"));

	for(index_t v = 0; v < NV; v++)
	{
		for(auto pi: lpatches[v])
		{
			patch_t & rp = patches[pi.first];
			os<<rp.xyz.col(pi.second)[0]<<" "<<rp.xyz.col(pi.second)[1]<<" "<<rp.xyz.col(pi.second)[2]<<" "<<pi.first<<" ";
		}
		os<<endl;
	}

	os.close();
}

void save_patches(vector<patch_t> & patches, size_t M)
{
	ofstream os(tmp_file_path("test-patch_wise_coordinates"));

	for(index_t p = 0; p < M; p++)
	{
		patch_t & rp = patches[p];
		for(index_t i = 0; i < rp.n; i++)
		{
			os<<rp.xyz.col(i)[0]<<" "<<rp.xyz.col(i)[1]<<" "<<rp.xyz.col(i)[2]<<" ";
		}
		os<<endl;
	}

	os.close();
}

void partial_mesh_reconstruction(size_t old_n_vertices, che * mesh, size_t M, vector<patch_t> & patches, vector<patches_map_t> & patches_map, a_mat & A, a_mat & alpha)
{
	#pragma omp parallel for
	for(index_t p = M; p < patches.size(); p++)
	{
		patch_t & rp = patches[p];

		if(rp.indexes)
		{
			a_vec x = rp.phi * A * alpha.col(p);

			rp.xyz.row(2) = x.t();

			rp.itransform();
		}
	}

	real_t h = 2;

	a_vec V(3);

	#pragma omp parallel for private(V)
	for(index_t v = old_n_vertices; v < mesh->n_vertices(); v++)
	{
		if(patches_map[v].size())
			V = non_local_means_vertex(alpha, v, patches, patches_map, h);
		else
		{
			V(0) = mesh->gt(v).x;
			V(1) = mesh->gt(v).y;
			V(2) = mesh->gt(v).z;
		}

		mesh->get_vertex(v) = *((vertex *) V.memptr());
	}

}

a_vec non_local_means_vertex(a_mat & alpha, const index_t & v, vector<patch> & patches, vector<vpatches_t> & patches_map, const real_t & h)
{
	a_vec n_a_vec(3, arma::fill::zeros);
	real_t sum = 0;

	real_t * w = new real_t[patches_map[v].size()];

	index_t i = 0;
	real_t d = 0;

	for(auto p: patches_map[v])
	{
		d = 0;
		for(auto q: patches_map[v])
			d += norm(alpha.col(p.first) - alpha.col(q.first));
		d /= patches_map[v].size();

		w[i] = exp(- d * d / h);
		sum += w[i];
		i++;
	}

	i = 0;
	for(auto p: patches_map[v])
	{
		w[i] /= sum;
		n_a_vec += w[i] * patches[p.first].xyz.col(p.second);
		i++;
	}

	delete [] w;

	return n_a_vec;
}

[[deprecated]]
void mesh_reconstruction(che * mesh, size_t M, vector<patch_t> & patches, vector<patches_map_t> & patches_map, a_mat & A, a_mat & alpha, const index_t & v_i)
{
	a_mat V(3, mesh->n_vertices(), arma::fill::zeros);

	#pragma omp parallel for
	for(index_t p = 0; p < M; p++)
	{
		patch_t & rp = patches[p];

		if(rp.phi.n_rows)
		{
			a_vec x = rp.phi * A * alpha.col(p);

			rp.xyz.row(2) = x.t();
			rp.itransform();
		}
	}

	real_t h = 0.2;

	#pragma omp parallel for
	for(index_t v = v_i; v < mesh->n_vertices(); v++)
	{
		if(patches_map[v].size())
			V.col(v) = non_local_means_vertex(alpha, v, patches, patches_map, h);
		else
		{
			V(0, v) = mesh->gt(v).x;
			V(1, v) = mesh->gt(v).y;
			V(2, v) = mesh->gt(v).z;
		}
	}

	
	// ------------------------------------------------------------------------

	vertex * new_vertices = (vertex *) V.memptr();

	real_t error = 0;
	#pragma omp parallel for reduction(+: error)
	for(index_t v = v_i; v < mesh->n_vertices(); v++)
		error += *(new_vertices[v] - mesh->get_vertex(v));

	gproshan_debug_var(mesh->n_vertices());
	error /= mesh->n_vertices();
	gproshan_debug_var(error);

	gproshan_debug_var(v_i);
	mesh->set_vertices(new_vertices + v_i, mesh->n_vertices() - v_i, v_i);
}

[[deprecated]]
a_vec non_local_means_vertex(a_mat & alpha, const index_t & v, vector<patch_t> & patches, vector<patches_map_t> & patches_map, const real_t & h)
{
	a_vec n_a_vec(3, arma::fill::zeros);
	real_t sum = 0;

	real_t * w = new real_t[patches_map[v].size()];

	index_t i = 0;
	real_t d = 0;

	for(auto p: patches_map[v])
	{
		d = 0;
		for(auto q: patches_map[v])
			d += norm(alpha.col(p.first) - alpha.col(q.first));
		d /= patches_map[v].size();

		w[i] = exp(- d * d / h);
		sum += w[i];
		i++;
	}

	for(auto p: patches_map[v])
	{
		w[i] /= sum;
		n_a_vec += w[i] * patches[p.first].xyz.col(p.second);
		i++;
	}

	delete [] w;

	return n_a_vec/i;
}

a_vec simple_means_vertex( const index_t & v, vector<patch> & patches, vector<vpatches_t> & patches_map)
{
	a_vec n_a_vec(3, arma::fill::zeros);

	for(auto p: patches_map[v])
		n_a_vec += patches[p.first].xyz.col(p.second);
	return n_a_vec/ patches_map[v].size();
}


} // namespace gproshan::mdict

