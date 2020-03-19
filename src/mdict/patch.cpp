#include "patch.h"

#include "dictionary.h"

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
typedef Data_Kernel::Vector_3 DVector;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;


size_t patch::expected_nv = 3 * dictionary::T * (dictionary::T + 1);

void patch::init(che * mesh, const index_t & v, const size_t & n_toplevels, const real_t & radio, index_t * _toplevel)
{
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	gather_vertices(mesh, v, radio, toplevel);

	if(!_toplevel) delete [] toplevel;
}	

// xyz = E.t * (xyz - avg)
void patch::transform()
{
	xyz.each_col() -= x;
	xyz = T.t() * xyz;
}

void patch::itransform()
{
	xyz = T * xyz;
	xyz.each_col() += x;
}

void patch::reset_xyz(che * mesh, vector<vpatches_t> & vpatches, const index_t & p, const fmask_t & mask)
{
	size_t m = vertices.size();
	if(mask)
	{
		m = 0;
		for(index_t i = 0; i < vertices.size(); i++)
			if(mask(vertices[i])) m++;
	}

	xyz.set_size(3, m);
	for(index_t j = 0, i = 0; i < vertices.size(); i++)
	{
		if(!mask || mask(vertices[i]))
		{
			const vertex & v = mesh->gt(vertices[i]);
			xyz(0, j) = v.x;
			xyz(1, j) = v.y;
			xyz(2, j) = v.z;

			vpatches[vertices[i]].push_back({p, j++});
		}
	}
}

void patch::gather_vertices(che * mesh, const index_t & v, const size_t & n_toplevels, index_t * toplevel)
{
	if(vertices.size()) vertices.clear();

	vertices.reserve(expected_nv);
	memset(toplevel, -1, sizeof(index_t) * mesh->n_vertices());
	
	link_t link;
	toplevel[v] = 0;
	vertices.push_back(v);
	for(index_t i = 0; i < vertices.size(); i++)
	{
		const index_t & v = vertices[i];
		if(toplevel[v] == n_toplevels)
			break;
		
		mesh->link(link, v);
		for(const index_t & he: link)
		{
			const index_t & u = mesh->vt(he);
			if(toplevel[u] == NIL)
			{
				vertices.push_back(u);
				toplevel[u] = toplevel[v] + 1;
			}
		}

		link.clear();	
	}	
}

void patch::gather_vertices(che * mesh, const index_t & v, const real_t & radio, index_t * toplevel)
{
	assert(x.n_elem == 3 && T.n_rows == 3 && T.n_cols == 3);
	if(vertices.size()) vertices.clear();
	
	vector<index_t> qvertices;
	qvertices.reserve(expected_nv);
	
	vertices.reserve(expected_nv);
	memset(toplevel, -1, sizeof(index_t) * mesh->n_vertices());
	
	size_t count_toplevel = 0;
	size_t current_toplevel = 0;

	a_vec p(3);
	link_t link;
	toplevel[v] = 0;
	qvertices.push_back(v);
	for(index_t i = 0; i < qvertices.size(); i++)
	{
		const index_t & v = qvertices[i];
		p(0) = mesh->gt(v).x;
		p(1) = mesh->gt(v).y;
		p(2) = mesh->gt(v).z;
		p = T.t() * (p - x);
		p(2) = 0;
		
		if(vertices.size() > expected_nv) break;
		if(toplevel[v] != current_toplevel)
		{
			if(count_toplevel == 0) break;
			current_toplevel++;
			count_toplevel = 0;
		}
		if(norm(p) <= radio)
		{
			vertices.push_back(v);
			count_toplevel++;
		}
		
		mesh->link(link, v);
		for(const index_t & he: link)
		{
			const index_t & u = mesh->vt(he);
			if(toplevel[u] == NIL)
			{
				qvertices.push_back(u);
				toplevel[u] = toplevel[v] + 1;
			}
		}

		link.clear();	
	}	
}

/// Compute the principal directions of the patch, centering in the vertex \f$v\f$.
/// See: https://doc.cgal.org/latest/Jet_fitting_3/index.html
void patch::jet_fit_directions(che * mesh, const index_t & v)
{
	size_t d_fitting = 2;
	size_t d_monge = 2;
	size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	assert(vertices.size() > min_points);

	vector<DPoint> in_points;
	in_points.reserve(vertices.size());
	for(const index_t & u: vertices)
		in_points.push_back(DPoint(mesh->gt(u).x, mesh->gt(u).y, mesh->gt(u).z));
	
	My_Monge_form monge_form;
	My_Monge_via_jet_fitting monge_fit;
	monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

	vertex normal = mesh->normal(v);
	monge_form.comply_wrt_given_normal(DVector(normal.x, normal.y, normal.z));
	
	x.set_size(3);
	x(0) = mesh->gt(v).x;
	x(1) = mesh->gt(v).y;
	x(2) = mesh->gt(v).z;

	T.set_size(3, 3);
	T(0, 0) = monge_form.maximal_principal_direction()[0];
	T(1, 0) = monge_form.maximal_principal_direction()[1];
	T(2, 0) = monge_form.maximal_principal_direction()[2];
	T(0, 1) = monge_form.minimal_principal_direction()[0];
	T(1, 1) = monge_form.minimal_principal_direction()[1];
	T(2, 1) = monge_form.minimal_principal_direction()[2];
	T(0, 2) = monge_form.normal_direction()[0];
	T(1, 2) = monge_form.normal_direction()[1];
	T(2, 2) = monge_form.normal_direction()[2];
}


} // namespace gproshan::mdict

