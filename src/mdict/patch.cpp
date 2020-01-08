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

#include <queue>

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

void patch::init(che * mesh, const index_t & v, const size_t & n_toplevels, const distance_t & radio, index_t * _toplevel)
{
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	gather_vertices(mesh, v, radio, toplevel);

	if(!_toplevel) delete [] toplevel;
}	

void patch::init_disjoint(che * mesh, const index_t & v, const size_t & n_toplevels, vector<index_t> & _vertices, index_t * _toplevel)
{
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	//vertices = _vertices;
	vertices = std::move(_vertices);

	if(!_toplevel) delete [] toplevel; // If it is null
}

void patch::init_radial_disjoint(che * mesh, const distance_t & radio, const index_t & v, const size_t & n_toplevels, vector<index_t> & _vertices, index_t * _toplevel)
{
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);

	for(index_t i = 0; i < vertices.size(); i++)
	{
		vertex n = mesh->normal(v);// normal at the center
		vertex p = mesh->get_vertex(vertices[i]); 
		vertex c = mesh->get_vertex(v); // central vertices

		p = p - c ;
		p = p - ((p,n)*n);
		if(*p <= radio)	_vertices.push_back(vertices[i]);

	}

	vertices = std::move(_vertices);
	if(!_toplevel) delete [] toplevel; // If it is null
}


// xyz = E.t * (xyz - avg)
void patch::transform()
{
	xyz.each_col() -= x;
	xyz = T.t() * xyz;
//	xyz.row(2).zeros();
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
		//p idx patche where belongs to
			//j: local index 
			//i: global index
			//if(vpatches[vertices[i]].size() == 0)
			vpatches[vertices[i]].push_back({p, j++});
		}
	}
}

void patch::reset_xyz_disjoint(che * mesh, distance_t * dist, size_t M, vector<vpatches_t> & vpatches, const index_t & p,  const fmask_t & mask)
{
	size_t m = vertices.size();
	if(mask)
	{
		m = 0;
		for(index_t i = 0; i < vertices.size(); i++)
			if(mask(i)) { dist[vertices[i] ] = float(p + 1) / M; m++;  } else {dist[vertices[i] ] = INFINITY; };
		
		gproshan_debug(number vertices considered);
		gproshan_debug_var(m);
		gproshan_debug(number vertices masked);
		gproshan_debug_var(vertices.size() - m);
	}
	
	xyz.set_size(3, m);
	for(index_t  j = 0, i = 0; i < vertices.size(); i++)
	{
		if(!mask || mask(i))
		{
			const vertex & v = mesh->gt(vertices[i]);
			xyz(0, j) = v.x;
			xyz(1, j) = v.y;
			xyz(2, j) = v.z;
			
			vpatches[vertices[i]].push_back({p, j++});
		}
	}
}

const a_vec patch::normal()
{
	return T.col(2);
}

void patch::save(const real_t & radio, const size_t & imsize, CImgList<real_t> & imlist)
{
	// Create images with the patches info

	//building the grid
	CImg<real_t> img(imsize, imsize);
	size_t x, y;
	img.fill(0);
	// for each x y plus 1, multiply by delta and floor, get i and j 
	for(index_t i = 0; i < vertices.size(); i++)
	{
		x = floor((xyz.col(i)[0] + radio) * (imsize - 1) / (2 * radio));
		y = floor((xyz.col(i)[1] + radio) * (imsize - 1) / (2 * radio));
		img(x,y) = xyz.col(i)[2];
	}
	
	img.resize(128, 128);
	imlist.insert(img.normalize(0, 255));
	//img.save("tmp/images/test_image.jpg");
	
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

void patch::gather_vertices(che * mesh, const index_t & v, const distance_t & radio, index_t * toplevel)
{
	assert(x.n_elem == 3 && T.n_rows == 3 && T.n_cols == 3);
	
	if(vertices.size()) vertices.clear();
	vertices.reserve(expected_nv);
	
	priority_queue<pair<distance_t, index_t> > qvertices;

	memset(toplevel, -1, sizeof(index_t) * mesh->n_vertices());
	
	a_vec p(3);
	link_t link;

	toplevel[v] = 0;
	qvertices.push({0, v});

	while(!qvertices.empty())
	{
		index_t v = qvertices.top().second;
		qvertices.pop();
		
		vertices.push_back(v);
		
		mesh->link(link, v);
		for(const index_t & he: link)
		{
			const index_t & u = mesh->vt(he);
			if(toplevel[u] == NIL)
			{
				p(0) = mesh->gt(u).x;
				p(1) = mesh->gt(u).y;
				p(2) = mesh->gt(u).z;
				p = T.t() * (p - x);

				toplevel[u] = toplevel[v] + 1;
				
				if(norm(p) < radio)
					qvertices.push({-norm(p), u});
			}
			link.clear();
		}
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

real_t patch::get_min_z()
{
	return xyz.row(2).min();
}

real_t patch::get_max_z()
{
	return xyz.row(2).max();
}

void patch::update_heights(real_t & min, real_t & max, bool flag)
{
	real_t tmp;
	if(flag)
	{
		for(index_t i = 0; i < xyz.n_cols; i++)
			{
				xyz(2, i) = (xyz(2, i) - min) / (max - min);
			}
	}	
	else
	{
		for(index_t i = 0; i < vertices.size(); i++)
		{
			tmp = xyz.col(i)[2];
			tmp = (max - min) * tmp + min;
			xyz.col(i)[2] = tmp;
		}
	}
	
}

void patch::save_z(ostream & os)
{
	index_t i;
	for( i = 0; i < vertices.size()-1; i++)
	{
		os<<xyz.col(i)[2]<<"\t";
	}
	os<<xyz.col(i)[2]<<"\n";
}


} // namespace gproshan::mdict

