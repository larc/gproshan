#include "patch.h"

#include "dictionary.h"
#include "che_sphere.h"
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
#define PI 3.14159265
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

void patch::init(che * mesh, const index_t & v, const size_t & n_toplevels, const distance_t & radio_, index_t * _toplevel)
{
	radio = radio_;
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	gather_vertices(mesh, v, radio_, toplevel);

	if(!_toplevel) delete [] toplevel;
}	

void patch::init_disjoint(che * mesh, const index_t & v, const size_t & n_toplevels, vector<index_t> & _vertices, index_t * _toplevel)
{
	radio = 1;
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices()];
	
	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	//vertices = _vertices;
	vertices = std::move(_vertices);

	if(!_toplevel) delete [] toplevel; // If it is null
}

bool  patch::exists(index_t idx)
{
	for(size_t i=1; i < vertices.size(); i++)
	{
		if(vertices[i] == idx)
			return true;
	}
	return false;
}

index_t patch::find(index_t * indexes, size_t nc, index_t idx_global)
{
	for(size_t i=0; i<nc; i++)
		if(indexes[i] == idx_global) return i;
	return -1;
}
bool patch::add_vertex_by_faces(vector<vertex> & N, index_t * indexes, size_t nc, double thr_angle, const geodesics & geo, che * mesh, const index_t & v, double &sum, double deviation)
{

	index_t a, b, i = 0;
	vertex min_he;
	double angle = PI;
	double tmp_angle;
	bool added = false;

	for_star(he, mesh, v)
	{
	
		a = mesh->vt(next(he)); //index of the next vertex index_t
		b = mesh->vt(prev(he)); 
		/*gproshan_debug_var(geo[b]);
		gproshan_debug_var(geo[a]);
		gproshan_debug_var(geo[v ]);
*/
		// If is an adjacent face
		if( geo[a] < geo[v] || geo[b] < geo[v] )
		{
			if(geo[a] < geo[v])
			{
				i = find(indexes, nc,a);
				//gproshan_debug_var(a);
				//gproshan_debug_var(i); 
			}
			else
			{
				i = find(indexes, nc, b); 
				//gproshan_debug_var(b);
				//gproshan_debug_var(i);
			}
			tmp_angle = acos( (mesh->normal_he(he), N[i]) );
		/*	gproshan_debug_var(tmp_angle);
			gproshan_debug_var(angle);
			gproshan_debug_var(thr_angle);
			gproshan_debug_var(N[i]);*/
			//gproshan_debug_var(mesh->normal_he(he));

			if ( angle >  tmp_angle  && tmp_angle < thr_angle &&  acos( (mesh->normal_he(he), N[0]) ) < deviation ) // Fullfill conditions
			{
				angle = tmp_angle;
				
				min_he = mesh->normal_he(he); 
				if( !exists(v) ) vertices.push_back(v);
				added = true;
			}
				
		}
		//gproshan_debug_var(N[i]);
	
	}
	//gproshan_debug_var(min_he);
	sum +=  acos( (min_he, N[i]) );

	N.push_back(min_he);
	return added;
}


void patch::init_radial_disjoint(che * mesh, const distance_t & radio_, const index_t & v, distance_t & euc_radio)
{

	//radio = radio_;
	radio = -INFINITY;
	//che_sphere my_sphere(1,12);
	//string sphere_file = tmp_file_path("sphere");
	//che_off::write_file(&my_sphere, sphere_file);

	normal_fit_directions(mesh, v);

	geodesics geo(mesh, {v}, geodesics::FM,  NULL, false, 0, radio_);
	index_t * indexes = new index_t[geo.n_sorted_index()];
	geo.copy_sorted_index(indexes, geo.n_sorted_index());

	a_vec vn = T.col(2);// normal at the center
	vertex n;
	n.x = vn(0); n.y = vn(1); n.z = vn(2);
	vertex p, c;
	vertices.push_back(v);
	euc_radio = -INFINITY;
	
	vector<vertex> N;
	N.push_back(n);
	//double angle;
	double sum_angle = 0;
	double delta = 0;

	//vertex prev_n, curr_n;
	//gproshan_debug_var(geo.n_sorted_index());

	for(index_t i=1; i<geo.n_sorted_index(); i++)
	{
		
		p = mesh->get_vertex(indexes[i]); 
		c = mesh->get_vertex(v); // central vertices

		p = p - c ;
		p = p - ((p,n)*n);
		//if(*p <= radio )
		//gproshan_debug_var(indexes[i]);
	//	sum_angle = acos( (n, mesh->normal(indexes[i]) ) );
		//angle = acos( (n, mesh->normal(indexes[i]) ) ) ;
		
		//prev_n = curr_n;
		//curr_n = normal_trim(geo, mesh, indexes[i]);
		//gproshan_debug_var(acos( (prev_n, curr_n) ));
		//if( angle < PI/2.5 && (sum_angle) <= delta * PI)
		//penalize gowing, I want them to grow only if they do not vary so much
		//if( angle < PI/2.5 && acos( (mesh->normal(indexes[i-1]), mesh->normal(indexes[i]) ) ) <= PI/8) // find borders
	
		// add one new candidate vertex //first regulates variation, // second regulates size of the patch
		if( add_vertex_by_faces(N, indexes, geo.n_sorted_index(), PI/6, geo, mesh, indexes[i], sum_angle, PI/2.5 ) && sum_angle < PI  )
		// pi is too much
		//if( angle < PI/2.5 && acos( (prev_n, curr_n) ) <= PI/7) // find borders
		{
			//gproshan_debug_var(vertices.size());
			
			if(*p > radio)
			{
				radio = *p;
			}
			//compute euclidean radio
			p = mesh->get_vertex(indexes[i]);
			if(*(p - c) > euc_radio)
				euc_radio = *(p - c);
			//gproshan_debug_var(euc_radio);
			//vertices.push_back(indexes[i]);
			//gproshan_debug_var(geo[indexes[i]]);
			//sum_angle += angle;
		//	delta += 0.001	;
			// sharp meshes 0.001
			// smooth meshes 0.035 at max
		}
		else
		{
			break;
		}
	//	gproshan_debug_var(acos( (mesh->normal(indexes[i-1]), mesh->normal(indexes[i]) ) ));
	}
	// Refit the points and update the radius
	size_t d_fitting = 2;
	size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	if(vertices.size() > min_points)
	{
		jet_fit_directions(mesh, v);
		n.x  = T(0, 2); n.y  = T(1, 2); n.z  = T(2, 2);
		radio = -INFINITY;

		for(index_t i=1; i < vertices.size(); i++)
		{
			p = mesh->get_vertex(indexes[i]); 
			c = mesh->get_vertex(v); // central vertices

			p = p - c ;
			p = p - ((p,n)*n);

			if(*p > radio)
			{
				radio = *p;
			}

		}

	}

	 
	//gproshan_debug_var(sum_angle);
/*	gproshan_debug_var(PI/2.5);
	gproshan_debug_var(sum_angle);
	gproshan_debug_var(PI/(delta+0.05)*(vertices.size()-1));*/
	delete indexes;
	//gproshan_debug_var(v);

}


void patch::update_radial_disjoint(che * mesh, const index_t & v, vector<index_t> & _vertices)
{
	a_vec vn = T.col(2);// normal at the center
	vertex n;
	n.x = vn(0); n.y = vn(1); n.z = vn(2);
	vertex p, c;

	if(vertices.size() >= 7)
	{
		radio = 0;
		vertices = std::move(_vertices);
		for(size_t i=1; i < vertices.size(); i++)
		{
			p = mesh->get_vertex(vertices[i]); 
			c = mesh->get_vertex(v); // central vertices

			p = p - c ;
			p = p - ((p,n)*n);
			if(*p > radio)
			{
				radio = *p;
			}
		}
		jet_fit_directions(mesh, v);
		
	}
	else
	{
		gproshan_debug_var(v);
		gproshan_debug_var(_vertices.size());
	}
	
}

void patch::init_curvature_growing(che * mesh, const index_t & v, a_mat & normals)
{
	radio = -INFINITY;
	geodesics geo(mesh, {v}, geodesics::FM,  NULL, false, mesh->n_vertices());
	index_t * indexes = new index_t[geo.n_sorted_index()];
	geo.copy_sorted_index(indexes, geo.n_sorted_index());
	a_vec vn = normals.col(v);
	vertex n = mesh->normal( v ) ;
	//n.x = vn(0); n.y = vn(1); n.z = vn(2);
//	gproshan_debug_var(geo.n_sorted_index());

	vertices.push_back(v);

	for(size_t i=1; i<geo.n_sorted_index(); i++)
	{
	//	gproshan_debug_var(geo[indexes[i]]);
		// normal at the center
	
		vertex p = mesh->get_vertex(indexes[i]); 
		vertex c = mesh->get_vertex(v); // central vertices

		//gather the good ones
		/*gproshan_debug_var(i);
		gproshan_debug_var(acos((n, mesh->normal(indexes[i]))) );
		gproshan_debug_var(PI/2);*/
		//
		if( acos((n, mesh->normal(indexes[i]) )) <= PI/2 )		
		{
			vertices.push_back(indexes[i]);
		}
		else
		{
			break;
		}
			

	}
	
	//gproshan_debug_var(vertices.size());
	size_t d_fitting = 2;
	size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	if(vertices.size() <= min_points )
	{
		vertices.clear();
	}
	else
	{
		jet_fit_directions(mesh, v);
		size_t n_vertices = vertices.size();
		vertices.clear();
		vertices.push_back(v);

		vn = T.col(2);// normal at the center
		n.x = vn(0); n.y = vn(1); n.z = vn(2);

		for(size_t i=1; i < n_vertices; i++)
		{
			vertex p = mesh->get_vertex(indexes[i]); 
			vertex c = mesh->get_vertex(v); // central vertices

			p = p - c ;
			p = p - ((p,n)*n);
		
			if( acos( (n, mesh->normal(indexes[i]) ) ) <= PI/2 )	
			//if( (n, mesh->normal(indexes[i])) >= 0 )	// zerear los otroso
			{
				if(*p > radio)
					radio = *p;
				vertices.push_back(indexes[i]);
			}
			else
			{
				break;
			}
				
		}
	//	gproshan_debug(after);
	//	gproshan_debug_var(vertices.size());
		if(vertices.size() <= min_points )
		{
			vertices.clear();
		}
			

	}
		

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
		/*gproshan_debug(number vertices considered);
		gproshan_debug_var(m);
		gproshan_debug(number vertices masked);
		gproshan_debug_var(vertices.size() - m);*/
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

void patch::scale_xyz(const real_t & radio_f)
{
	real_t factor = radio_f/radio;
	xyz = factor * xyz;

}

void patch::iscale_xyz(const real_t & radio_f)
{
	real_t factor = radio_f/radio;
	xyz =  xyz / factor;
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

void patch::normal_fit_directions(che * mesh, const index_t & v)
{
	x.set_size(3);
	x(0) = mesh->gt(v).x;
	x(1) = mesh->gt(v).y;
	x(2) = mesh->gt(v).z;

	
	vertex nz = mesh->normal(v);
	vertex nx = mesh->gt_vt_next_evt(v);
//	GT[VT[next(EVT[v]]]
	vertex c = mesh->get_vertex(v);
	vertex ny;
	nx = nx - c ;
	nx = nx - ((nx,nz)*nz);

	ny = (nz * nx);
	nx.unit();
	ny.unit();
	nz.unit();


	T.set_size(3, 3);
	T(0, 0) = nx[0];
	T(1, 0) = nx[1];
	T(2, 0) = nx[2];
	T(0, 1) = ny[0];
	T(1, 1) = ny[1];
	T(2, 1) = ny[2];
	T(0, 2) = nz[0];
	T(1, 2) = nz[1];
	T(2, 2) = nz[2];

	T = normalise(T);

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

void patch::compute_avg_distance()
{
	avg_dist = INFINITY;
	vector<double> distances;
	for(size_t i = 0; i < vertices.size(); i++)
		for(size_t j = i+1; j < vertices.size(); j++)
		{
			a_vec a = xyz.col(i);
			a_vec b = xyz.col(j);
			a(2) = 0;
			b(2) = 0;
			distances.push_back(norm(a - b));
		}
	sort(distances.begin(), distances.end());
	size_t n_elem = distances.size();
	if(distances.size()%2 ==0)
	{
		avg_dist = (distances[n_elem/2] + distances[(n_elem/2 -1) ])/2;
	}
	else 
	{	
		avg_dist = distances[n_elem/2];
	}	
			//avg_dist = avg_dist + norm(xyz.col(i)- xyz.col(j));
}

bool patch::is_covered( bool * covered)
{
	for(index_t i = 0; i < vertices.size(); i++)
		if(!covered[i]) return false;
	
	return true;
}

	//avg_dist = avg_dist/vertices.size();
} // namespace gproshan::mdict

