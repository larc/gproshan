#include <gproshan/mdict/patch.h>

#include <gproshan/mdict/msparse_coding.h>
#include <gproshan/geodesics/geodesics.h>

#include <queue>
#include <random>

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


size_t patch::expected_nv = 3 * msparse_coding::T * (msparse_coding::T + 1);
real_t patch::nyquist_factor = 0.5;

void patch::init(che * mesh, const index_t & v, const size_t & n_toplevels, const real_t & radio_, index_t * _toplevel)
{
	radio = radio_;
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices];

	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	gather_vertices(mesh, v, radio_, toplevel);

	if(!_toplevel) delete [] toplevel;
}

void patch::init_disjoint(che * mesh, const index_t & v, const size_t & n_toplevels, std::vector<index_t> & _vertices, index_t * _toplevel)
{
	radio = 1;
	index_t * toplevel = _toplevel ? _toplevel : new index_t[mesh->n_vertices];

	gather_vertices(mesh, v, n_toplevels, toplevel);
	jet_fit_directions(mesh, v);
	//vertices = _vertices;
	vertices = std::move(_vertices);

	if(!_toplevel) delete [] toplevel; // If it is null
}

bool patch::exists(index_t idx)
{
	for(size_t i=1; i < vertices.size(); ++i)
	{
		if(vertices[i] == idx)
			return true;
	}
	return false;
}

index_t patch::find(const index_t * indexes, size_t nc, index_t idx_global)
{
	for(size_t i=0; i<nc; ++i)
		if(indexes[i] == idx_global) return i;

	return -1;
}

bool patch::add_vertex_by_trigs(vertex & n, std::vector<vertex> & N, double thr_angle, const real_t * geo, che * mesh, const index_t & v, real_t & area, real_t & proj_area, real_t deviation)
{
	// it needs to return both vertices
	// it needs to filter repeated indexes.
	// p should be the maximun

	index_t a, b, i = 0;
	vertex min_he;
	double area_face = 0, proj_area_face = 0;
	double angle = M_PI;
	double tmp_angle;
	bool added = false;
	vertex pav, pbv, va, vb,vv;

	for(const index_t & he: mesh->star(v))
	{
		a = mesh->halfedge(he_next(he)); //index of the next vertex index_t
		b = mesh->halfedge(he_prev(he));
		va = mesh->point(a);
		vb = mesh->point(b);
		vv = mesh->point(v);
		// If is an adjacent face
		assert(a < mesh->n_vertices);
		assert(b < mesh->n_vertices);
		assert(v < mesh->n_vertices);

		if(geo[a] < geo[v] || geo[b] < geo[v] )
		{
			if(geo[a] < geo[v])
				i = find(vertices.data(), vertices.size(), a);
			else
				i = find(vertices.data(), vertices.size(), b);

			tmp_angle = acos(dot(mesh->normal_he(he), N[i]));

			if(angle > tmp_angle && tmp_angle < thr_angle && acos(dot(mesh->normal_he(he), N[0])) < deviation) // Fullfill conditions
			{
				angle = tmp_angle;
				area_face = mesh->area_trig(he / 3);

				// compute projected area
				pav = va - vv + (dot(n, vv) - dot(n, va)) * n;
				pbv = vb - vv + (dot(n, vv) - dot(n, vb)) * n;
				proj_area_face = norm(cross(pav, pbv)) / 2;

				min_he = mesh->normal_he(he);
				added = true;
			}
		}
	}
	//p = mesh->point(indexes[i]);
	//p = p - c ;
	//p = p - ((p,n)*n);

	area += area_face;
	proj_area += proj_area_face;

	if(added)
	{
		vertices.push_back(v);
		N.push_back(min_he);
	}

	return added;
}

void patch::init_random(const vertex & c, const a_mat & T, const real_t & radio, const real_t & max_radio, const real_t & percent, const real_t & fr)
{
	this->radio = radio;
	this->T = T;

	x.resize(3);
	x(0) = c.x();
	x(1) = c.y();
	x(2) = c.z();


	std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0, 1);
	//std::normal_distribution<double> dis(0,0.5);

	// free the parameters to the interface
	// fix the point cloud viewer
	size_t n_points = (radio / max_radio) * percent; // change this using a sigmoid function
	xyz.resize(3, n_points);

	xyz(0, 0) = 0;
	xyz(1, 0) = 0;
	xyz(2, 0) = 0;

	for(size_t i = 1; i < n_points; ++i)
	{
		double a = abs(dis(gen)) * 2 * M_PI;
		double r = fr * abs(dis(gen));

		xyz(0, i) = r * cos(a);
		xyz(1, i) = r * sin(a);
		xyz(2, i) = 0;
	}
}

void patch::recover_radial_disjoint(che * mesh, const real_t & radio_, const index_t & v)
{
	// for small meshes 6000 0.e-5
	// for others 2.e-5
	geodesics::params params;
	params.radio = radio_ + 1e-5;//numeric_limits<real_t>::epsilon();

	geodesics geo(mesh, {v}, params);

	index_t * indexes = new index_t[geo.n_sorted_index()];
	geo.copy_sorted_index(indexes, geo.n_sorted_index());

	vertices.push_back(v);

	for(index_t i=1; i<geo.n_sorted_index(); ++i)
	{
		vertices.push_back(indexes[i]);
	}

	size_t d_fitting = 2;
	vertex p, c, n;
	size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	if(vertices.size() > min_points)
	{
		jet_fit_directions(mesh, v);
		n.x() = T(0, 2); n.y() = T(1, 2); n.z() = T(2, 2);
		radio = -INFINITY;

		for(index_t i=1; i < vertices.size(); ++i)
		{
			p = mesh->point(indexes[i]);
			c = mesh->point(v); // central vertices

			p = p - c ;
			p = p - dot(p, n) * n;

			if(norm(p) > radio)
			{
				radio = norm(p);
			}
		}
	}
	else
	{
		gproshan_debug_var(vertices.size());
	}

}

void patch::init_radial_disjoint(	real_t & euc_radio,
									real_t & geo_radio,
									che * mesh,
									const index_t & v,
									const real_t & delta,
									const real_t & sum_thres,
									const real_t & area_thres,
									const real_t & area_mesh
									)
{
	radio = -INFINITY;
	min_nv = 128;

	euc_radio = -INFINITY;

	normal_fit_directions(mesh, v);

	a_vec vn = T.col(2);
	vertex n = { vn(0), vn(1), vn(2) };

	vertices.push_back(v);

	std::vector<vertex> N;
	N.push_back(n);

	real_t area = 0;
	real_t proj_area = std::numeric_limits<real_t>::epsilon();
	real_t ratio;

	vertex c = mesh->point(v);

	geodesics::params params;
	params.dist_alloc = new real_t[mesh->n_vertices];
	params.fun = [&](const index_t & u) -> bool
	{
		if(u == v) return true;

		ratio = area / proj_area;

		if(add_vertex_by_trigs(n, N, delta, params.dist_alloc, mesh, u, area, proj_area, M_PI / 2.5 ) && (ratio < sum_thres || (area / area_mesh) < area_thres) )
		{
			euc_radio = std::max(euc_radio, norm(mesh->point(u) - c));
			return true;
		}

		return false;
	};

	geodesics geo(mesh, {v}, params);
	delete [] params.dist_alloc;

	// Refit the points and update the radius
	size_t d_fitting = 2;
	size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	if(vertices.size() > min_points)
		jet_fit_directions(mesh, v);
	else
		normal_fit_directions(mesh,v);

	n.x() = T(0, 2); n.y() = T(1, 2); n.z() = T(2, 2);
	radio = -INFINITY;

	vertex p;
	for(auto & vi: vertices)
	{
		p = mesh->point(vi);

		p = p - c ;
		p = p - dot(p, n) * n;

		radio = std::max(radio, norm(p));
	}

	geo_radio = geo[vertices.back()];
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

void patch::reset_xyz(che * mesh, std::vector<vpatches_t> & vpatches, const index_t & p, const fmask_t & mask)
{
	size_t m = vertices.size();

	if(mask)
	{
		m = 0;
		for(index_t i = 0; i < vertices.size(); ++i)
			if(mask(vertices[i])) ++m;
	}

	xyz.set_size(3, m);
	for(index_t j = 0, i = 0; i < vertices.size(); ++i)
	{
		if(!mask || mask(vertices[i]))
		{
			const vertex & v = mesh->point(vertices[i]);
			xyz(0, j) = v.x();
			xyz(1, j) = v.y();
			xyz(2, j) = v.z();
		//p idx patche where belongs to
			//j: local index
			//i: global index
			//if(vpatches[vertices[i]].size() == 0)
			//vpatches[vertices[i]].push_back({p, j++});
			vpatches[vertices[i]][p] = j++;
		}
	}



}

double area_tri(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}

void patch::remove_extra_xyz_disjoint(size_t & max_points)
{
	if(vertices.size() > max_points)
	{
		arma::uvec xi;
		xi.zeros(max_points);
		for (size_t i=1; i< max_points; ++i) xi[i] = i;
		xi = arma::shuffle(xi);
		xyz = xyz.cols(xi);
	}

}
void patch::add_extra_xyz_disjoint(che * mesh, std::vector<vpatches_t> & vpatches, const index_t & p)
{

	size_t m = std::max (vertices.size(), min_nv);


	size_t j = vertices.size();
	std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0, 1);

	//if(p == 76)	gproshan_debug_var(xyz);

	while(j < m )
	{

		// add new vertices
		// create a random point
		real_t a = abs(dis(gen)) * 2 * M_PI;
		real_t r = abs(dis(gen));
		a_vec np = { r * cos(a), r * sin(a), 0 };

		//gproshan_debug_var(np);
		// find the closest point
		index_t min_v;
		double min_d = INFINITY;
		for(index_t v: vertices)
		{
			a_vec aux = xyz.col(vpatches[v][p]);
			aux(2) = 0;

			if(norm(np - aux) < min_d)
			{
				min_d = norm(np - aux);
				min_v = v;
			}
		}

		// forstar to find closest trinagle
		a_mat abc(3,3);
		for(const index_t & he: mesh->star(min_v))
		{
			//discard triangles outside the patch
			vpatches_t & ma = vpatches[mesh->halfedge(he_next(he))];
			vpatches_t & mb = vpatches[mesh->halfedge(he_prev(he))];

			if(ma.find(p) != ma.end() && mb.find(p) != mb.end())
			{
				arma::uvec xi = { vpatches[min_v][p], ma[p], mb[p] };
				abc = xyz.cols(xi);

				//gproshan_debug_var(np);
				// verify if this is inside a triangle

				double A = area_tri(abc(0,0),	abc(1,0),	abc(0,1),	abc(1,1),	abc(0,2),	abc(1,2) );
				double A1 = area_tri(np(0),		np(1),		abc(0,1),	abc(1,1),	abc(0,2),	abc(1,2) );
				double A2 = area_tri(abc(0,0), 	abc(1,0),	np(0),		np(1),		abc(0,2),	abc(1,2) );
				double A3 = area_tri(abc(0,0), 	abc(1,0),	abc(0,1),	abc(1,1),	np(0),		np(1) );



				if(abs(A - (A1 + A2 + A3)) < std::numeric_limits<real_t>::epsilon())
				{
					a_mat proj_abc = abc.tail_cols(2).each_col() - abc.col(0);
					np -= abc.col(0);

					a_vec coef = arma::inv(proj_abc.head_rows(2)) * np.head(2);
					np = proj_abc * coef + abc.col(0);

					if(!isnan(np(2)))
					{
						xyz(0, j) = np(0);
						xyz(1, j) = np(1);
						xyz(2, j) = np(2);
						++j;
						/*if(p == 76)
						{
							gproshan_debug_var(np);
							gproshan_debug_var(abc);
							gproshan_debug_var(coef);
						}*/
						break;
					}
					/*if(!isnan(z)) //z = 0;*/

				}
			}

		}
	}
}

void patch::reset_xyz_disjoint(che * mesh, real_t * dist, size_t M, std::vector<vpatches_t> & vpatches, const index_t & p, const fmask_t & mask)
{
	size_t m = vertices.size();
	if(mask)
	{
		m = 0;
		for(index_t i = 0; i < vertices.size(); ++i)
			if(mask(i))
			{
				dist[vertices[i]] = float(p + 1) / M;
				++m;
			}
			else
			{
				dist[vertices[i]] = INFINITY;
			};
		/*gproshan_debug(number vertices considered);
		gproshan_debug_var(m);
		gproshan_debug(number vertices masked);
		gproshan_debug_var(vertices.size() - m);*/
	}

	m = std::max(vertices.size(), min_nv);
	xyz.set_size(3, m);

	index_t i = 0;
	for(auto & vi: vertices)
	{
		if(!mask || mask(vi))
		{
			const vertex & v = mesh->point(vi);
			xyz(0, i) = v.x();
			xyz(1, i) = v.y();
			xyz(2, i) = v.z();

			vpatches[vi][p] = i++;
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
	xyz = xyz / factor;
}

const a_vec patch::normal()
{
	return T.col(2);
}

void patch::gather_vertices(che * mesh, const index_t & v, const size_t & n_toplevels, index_t * toplevel)
{
	if(vertices.size()) vertices.clear();

	vertices.reserve(expected_nv);
	memset(toplevel, -1, sizeof(index_t) * mesh->n_vertices);

	toplevel[v] = 0;
	vertices.push_back(v);

	for(index_t i = 0; i < vertices.size(); ++i)
	{
		const index_t & v = vertices[i];
		if(toplevel[v] == n_toplevels)
			break;

		for(const index_t & u: mesh->link(v))
			if(toplevel[u] == NIL)
			{
				vertices.push_back(u);
				toplevel[u] = toplevel[v] + 1;
			}
	}
}

void patch::gather_vertices(che * mesh, const index_t & v, const real_t & radio, index_t * toplevel)
{
	assert(x.n_elem == 3 && T.n_rows == 3 && T.n_cols == 3);

	if(vertices.size()) vertices.clear();
	vertices.reserve(expected_nv);

	std::priority_queue<std::pair<real_t, index_t> > qvertices;

	memset(toplevel, -1, sizeof(index_t) * mesh->n_vertices);

	a_vec p(3);

	toplevel[v] = 0;
	qvertices.push({0, v});

	while(!qvertices.empty())
	{
		index_t v = qvertices.top().second;
		qvertices.pop();

		vertices.push_back(v);

		for(const index_t & u: mesh->link(v))
		{
			if(toplevel[u] == NIL)
			{
				p(0) = mesh->point(u).x();
				p(1) = mesh->point(u).y();
				p(2) = mesh->point(u).z();
				p = T.t() * (p - x);

				toplevel[u] = toplevel[v] + 1;

				if(norm(p) < radio)
					qvertices.push({-norm(p), u});
			}
		}
	}
}

/// Compute the principal directions of the patch, centering in the vertex \f$v\f$.
/// See: https://doc.cgal.org/latest/Jet_fitting_3/index.html
void patch::jet_fit_directions(che * mesh, const index_t & v)
{
	size_t d_fitting = 2;
	size_t d_monge = 2;
	//size_t min_points = (d_fitting + 1) * (d_fitting + 2) / 2;
	//assert(vertices.size() > min_points);

	std::vector<DPoint> in_points;
	in_points.reserve(vertices.size());
	for(const index_t & u: vertices)
		in_points.push_back(DPoint(mesh->point(u).x(), mesh->point(u).y(), mesh->point(u).z()));

	My_Monge_form monge_form;
	My_Monge_via_jet_fitting monge_fit;
	monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

	vertex normal = mesh->normal(v);
	monge_form.comply_wrt_given_normal(DVector(normal.x(), normal.y(), normal.z()));

	x.set_size(3);
	x(0) = mesh->point(v).x();
	x(1) = mesh->point(v).y();
	x(2) = mesh->point(v).z();

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
	x(0) = mesh->point(v).x();
	x(1) = mesh->point(v).y();
	x(2) = mesh->point(v).z();


	vertex nz = mesh->normal(v);
	vertex nx = mesh->vertex_he(he_next(mesh->evt(v)));
//	GT[VT[he_next(EVT[v]]]
	vertex c = mesh->point(v);
	vertex ny;
	nx = nx - c ;
	nx = nx - dot(nx, nz) * nz;

	ny = cross(nz, nx);
	nx = normalize(nx);
	ny = normalize(ny);
	nz = normalize(nz);


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
		for(index_t i = 0; i < xyz.n_cols; ++i)
			{
				xyz(2, i) = (xyz(2, i) - min) / (max - min);
			}
	}
	else
	{
		for(index_t i = 0; i < vertices.size(); ++i)
		{
			tmp = xyz.col(i)[2];
			tmp = (max - min) * tmp + min;
			xyz.col(i)[2] = tmp;
		}
	}

}

void patch::save_z(std::ostream & os)
{
	index_t i;
	for( i = 0; i < vertices.size()-1; ++i)
	{
		os<<xyz.col(i)[2]<<"\t";
	}
	os<<xyz.col(i)[2]<<"\n";
}

void patch::compute_avg_distance(che * mesh, std::vector<vpatches_t> & vpatches, const index_t & p)
{
	avg_dist = INFINITY;
	std::vector<double> distances;

	for(size_t i = 0; i < vertices.size(); ++i)
	{
		for(const index_t & u: mesh->link(vertices[i]))
		{
			for(auto itp: vpatches[u])
			{
				if(itp.first == p)
				{
					a_vec a = xyz.col(i);
					a_vec b = xyz.col(itp.second);
					a(2) = 0;
					b(2) = 0;
					distances.push_back(norm(a - b));
					break;
				}
			}
		}
	}
	/*
		for(size_t j = i+1; j < vertices.size(); ++j) // replace for 1 ring
		{
			a_vec a = xyz.col(i);
			a_vec b = xyz.col(j);
			a(2) = 0;
			b(2) = 0;
			distances.push_back(norm(a - b));
		}
	*/
	std::sort(distances.begin(), distances.end());
	size_t n_elem = distances.size();
	if(distances.size()%2 ==0)
	{
		avg_dist = (distances[n_elem/2] + distances[(n_elem/2 -1)])/2;
	}
	else
	{
		avg_dist = distances[n_elem/2];
	}
			//avg_dist = avg_dist + norm(xyz.col(i)- xyz.col(j));
}

bool patch::is_covered( bool * covered)
{
	for(index_t i = 0; i < vertices.size(); ++i)
		if(!covered[i]) return false;

	return true;
}


} // namespace gproshan::mdict

