#include "che_fill_hole.h"

#include "che_off.h"
#include "laplacian.h"
#include "viewer/viewer.h"

#include <queue>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


bool operator<(const border_t & a, const border_t & b)
{
	return a.theta > b.theta;
}

a_vec normal_face(const vector<a_vec> & tmp_vertices, const index_t & a_v, const index_t & b_v, const index_t & c_v)
{
	a_vec a = tmp_vertices[c_v] - tmp_vertices[a_v];
	a_vec b = tmp_vertices[b_v] - tmp_vertices[a_v];
	return normalise(cross(a,b));
}

che * mesh_simple_fill_hole(che * mesh, const vector<index_t> & border_vertices, const size_t & max_iter = 1000)
{
	vector<vertex> vertices;
	vertex normal, normal_v, edge_v, v;
	a_mat E(3,3);
	a_vec ve(3);

	vertices.reserve(border_vertices.size());

	debug_me(unfolding hole)

	for(const index_t & b: border_vertices)
	{
		v = mesh->gt(b);
		normal_v = mesh->normal(b);
		edge_v = mesh->gt_vt_next_evt(b) - v;
		edge_v -= (normal_v, edge_v) * normal_v;

		E(0, 2) = normal_v.x;
		E(1, 2) = normal_v.y;
		E(2, 2) = normal_v.z;

		E(0, 0) = edge_v.x;
		E(1, 0) = edge_v.y;
		E(2, 0) = edge_v.z;

		E.col(1) = cross(E.col(2), E.col(0));
		E = normalise(E);

		ve(0) = v.x; ve(1) = v.y; ve(2) = v.z;

		ve = E.t() * ve;
		vertices.push_back(*((vertex *) ve.memptr()));
		viewer::other_vertices.push_back(vertices.back());
	}

	return fill_hole_front_angles(vertices, mesh->mean_edge(), normal, max_iter);
}

che * mesh_fill_hole(che * mesh, const vector<index_t> & border_vertices, const size_t & max_iter, const vector<pair<index_t, index_t> > & split_indices = {})
{
	vector<vertex> vertices[2];
	vector<index_t> merge_vertices[2];
	vertex normal;

	size_t size = border_vertices.size();
	index_t i, j, n_v;
	che * hole = NULL;
	che * aux_hole;

	index_t * vmap_border = new index_t[size];

	index_t c = 1;
	distance_t mean_edge = mesh->mean_edge();

	auto gen_vertices = [&mean_edge](vector<index_t> & merge_vertices, vector<vertex> & vertices, const vertex & va, const vertex & vb, const index_t & delta_v = 0)
	{
		real_t L = *(va - vb);
		size_t N = L / mean_edge;
		L = N;

		while(--N)
		{
			merge_vertices.push_back(vertices.size() + delta_v);
			vertices.push_back(vb + (N / L) * (va - vb));
		}
	};

	auto add_border_vertices = [&](const index_t & i, const index_t & j, const index_t & delta_v = 0) -> index_t
	{
		index_t end_v = j < i ? j + size : j;
		for(index_t v = i; v <= end_v; v++)
		{
			vmap_border[v % size] = vertices[c].size() + delta_v;
			vertices[c].push_back(mesh->gt(border_vertices[v % size]));
			normal += mesh->normal(border_vertices[v % size]);
		}
		return end_v - i + 1;
	};

	for(auto & p: split_indices)
	{
		if(!hole)
		{
			i = p.first;
			j = p.second;

			normal = 0;
			n_v = add_border_vertices(i, j);
			normal /= n_v;

			merge_vertices[c].push_back(n_v - 1);
			gen_vertices(merge_vertices[c], vertices[c], vertices[c].back(), vertices[c].front());
			merge_vertices[c].push_back(0);

			hole = fill_hole_front_angles(vertices[c], mean_edge, normal, max_iter);
		}
		else // fill rest partial holes
		{
			reverse(merge_vertices[!c].begin(), merge_vertices[!c].end());
			for(index_t v: merge_vertices[!c])
				vertices[c].push_back(hole->gt(v));

			normal = 0;
			n_v = 0;

			if(j != p.second)
			{
				n_v = add_border_vertices((j + 1) % size, p.second, hole->n_vertices() - merge_vertices[!c].size());
				merge_vertices[c].push_back(hole->n_vertices() + n_v - 1);
			}
			else merge_vertices[c].push_back(merge_vertices[!c].back());

			gen_vertices(merge_vertices[c], vertices[c], mesh->gt(border_vertices[p.second]), mesh->gt(border_vertices[p.first]), hole->n_vertices() - merge_vertices[!c].size());

			if(i != p.first)
			{
				merge_vertices[c].push_back(vertices[c].size() + hole->n_vertices() - merge_vertices[!c].size());
				n_v += add_border_vertices(p.first, i > 0 ? i - 1 : size - 1 , hole->n_vertices() - merge_vertices[!c].size());
			}
			else merge_vertices[c].push_back(merge_vertices[!c].front());

			normal /= n_v;

			aux_hole = fill_hole_front_angles(vertices[c], mean_edge, normal, max_iter);
			hole->merge(aux_hole, merge_vertices[!c]);

			i = p.first;
			j = p.second;

			delete aux_hole;
		}

		vertices[!c].clear();
		merge_vertices[!c].clear();
		c = !c;
	}

	if(!hole)
	{
		normal = 0;

		for(auto & b: border_vertices)
		{
			normal += mesh->normal(b);
			vertices[c].push_back(mesh->gt(b));
		}
		normal /= vertices[c].size();

		hole = fill_hole_front_angles(vertices[c], mesh->mean_edge(), normal, max_iter);
	}
	else
	{
		reverse(merge_vertices[!c].begin(), merge_vertices[!c].end());

		for(index_t v: merge_vertices[!c])
			vertices[c].push_back(hole->gt(v));

		i = i > 0 ? i - 1 : size - 1;
		j = (j + 1) % size;

		normal = 0;
		n_v = add_border_vertices(j, i, hole->n_vertices() - merge_vertices[!c].size());
		normal /= n_v;

		aux_hole = NULL;
		aux_hole = fill_hole_front_angles(vertices[c], mesh->mean_edge(), normal, max_iter);

		hole->merge(aux_hole, merge_vertices[!c]);
		hole->set_head_vertices(vmap_border, size);

		delete aux_hole;
		delete [] vmap_border;
	}


	if(hole && !hole->is_manifold())
	{
		che_off::write_file(hole, PATH_TEST + "fill_holes/fatal_error.off");
		delete hole;
		return NULL;
	}

	return hole;
}

void split_border(vector<pair<index_t, index_t> > & split_indices, che * mesh, const vector<index_t> & border_vertices)
{
	size_t n = border_vertices.size();
	a_mat data(3, n);
	a_mat means;

	vertex normal;
	for(index_t i = 0; i < n; i++)
	{
		normal = mesh->normal(border_vertices[i]);
		data(0, i) = normal.x;
		data(1, i) = normal.y;
		data(2, i) = normal.z;
		/*
		data(3, i) = mesh->gt(border_vertices[i]).x;
		data(4, i) = mesh->gt(border_vertices[i]).y;
		data(5, i) = mesh->gt(border_vertices[i]).z;
		*/
	}

	//index_t * clusters = NULL;

	index_t k = 2;
	if(kmeans(means, data, k, arma::random_subset, 50, false))
	{
	//	clusters = new index_t[n];
		index_t a, b;
		a = NIL; b = NIL; // review this
		for(index_t i = 0; i < n; i++)
		{
			distance_t d, d_min = INFINITY;
			for(index_t c = 0; c < k; c++)
			{
				d = norm(data.col(i) - means.col(c));
				if(d < d_min)
				{
					d_min = d;
				//	clusters[i] = c;
					b = c;
				}
			}
			if(b != a)
			{
				cerr << b << " " << i << endl;
				a = b;
//				viewer::select_vertices.push_back(border_vertices[i]);
			}
		}
	}
}

vector<index_t> * fill_all_holes(che * mesh, const size_t & max_iter)
{
	vector<index_t> * border_vertices;
	che ** holes;

	tie(border_vertices, holes) = fill_all_holes_meshes(mesh, max_iter);
	if(holes)
	{
		for(index_t b = 0; b < mesh->n_borders(); b++)
			if(holes[b]) delete holes[b];
	}
	delete [] holes;
	return border_vertices;
}

tuple<vector<index_t> *, che **> fill_all_holes_meshes(che * mesh, const size_t & max_iter)
{
	vector<index_t> * border_vertices = NULL;
	che ** holes = NULL;

	const size_t n_borders = mesh->n_borders();
	if(!n_borders) return make_tuple(border_vertices, holes);

	border_vertices = new vector<index_t>[n_borders];
	holes = new che*[n_borders];

	debug_me(inpainting)

	for(index_t b = 0; b < n_borders; b++)
		mesh->border(border_vertices[b], b);

	debug_me(inpainting)
	for(index_t b = 0; b < n_borders; b++)
	{
		debug(b)
//		vector<pair<index_t, index_t> > split_indices;
//		split_border(split_indices, mesh, border_vertices[b]);
//		holes[b] = mesh_fill_hole(mesh, border_vertices[b], max_iter, { {77, 106}, {67, 106}, {38, 11} });
		holes[b] = mesh_fill_hole(mesh, border_vertices[b], max_iter);
	debug_me(inpainting)
		//holes[b]->write_file(PATH_TEST + string("fill_holes/partial") + "_" + to_string(b) + "_" + mesh->name() + ".off");
	debug_me(inpainting)
	}

	debug_me(inpainting)
	for(index_t b = 0; b < n_borders; b++)
		if(holes[b])
		{
	debug_me(inpainting)
			mesh->merge(holes[b], border_vertices[b]);
	debug_me(inpainting)
		}


	debug(mesh->n_borders())
	return make_tuple(border_vertices, holes);
}

che * fill_hole_front_angles_test(che * mesh, vector<index_t> & front_vertices, size_t p_iter, bool & is_grow)
{
	debug_me(filling holes)
	distance_t perimeter = 0.0, init_perimeter = 0.0;

	real_t lenght = mesh->mean_edge();
	priority_queue<border_t> front;

	vector<vertex> vertices;
	vector<index_t> faces;

	for(index_t v: front_vertices)
		vertices.push_back(mesh->gt(v));

	vector<a_vec> tmp_vertices(vertices.size());
	vector<a_vec> tmp_normals(vertices.size());

	vertex normal;
	for(index_t v = 0; v < vertices.size(); v++)
	{
		normal = mesh->normal(front_vertices[v]);
		if(is_grow) normal = -normal;

		tmp_normals[v].resize(3);
		tmp_normals[v](0) = normal.x;
		tmp_normals[v](1) = normal.y;
		tmp_normals[v](2) = normal.z;

		tmp_vertices[v].resize(3);
		tmp_vertices[v](0) = vertices[v].x;
		tmp_vertices[v](1) = vertices[v].y;
		tmp_vertices[v](2) = vertices[v].z;

		if(v) init_perimeter += *(vertices[v] - vertices[v - 1]);
	}


	init_perimeter += *(vertices.back() - vertices.front());
	perimeter = init_perimeter;

//	lenght = perimeter / vertices.size();

	bool o = is_grow;

	vector<bool> is_border(vertices.size());
	vector<array<index_t, 2> > neighbors(vertices.size());

	index_t v, p_v, n_v;
	for(v = 0; v < vertices.size(); v++)
	{
		n_v = (v + 1) % vertices.size();
		p_v = v > 0 ? v - 1: vertices.size() - 1;

		is_border[v] = true;
		neighbors[v][!o] = p_v;
		neighbors[v][o] = n_v;

		front.push(border_t(tmp_vertices, v, neighbors[v], o, tmp_normals[v]));
	//	viewer::vectors.push_back(vertex(tmp_vertices[v](0), tmp_vertices[v](1), tmp_vertices[v](2)));
	//	a_vec normal = tmp_vertices[v] + lenght * 3 * normalise(tmp_normals[v]);
	//	viewer::vectors.push_back(vertex(normal(0), normal(1), normal(2)));
	}

	border_t top;

	angle_t a75 = 75.0 * M_PI / 180;
	angle_t a135 = 135.0 * M_PI / 180;

	a_vec m_vec;
	a_vec m_normal;
	while(!front.empty() && p_iter--)
	{

		while(!front.empty() &&
						!is_border[front.top().v] &&
						neighbors[front.top().v][0] == NIL &&
						neighbors[front.top().v][1] == NIL )
			front.pop();

		if(front.empty()) break;

		top = front.top();
		front.pop();

		v = top.v;
		p_v = neighbors[v][!o];
		n_v = neighbors[v][o];

		if(p_v == n_v) break;
		if(!is_border[p_v] || !is_border[n_v]) break;

		border_t b_p(tmp_vertices, p_v, neighbors[p_v], o, tmp_normals[p_v]);
		border_t b_n(tmp_vertices, n_v, neighbors[n_v], o, tmp_normals[n_v]);

		bool close_vertex = false;
		if(b_p.theta < M_PI && norm(tmp_vertices[n_v] - tmp_vertices[neighbors[p_v][!o]]) < 1.5 * lenght)
			close_vertex = true;
		if(b_n.theta < M_PI && norm(tmp_vertices[p_v] - tmp_vertices[neighbors[n_v][o]]) < 1.5 * lenght)
			close_vertex = true;

		if( top.theta <= M_PI )
		{

			perimeter -= norm(tmp_vertices[v] - tmp_vertices[p_v]);
			perimeter -= norm(tmp_vertices[n_v] - tmp_vertices[v]);
		}

		lenght = ( norm(tmp_vertices[v] - tmp_vertices[p_v]) + norm(tmp_vertices[n_v] - tmp_vertices[v]) ) / 2;
		m_normal = ( tmp_normals[v] + tmp_normals[p_v] + tmp_normals[n_v] ) / 3;
		//m_normal = normalise(m_normal);

		if(top.theta <= a75 || close_vertex)
		{
			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(p_v);

			is_border[v] = false;

			neighbors[v] = {NIL, NIL};
			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = n_v;
			neighbors[n_v][!o] = p_v;
			neighbors[n_v][o] = neighbors[n_v][o];

			front.push(border_t(tmp_vertices, p_v, neighbors[p_v], o, tmp_normals[p_v]));
			front.push(border_t(tmp_vertices, n_v, neighbors[n_v], o, tmp_normals[n_v]));

			perimeter += norm(tmp_vertices[n_v] - tmp_vertices[p_v]);
		}
		else if(top.theta <= a135)
		{
			index_t m_v = tmp_vertices.size();

			m_vec = top.new_vertex(tmp_vertices, 0.5, lenght, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);

			faces.push_back(m_v);
			faces.push_back(v);
			faces.push_back(p_v);

			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(m_v);

//			m_normal = normalise(normal_face(tmp_vertices, m_v, v, p_v) + normal_face(tmp_vertices, n_v, v, m_v));
			tmp_normals.push_back(m_normal);

			is_border[v] = false;
			is_border.push_back(true);

			neighbors[v] = {NIL, NIL};
			neighbors.push_back({NIL, NIL});

			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = m_v;
			neighbors[m_v][!o] = p_v;
			neighbors[m_v][o] = n_v;
			neighbors[n_v][!o] = m_v;
			neighbors[n_v][o] = neighbors[n_v][o];

			front.push(border_t(tmp_vertices, p_v, neighbors[p_v], o, tmp_normals[p_v]));
			front.push(border_t(tmp_vertices, m_v, neighbors[m_v], o, m_normal));
			front.push(border_t(tmp_vertices, n_v, neighbors[n_v], o, tmp_normals[n_v]));

			perimeter += norm(tmp_vertices[m_v] - tmp_vertices[p_v]);
			perimeter += norm(tmp_vertices[n_v] - tmp_vertices[m_v]);
		}
		else if(top.theta <= M_PI)
		{
			index_t m_v = tmp_vertices.size();

			m_vec = top.new_vertex(tmp_vertices, 1./3, lenght, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);
			m_vec = top.new_vertex(tmp_vertices, 2./3, lenght, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);

			faces.push_back(m_v);
			faces.push_back(v);
			faces.push_back(p_v);

			faces.push_back(m_v + 1);
			faces.push_back(v);
			faces.push_back(m_v);

			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(m_v + 1);

//			m_normal = normalise(normal_face(tmp_vertices, m_v, v, p_v) + normal_face(tmp_vertices, m_v + 1, v, m_v));
			tmp_normals.push_back(m_normal);
//			m_normal = normalise(normal_face(tmp_vertices, m_v + 1, v, m_v) + normal_face(tmp_vertices, n_v, v, m_v + 1));
			tmp_normals.push_back(m_normal);

			is_border[v] = false;
			is_border.push_back(true);
			is_border.push_back(true);

			neighbors[v] = {NIL, NIL};
			neighbors.push_back({NIL, NIL});
			neighbors.push_back({NIL, NIL});

			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = m_v;
			neighbors[m_v][!o] = p_v;
			neighbors[m_v][o] = m_v + 1;
			neighbors[m_v + 1][!o] = m_v;
			neighbors[m_v + 1][o] = n_v;
			neighbors[n_v][!o] = m_v + 1;
			neighbors[n_v][o] = neighbors[n_v][o];

			front.push(border_t(tmp_vertices, p_v, neighbors[p_v], o, tmp_normals[p_v]));
			front.push(border_t(tmp_vertices, m_v, neighbors[m_v], o, m_normal));
			front.push(border_t(tmp_vertices, m_v + 1, neighbors[m_v + 1], o, m_normal));
			front.push(border_t(tmp_vertices, n_v, neighbors[n_v], o, tmp_normals[n_v]));

			perimeter += norm(tmp_vertices[p_v] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1 ] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1 ] - tmp_vertices[n_v]);
		}
	}

	if( init_perimeter < perimeter )
	{
		is_grow = true;
	//	return NULL;
	}

	vertices.clear();

	for(a_vec r: tmp_vertices)
		vertices.push_back(vertex(r[0], r[1], r[2]));

	for(index_t v = 0; false && v < tmp_vertices.size(); v++)
	{
	//	viewer::vectors.push_back(vertex(tmp_vertices[v](0), tmp_vertices[v](1), tmp_vertices[v](2)));
		a_vec normal = tmp_vertices[v] + lenght * 3 * normalise(tmp_normals[v]);
	//	viewer::vectors.push_back(vertex(normal(0), normal(1), normal(2)));
	}

	debug(perimeter)
//	debug_me(filling holes)
//	debug(vertices.size())
//	debug(faces.size())
	return faces.size() == 0 ? NULL : new che(vertices.data(), vertices.size(), faces.data(), faces.size() / 3);
}

che * fill_hole_front_angles(vector<vertex> & vertices, const real_t & lenght, const vertex & normal, const size_t & max_iter, bool is_grow)
{
	size_t p_iter = max_iter;
	distance_t perimeter = 0.0;
	distance_t init_perimeter = 0.0;

	priority_queue<border_t> front;
	vector<index_t> faces;

	// PCA --------------------------------------------------------------------------

	a_mat V(3, vertices.size());
	for(index_t v = 0; v < vertices.size(); v++)
	{
		V(0,v) = vertices[v][0];
		V(1,v) = vertices[v][1];
		V(2,v) = vertices[v][2];

		if(v) init_perimeter += norm(V.col(v) - V.col(v-1));
	}

	init_perimeter += *(vertices.back() - vertices.front());
	perimeter = init_perimeter;
	//debug(perimeter)

	a_vec avg = mean(V, 1);
	V.each_col() -= avg;

	a_vec orientation(3);
	orientation(0) = normal.x;
	orientation(1) = normal.y;
	orientation(2) = normal.z;

	a_mat E;
	a_vec eigval;
	eig_sym(eigval, E, V * V.t());
	E.swap_cols(0, 2);
	//debug(E)


	//debug(E * normalise(orientation)) 
	//debug(dot(E.col(2), orientation)) 
	E.col(2) = normalise(dot(orientation, E.col(2)) * E.col(2));
	E.col(1) = normalise(cross(E.col(2), E.col(0)));
	//debug(E)
	//debug(dot(orientation, E.col(2)))

	V = E.t() * V;
//	V.each_col([](a_vec & v){v(2) = 0;});

	bool o = is_grow;

	a_vec a = V.col(0);
	a_vec b = V.col(1);
	a(2) = b(2) = 0;
	a = normalise(a);
	b = normalise(b);

	// END PCA ----------------------------------------------------------------------

	vector<a_vec> tmp_vertices(vertices.size());
	vector<bool> is_border(vertices.size());
	vector<array<index_t, 2> > neighbors(vertices.size());

	index_t v, p_v, n_v;
	for(v = 0; v < vertices.size(); v++)
		tmp_vertices[v] = V.col(v);

	auto push_front = [&front](const border_t & b)
	{
		if(b.theta <= M_PI) front.push(b);
	};

	for(v = 0; v < vertices.size(); v++)
	{
		n_v = (v + 1) % vertices.size();
		p_v = v > 0 ? v - 1: vertices.size() - 1;

		is_border[v] = true;
		neighbors[v][!o] = p_v;
		neighbors[v][o] = n_v;

		border_t aux(tmp_vertices, v, neighbors[v], o);
		if(p_v != NIL && n_v != NIL)
			push_front(aux);
	}

	border_t top;

	angle_t a75 = 75.0 * M_PI / 180;
	angle_t a135 = 135.0 * M_PI / 180;

	a_vec m_vec;
	while(!front.empty() && p_iter-- && p_iter < 2000)
	{
		while(!front.empty() &&
						(!is_border[front.top().v] ||
						neighbors[front.top().v][0] == NIL ||
						neighbors[front.top().v][1] == NIL ) )
			front.pop();

		if(front.empty()) break;

		top = front.top();
		front.pop();

		v = top.v;
		p_v = neighbors[v][!o];
		n_v = neighbors[v][o];

		if(p_v == n_v) break;
		if(!is_border[p_v] || !is_border[n_v]) break;

		border_t b_p(tmp_vertices, p_v, neighbors[p_v], o);
		border_t b_n(tmp_vertices, n_v, neighbors[n_v], o);

		bool close_vertex = false;
		if(b_p.theta <= M_PI && norm(tmp_vertices[n_v] - tmp_vertices[neighbors[p_v][!o]]) < 1.5 * lenght)
			close_vertex = true;
		if(b_n.theta <= M_PI && norm(tmp_vertices[p_v] - tmp_vertices[neighbors[n_v][o]]) < 1.5 * lenght)
			close_vertex = true;

		if( top.theta <= M_PI )
		{
			perimeter -= norm(tmp_vertices[v] - tmp_vertices[p_v]);
			perimeter -= norm(tmp_vertices[n_v] - tmp_vertices[v]);
		}

		if(top.theta <= a75 || close_vertex)
		{
			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(p_v);

			is_border[v] = false;

			//LAST TRIANGLE
			if(neighbors[n_v][o] == p_v)
			{
				perimeter = 0;
				break;
			}

			neighbors[v] = {NIL, NIL};
			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = n_v;
			neighbors[n_v][!o] = p_v;
			neighbors[n_v][o] = neighbors[n_v][o];


			push_front(border_t(tmp_vertices, p_v, neighbors[p_v], o));
			push_front(border_t(tmp_vertices, n_v, neighbors[n_v], o));

			perimeter += norm(tmp_vertices[n_v] - tmp_vertices[p_v]);
		}
		else if(top.theta <= a135)
		{
			index_t m_v = tmp_vertices.size();

			m_vec = top.new_vertex(tmp_vertices, 0.5, lenght, neighbors[v], o);
			tmp_vertices.push_back(m_vec);

			faces.push_back(m_v);
			faces.push_back(v);
			faces.push_back(p_v);

			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(m_v);

			is_border[v] = false;
			is_border.push_back(true);

			neighbors[v] = {NIL, NIL};
			neighbors.push_back({NIL, NIL});

			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = m_v;
			neighbors[m_v][!o] = p_v;
			neighbors[m_v][o] = n_v;
			neighbors[n_v][!o] = m_v;
			neighbors[n_v][o] = neighbors[n_v][o];

			push_front(border_t(tmp_vertices, p_v, neighbors[p_v], o));
			push_front(border_t(tmp_vertices, m_v, neighbors[m_v], o));
			push_front(border_t(tmp_vertices, n_v, neighbors[n_v], o));


			perimeter += norm(tmp_vertices[m_v] - tmp_vertices[p_v]);
			perimeter += norm(tmp_vertices[n_v] - tmp_vertices[m_v]);
		}
		else if(top.theta <= M_PI)
		{
			index_t m_v = tmp_vertices.size();

			m_vec = top.new_vertex(tmp_vertices, 1./3, lenght, neighbors[v], o);
			tmp_vertices.push_back(m_vec);
			m_vec = top.new_vertex(tmp_vertices, 2./3, lenght, neighbors[v], o);
			tmp_vertices.push_back(m_vec);

			faces.push_back(m_v);
			faces.push_back(v);
			faces.push_back(p_v);

			faces.push_back(m_v + 1);
			faces.push_back(v);
			faces.push_back(m_v);

			faces.push_back(n_v);
			faces.push_back(v);
			faces.push_back(m_v + 1);

			is_border[v] = false;
			is_border.push_back(true);
			is_border.push_back(true);

			neighbors[v] = {NIL, NIL};
			neighbors.push_back({NIL, NIL});
			neighbors.push_back({NIL, NIL});

			neighbors[p_v][!o] = neighbors[p_v][!o];
			neighbors[p_v][o] = m_v;
			neighbors[m_v][!o] = p_v;
			neighbors[m_v][o] = m_v + 1;
			neighbors[m_v + 1][!o] = m_v;
			neighbors[m_v + 1][o] = n_v;
			neighbors[n_v][!o] = m_v + 1;
			neighbors[n_v][o] = neighbors[n_v][o];

			push_front(border_t(tmp_vertices, p_v, neighbors[p_v], o));
			push_front(border_t(tmp_vertices, m_v, neighbors[m_v], o));
			push_front(border_t(tmp_vertices, m_v + 1, neighbors[m_v + 1], o));
			push_front(border_t(tmp_vertices, n_v, neighbors[n_v], o));


			perimeter += norm(tmp_vertices[p_v] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1 ] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1 ] - tmp_vertices[n_v]);
		}
	}

	if(init_perimeter < perimeter)
	{
		if(!is_grow)
			return fill_hole_front_angles(vertices, lenght, -normal, max_iter, !is_grow);
		else return NULL;
	}
/*
	a_vec axis;
	axis = avg + E.col(0) * lenght * 3;
	viewer::vectors.push_back(vertex(avg(0), avg(1), avg(2)));
	viewer::vectors.push_back(vertex(axis(0), axis(1), axis(2)));
	axis = avg + E.col(1) * lenght * 3;
	viewer::vectors.push_back(vertex(avg(0), avg(1), avg(2)));
	viewer::vectors.push_back(vertex(axis(0), axis(1), axis(2)));
	axis = avg + E.col(2) * lenght * 3;
	viewer::vectors.push_back(vertex(avg(0), avg(1), avg(2)));
	viewer::vectors.push_back(vertex(axis(0), axis(1), axis(2)));
	axis = avg + normalise(orientation) * lenght * 3;
	viewer::vectors.push_back(vertex(avg(0), avg(1), avg(2)));
	viewer::vectors.push_back(vertex(axis(0), axis(1), axis(2)));
*/
	vertices.clear();
	vertices.reserve(tmp_vertices.size());

	for(a_vec r: tmp_vertices)
	{
		r = E * r + avg;
		vertices.push_back(vertex(r[0], r[1], r[2]));
	}

	return faces.size() ? new che(vertices.data(), vertices.size(), faces.data(), faces.size() / 3) : NULL;
}

void get_real_tri(che * mesh, vector<index_t> & select_vertices, vector<vertex> & triangle, vector<size_t> & tri_sizes )
{
	// Drawing a triangle in the middle of the border
	size_t div = select_vertices.size() / 3;
	size_t r = select_vertices.size() % 3;

	//Defining the ranges
	size_t b = div;
	size_t c = 2 * div;
	size_t d = select_vertices.size();

	vertex tri[3];

	if(r == 2)
		c++;


	for (size_t i = 0; i < d ; i++)
	{
		if(i < b)
		{
			tri[0] += mesh->gt(select_vertices[i]);
			tri_sizes[0]++;
		}
		else if ( i >= b && i < c)
		{
			tri[1] += mesh->gt(select_vertices[i]);
			tri_sizes[1]++;
		}
		else
		{
			tri[2] += mesh->gt(select_vertices[i]);
			tri_sizes[2]++;
		}
	}

	distance_t weight = 1.8;

	distance_t wp = weight / select_vertices.size();
	distance_t aux = wp * tri_sizes[0];
	distance_t wo = (1 - aux) / ( tri_sizes[1] + tri_sizes[2] );

	triangle.push_back( (wp * tri[0]) + wo * (tri[1] + tri[2]) );

	aux = wp * tri_sizes[1];
	wo = (1 - aux) / ( tri_sizes[0] + tri_sizes[2] );

	triangle.push_back( (wp * tri[1]) + wo * (tri[0] + tri[2]) );

	aux = wp * tri_sizes[2];
	wo = (1 - aux) / ( tri_sizes[0] + tri_sizes[1] );

	triangle.push_back( (wp * tri[2]) + wo * (tri[0] + tri[1]) );
}

che * fill_hole_center_triangle(che * mesh, vector<index_t> & select_vertices, index_t index)
{
	size_t n_vertices = select_vertices.size() + 3;
	size_t n_faces = select_vertices.size() + 4;

	vertex * vertices = new vertex[n_vertices];
	index_t * faces = new index_t[n_faces * che::P];

	vector<vertex> triangle;
	vector<size_t> tri_sizes(3,0);

	get_real_tri(mesh, select_vertices, triangle, tri_sizes);

	index_t i = 0;
	for(index_t v: select_vertices)
		vertices[i++] = mesh->gt(v);

	vertices[i++] = triangle[0];
	vertices[i++] = triangle[1];
	vertices[i++] = triangle[2];

	size_t tri_init = select_vertices.size();
	index_t f = 0;

	i = 0;
	for( ; i< tri_sizes[0]-1; i++)
	{
		faces[f++] = i;
		faces[f++] = tri_init;
		faces[f++] = i + 1;
	}

	i++;
	for( ; i < tri_sizes[0] + tri_sizes[1] - 1; i++)
	{
		faces[f++] = i;
		faces[f++] = tri_init + 1;
		faces[f++] = i + 1;
	}

	i++;
	for( ; i < select_vertices.size() - 1; i++)
	{
		faces[f++] = i;
		faces[f++] = tri_init + 2;
		faces[f++] = i + 1;
	}

	size_t aux_i = tri_sizes[0];

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init;
		faces[f++] = tri_init + 1;

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init + 1;
		faces[f++] = aux_i;

	aux_i = tri_sizes[0] + tri_sizes[1];

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init + 1;
		faces[f++] = tri_init + 2;

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init + 2;
		faces[f++] = aux_i;

	aux_i = select_vertices.size();

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init + 2;
		faces[f++] = tri_init;

		faces[f++] = aux_i - 1;
		faces[f++] = tri_init;
		faces[f++] = 0;

	faces[f++] = tri_init + 2;
	faces[f++] = tri_init + 1;
	faces[f++] = tri_init;

	che * new_off = new che(vertices, n_vertices, faces, n_faces);

	delete [] vertices;
	delete [] faces;

	return new_off;
}


} // namespace gproshan

