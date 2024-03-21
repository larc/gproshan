#include <gproshan/mesh/che_fill_hole.h>

#include <gproshan/mesh/che_off.h>
#include <gproshan/laplacian/laplacian.h>

#include <queue>


// geometry processing and shape analysis framework
namespace gproshan {


bool operator<(const border_t & a, const border_t & b)
{
	return a.theta > b.theta;
}

a_vec normal_face(const std::vector<a_vec> & tmp_vertices, const index_t a_v, const index_t b_v, const index_t c_v)
{
	a_vec a = tmp_vertices[c_v] - tmp_vertices[a_v];
	a_vec b = tmp_vertices[b_v] - tmp_vertices[a_v];
	return normalise(cross(a,b));
}

che * mesh_simple_fill_hole(che * mesh, const std::vector<index_t> & border_vertices, const size_t max_iter = 1000)
{
	std::vector<vertex> vertices;
	vertex normal, normal_v, edge_v, v;
	a_mat E(3, 3);
	a_vec ve(3);

	vertices.reserve(size(border_vertices));

	for(const index_t b: border_vertices)
	{
		v = mesh->point(b);
		normal_v = mesh->normal(b);
		edge_v = mesh->vertex_he(he_next(mesh->evt(b))) - v;
		edge_v -= dot(normal_v, edge_v) * normal_v;

		E(0, 2) = normal_v.x();
		E(1, 2) = normal_v.y();
		E(2, 2) = normal_v.z();

		E(0, 0) = edge_v.x();
		E(1, 0) = edge_v.y();
		E(2, 0) = edge_v.z();

		E.col(1) = cross(E.col(2), E.col(0));
		E = normalise(E);

		ve(0) = v.x(); ve(1) = v.y(); ve(2) = v.z();

		ve = E.t() * ve;
		vertices.push_back(*((vertex *) ve.memptr()));
	}

	return fill_hole_front_angles(vertices, mesh->mean_edge(), normal, max_iter);
}

che * mesh_fill_hole(che * mesh, const std::vector<index_t> & border_vertices, const size_t max_iter, const std::vector<std::pair<index_t, index_t> > & split_indices = {})
{
	std::vector<vertex> vertices[2];
	std::vector<index_t> merge_vertices[2];
	vertex normal;

	size_t nb = size(border_vertices);
	index_t i = NIL, j, n_v;
	che * hole = nullptr;
	che * aux_hole;

	index_t * vmap_border = new index_t[nb];

	index_t c = 1;
	real_t mean_edge = mesh->mean_edge();

	auto gen_vertices = [&mean_edge](std::vector<index_t> & merge_vertices, std::vector<vertex> & vertices, const vertex & va, const vertex & vb, const index_t delta_v = 0)
	{
		real_t L = length(va - vb);
		size_t N = L / mean_edge;
		L = N;

		while(--N)
		{
			merge_vertices.push_back(size(vertices) + delta_v);
			vertices.push_back(vb + (N / L) * (va - vb));
		}
	};

	auto add_border_vertices = [&](const index_t i, const index_t j, const index_t delta_v = 0) -> index_t
	{
		index_t end_v = j < i ? j + nb : j;
		for(index_t v = i; v <= end_v; ++v)
		{
			vmap_border[v % nb] = size(vertices[c]) + delta_v;
			vertices[c].push_back(mesh->point(border_vertices[v % nb]));
			normal += mesh->normal(border_vertices[v % nb]);
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
			reverse(begin(merge_vertices[!c]), end(merge_vertices[!c]));
			for(index_t v: merge_vertices[!c])
				vertices[c].push_back(hole->point(v));

			normal = 0;
			n_v = 0;

			if(j != p.second)
			{
				n_v = add_border_vertices((j + 1) % nb, p.second, hole->n_vertices - size(merge_vertices[!c]));
				merge_vertices[c].push_back(hole->n_vertices + n_v - 1);
			}
			else merge_vertices[c].push_back(merge_vertices[!c].back());

			gen_vertices(merge_vertices[c], vertices[c], mesh->point(border_vertices[p.second]), mesh->point(border_vertices[p.first]), hole->n_vertices - size(merge_vertices[!c]));

			if(i != p.first)
			{
				merge_vertices[c].push_back(std::size(vertices[c]) + hole->n_vertices - size(merge_vertices[!c]));
				n_v += add_border_vertices(p.first, i > 0 ? i - 1 : nb - 1, hole->n_vertices - size(merge_vertices[!c]));
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
			vertices[c].push_back(mesh->point(b));
		}
		normal /= size(vertices[c]);

		hole = fill_hole_front_angles(vertices[c], mesh->mean_edge(), normal, max_iter);
	}
	else
	{
		reverse(begin(merge_vertices[!c]), end(merge_vertices[!c]));

		for(index_t v: merge_vertices[!c])
			vertices[c].push_back(hole->point(v));

		i = i > 0 ? i - 1 : nb - 1;
		j = (j + 1) % nb;

		normal = 0;
		n_v = add_border_vertices(j, i, hole->n_vertices - size(merge_vertices[!c]));
		normal /= n_v;

		aux_hole = nullptr;
		aux_hole = fill_hole_front_angles(vertices[c], mesh->mean_edge(), normal, max_iter);

		hole->merge(aux_hole, merge_vertices[!c]);
		hole->set_head_vertices(vmap_border, nb);

		delete aux_hole;
		delete [] vmap_border;
	}


	if(hole && !hole->is_manifold())
	{
		che_off::write_file(hole, tmp_file_path("fill_holes_error.off"));
		delete hole;
		return nullptr;
	}

	return hole;
}

void split_border(std::vector<std::pair<index_t, index_t> > & , che * mesh, const std::vector<index_t> & border_vertices)
{
	size_t n = size(border_vertices);
	a_mat data(3, n);
	a_mat means;

	vertex normal;
	for(index_t i = 0; i < n; ++i)
	{
		normal = mesh->normal(border_vertices[i]);
		data(0, i) = normal.x();
		data(1, i) = normal.y();
		data(2, i) = normal.z();
		/*
		data(3, i) = mesh->point(border_vertices[i]).x();
		data(4, i) = mesh->point(border_vertices[i]).y();
		data(5, i) = mesh->point(border_vertices[i]).z();
		*/
	}

	//index_t * clusters = nullptr;

	index_t k = 2;
	if(kmeans(means, data, k, arma::random_subset, 50, false))
	{
	//	clusters = new index_t[n];
		index_t a, b;
		a = NIL; b = NIL; // review this
		for(index_t i = 0; i < n; ++i)
		{
			real_t d, d_min = INFINITY;
			for(index_t c = 0; c < k; ++c)
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
				std::cerr << b << " " << i << std::endl;
				a = b;
			}
		}
	}
}

std::vector<index_t> * fill_all_holes(che * mesh, const size_t max_iter)
{
	gproshan_error(holes);
	std::vector<index_t> * border_vertices;
	che ** holes;
	gproshan_error(holes);

	tie(border_vertices, holes) = fill_all_holes_meshes(mesh, max_iter);
	gproshan_error(holes);
	if(holes)
	{
		// FIX_BOUND
		/*
		for(index_t b = 0; b < mesh->n_borders; ++b)
			if(holes[b]) delete holes[b];
		*/
	}
	//delete [] holes;
	return border_vertices;
}

std::tuple<std::vector<index_t> *, che **> fill_all_holes_meshes(che * mesh, const size_t max_iter)
{
	std::vector<index_t> * border_vertices = nullptr;
	che ** holes = nullptr;

	std::vector<index_t> bounds = mesh->bounds();
	const size_t n_borders = size(bounds);

	if(!n_borders) return make_tuple(border_vertices, holes);

	border_vertices = new std::vector<index_t>[n_borders];
	holes = new che*[n_borders];

	gproshan_debug(inpainting);
	gproshan_error(holes);

	for(index_t b = 0; b < n_borders; ++b)
		border_vertices[b] = mesh->boundary(bounds[b]);

	gproshan_debug(inpainting);
	gproshan_error(holes);
	for(index_t b = 0; b < n_borders; ++b)
	{
	gproshan_error(holes);
		gproshan_debug_var(b);
//		std::vector<std::pair<index_t, index_t> > split_indices;
//		split_border(split_indices, mesh, border_vertices[b]);
//		holes[b] = mesh_fill_hole(mesh, border_vertices[b], max_iter, { {77, 106}, {67, 106}, {38, 11} });
		holes[b] = mesh_fill_hole(mesh, border_vertices[b], max_iter);
	gproshan_debug(inpainting);
		if(holes[b]) che_off::write_file(holes[b], tmp_file_path("fill_holes_" + std::to_string(b) + "_" + mesh->name() + ".off"));
	gproshan_debug(inpainting);
	gproshan_error(holes);
	}

	che * old = nullptr;
	for(index_t b = 0; b < n_borders; ++b)
		if(holes[b])
		{
			old = mesh;
			mesh = old->merge(holes[b], border_vertices[b]);
			delete old;
		}


	return make_tuple(border_vertices, holes);
}

che * fill_hole_front_angles_test(che * mesh, std::vector<index_t> & front_vertices, size_t p_iter, bool & is_grow)
{
	gproshan_debug(filling holes);
	real_t perimeter = 0.0, init_perimeter = 0.0;

	real_t length = mesh->mean_edge();
	std::priority_queue<border_t> front;

	std::vector<vertex> vertices;
	std::vector<index_t> trigs;

	for(index_t v: front_vertices)
		vertices.push_back(mesh->point(v));

	std::vector<a_vec> tmp_vertices(size(vertices));
	std::vector<a_vec> tmp_normals(size(vertices));

	vertex normal;
	for(index_t v = 0; v < size(vertices); ++v)
	{
		normal = mesh->normal(front_vertices[v]);
		if(is_grow) normal = -normal;

		tmp_normals[v].resize(3);
		tmp_normals[v](0) = normal.x();
		tmp_normals[v](1) = normal.y();
		tmp_normals[v](2) = normal.z();

		tmp_vertices[v].resize(3);
		tmp_vertices[v](0) = vertices[v].x();
		tmp_vertices[v](1) = vertices[v].y();
		tmp_vertices[v](2) = vertices[v].z();

		if(v) init_perimeter += norm(vertices[v] - vertices[v - 1]);
	}


	init_perimeter += norm(vertices.back() - vertices.front());
	perimeter = init_perimeter;

//	length = perimeter / size(vertices);

	bool o = is_grow;

	std::vector<bool> is_border(size(vertices));
	std::vector<std::array<index_t, 2> > neighbors(size(vertices));

	index_t v, p_v, n_v;
	for(v = 0; v < size(vertices); ++v)
	{
		n_v = (v + 1) % size(vertices);
		p_v = v > 0 ? v - 1: size(vertices) - 1;

		is_border[v] = true;
		neighbors[v][!o] = p_v;
		neighbors[v][o] = n_v;

		front.push(border_t(tmp_vertices, v, neighbors[v], o, tmp_normals[v]));
	}

	border_t top;

	real_t a75 = 75.0 * M_PI / 180;
	real_t a135 = 135.0 * M_PI / 180;

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
		if(b_p.theta < M_PI && norm(tmp_vertices[n_v] - tmp_vertices[neighbors[p_v][!o]]) < 1.5 * length)
			close_vertex = true;
		if(b_n.theta < M_PI && norm(tmp_vertices[p_v] - tmp_vertices[neighbors[n_v][o]]) < 1.5 * length)
			close_vertex = true;

		if( top.theta <= M_PI )
		{

			perimeter -= norm(tmp_vertices[v] - tmp_vertices[p_v]);
			perimeter -= norm(tmp_vertices[n_v] - tmp_vertices[v]);
		}

		length = ( norm(tmp_vertices[v] - tmp_vertices[p_v]) + norm(tmp_vertices[n_v] - tmp_vertices[v]) ) / 2;
		m_normal = ( tmp_normals[v] + tmp_normals[p_v] + tmp_normals[n_v] ) / 3;
		//m_normal = normalise(m_normal);

		if(top.theta <= a75 || close_vertex)
		{
			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

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
			index_t m_v = size(tmp_vertices);

			m_vec = top.new_vertex(tmp_vertices, 0.5, length, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);

			trigs.push_back(m_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(m_v);

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
			index_t m_v = size(tmp_vertices);

			m_vec = top.new_vertex(tmp_vertices, 1./3, length, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);
			m_vec = top.new_vertex(tmp_vertices, 2./3, length, neighbors[v], o, tmp_normals[v]);
			tmp_vertices.push_back(m_vec);

			trigs.push_back(m_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

			trigs.push_back(m_v + 1);
			trigs.push_back(v);
			trigs.push_back(m_v);

			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(m_v + 1);

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
			perimeter += norm(tmp_vertices[m_v + 1] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1] - tmp_vertices[n_v]);
		}
	}

	if( init_perimeter < perimeter )
	{
		is_grow = true;
	//	return nullptr;
	}

	vertices.clear();

	for(a_vec r: tmp_vertices)
		vertices.push_back({r[0], r[1], r[2]});

	for(index_t v = 0; false && v < size(tmp_vertices); ++v)
		a_vec normal = tmp_vertices[v] + length * 3 * normalise(tmp_normals[v]);

	gproshan_debug_var(perimeter);
//	gproshan_debug(filling holes);
//	gproshan_debug_var(size(vertices));
//	gproshan_debug_var(size(trigs));
	return size(trigs) == 0 ? nullptr : new che(vertices.data(), size(vertices), trigs.data(), size(trigs) / 3);
}

che * fill_hole_front_angles(std::vector<vertex> & vertices, const real_t length, const vertex & normal, const size_t max_iter, bool is_grow)
{
	size_t p_iter = max_iter;
	real_t perimeter = 0.0;
	real_t init_perimeter = 0.0;

	std::priority_queue<border_t> front;
	std::vector<index_t> trigs;

	// PCA --------------------------------------------------------------------------

	a_mat V(3, size(vertices));
	for(index_t v = 0; v < size(vertices); ++v)
	{
		V(0,v) = vertices[v][0];
		V(1,v) = vertices[v][1];
		V(2,v) = vertices[v][2];

		if(v) init_perimeter += norm(V.col(v) - V.col(v-1));
	}

	init_perimeter += norm(vertices.back() - vertices.front());
	perimeter = init_perimeter;
	//debug(perimeter)

	a_vec avg = mean(V, 1);
	V.each_col() -= avg;

	a_vec orientation(3);
	orientation(0) = normal.x();
	orientation(1) = normal.y();
	orientation(2) = normal.z();

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

	std::vector<a_vec> tmp_vertices(size(vertices));
	std::vector<bool> is_border(size(vertices));
	std::vector<std::array<index_t, 2> > neighbors(size(vertices));

	index_t v, p_v, n_v;
	for(v = 0; v < size(vertices); ++v)
		tmp_vertices[v] = V.col(v);

	auto push_front = [&front](const border_t & b)
	{
		if(b.theta <= M_PI) front.push(b);
	};

	for(v = 0; v < size(vertices); ++v)
	{
		n_v = (v + 1) % size(vertices);
		p_v = v > 0 ? v - 1: size(vertices) - 1;

		is_border[v] = true;
		neighbors[v][!o] = p_v;
		neighbors[v][o] = n_v;

		border_t aux(tmp_vertices, v, neighbors[v], o);
		if(p_v != NIL && n_v != NIL)
			push_front(aux);
	}

	border_t top;

	real_t a75 = 75.0 * M_PI / 180;
	real_t a135 = 135.0 * M_PI / 180;

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
		if(b_p.theta <= M_PI && norm(tmp_vertices[n_v] - tmp_vertices[neighbors[p_v][!o]]) < 1.5 * length)
			close_vertex = true;
		if(b_n.theta <= M_PI && norm(tmp_vertices[p_v] - tmp_vertices[neighbors[n_v][o]]) < 1.5 * length)
			close_vertex = true;

		if( top.theta <= M_PI )
		{
			perimeter -= norm(tmp_vertices[v] - tmp_vertices[p_v]);
			perimeter -= norm(tmp_vertices[n_v] - tmp_vertices[v]);
		}

		if(top.theta <= a75 || close_vertex)
		{
			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

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
			index_t m_v = size(tmp_vertices);

			m_vec = top.new_vertex(tmp_vertices, 0.5, length, neighbors[v], o);
			tmp_vertices.push_back(m_vec);

			trigs.push_back(m_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(m_v);

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
			index_t m_v = size(tmp_vertices);

			m_vec = top.new_vertex(tmp_vertices, 1./3, length, neighbors[v], o);
			tmp_vertices.push_back(m_vec);
			m_vec = top.new_vertex(tmp_vertices, 2./3, length, neighbors[v], o);
			tmp_vertices.push_back(m_vec);

			trigs.push_back(m_v);
			trigs.push_back(v);
			trigs.push_back(p_v);

			trigs.push_back(m_v + 1);
			trigs.push_back(v);
			trigs.push_back(m_v);

			trigs.push_back(n_v);
			trigs.push_back(v);
			trigs.push_back(m_v + 1);

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
			perimeter += norm(tmp_vertices[m_v + 1] - tmp_vertices[m_v]);
			perimeter += norm(tmp_vertices[m_v + 1] - tmp_vertices[n_v]);
		}
	}

	if(init_perimeter < perimeter)
	{
		if(!is_grow)
			return fill_hole_front_angles(vertices, length, -normal, max_iter, !is_grow);
		else return nullptr;
	}

	vertices.clear();
	vertices.reserve(size(tmp_vertices));

	for(a_vec r: tmp_vertices)
	{
		r = E * r + avg;
		vertices.push_back({r[0], r[1], r[2]});
	}

	return size(trigs) ? new che(vertices.data(), size(vertices), trigs.data(), size(trigs) / 3) : nullptr;
}

void get_real_tri(che * mesh, std::vector<index_t> & select_vertices, std::vector<vertex> & triangle, std::vector<size_t> & tri_sizes )
{
	// Drawing a triangle in the middle of the border
	size_t div = size(select_vertices) / 3;
	size_t r = size(select_vertices) % 3;

	//Defining the ranges
	size_t b = div;
	size_t c = 2 * div;
	size_t d = size(select_vertices);

	vertex tri[3];

	if(r == 2) ++c;


	for (size_t i = 0; i < d ; ++i)
	{
		if(i < b)
		{
			tri[0] += mesh->point(select_vertices[i]);
			++tri_sizes[0];
		}
		else if ( i >= b && i < c)
		{
			tri[1] += mesh->point(select_vertices[i]);
			++tri_sizes[1];
		}
		else
		{
			tri[2] += mesh->point(select_vertices[i]);
			++tri_sizes[2];
		}
	}

	real_t weight = 1.8;

	real_t wp = weight / size(select_vertices);
	real_t aux = wp * tri_sizes[0];
	real_t wo = (1 - aux) / ( tri_sizes[1] + tri_sizes[2] );

	triangle.push_back( (wp * tri[0]) + wo * (tri[1] + tri[2]) );

	aux = wp * tri_sizes[1];
	wo = (1 - aux) / ( tri_sizes[0] + tri_sizes[2] );

	triangle.push_back( (wp * tri[1]) + wo * (tri[0] + tri[2]) );

	aux = wp * tri_sizes[2];
	wo = (1 - aux) / ( tri_sizes[0] + tri_sizes[1] );

	triangle.push_back( (wp * tri[2]) + wo * (tri[0] + tri[1]) );
}

che * fill_hole_center_triangle(che * mesh, std::vector<index_t> & select_vertices, index_t )
{
	size_t n_vertices = size(select_vertices) + 3;
	size_t n_trigs = size(select_vertices) + 4;

	vertex * vertices = new vertex[n_vertices];
	index_t * trigs = new index_t[n_trigs * 3];

	std::vector<vertex> triangle;
	std::vector<size_t> tri_sizes(3,0);

	get_real_tri(mesh, select_vertices, triangle, tri_sizes);

	index_t i = 0;
	for(index_t v: select_vertices)
		vertices[i++] = mesh->point(v);

	vertices[i++] = triangle[0];
	vertices[i++] = triangle[1];
	vertices[i++] = triangle[2];

	size_t tri_init = size(select_vertices);
	index_t f = 0;

	i = 0;
	for( ; i< tri_sizes[0]-1; ++i)
	{
		trigs[f++] = i;
		trigs[f++] = tri_init;
		trigs[f++] = i + 1;
	}

	++i;
	for( ; i < tri_sizes[0] + tri_sizes[1] - 1; ++i)
	{
		trigs[f++] = i;
		trigs[f++] = tri_init + 1;
		trigs[f++] = i + 1;
	}

	++i;
	for( ; i < size(select_vertices) - 1; ++i)
	{
		trigs[f++] = i;
		trigs[f++] = tri_init + 2;
		trigs[f++] = i + 1;
	}

	size_t aux_i = tri_sizes[0];

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init;
		trigs[f++] = tri_init + 1;

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init + 1;
		trigs[f++] = aux_i;

	aux_i = tri_sizes[0] + tri_sizes[1];

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init + 1;
		trigs[f++] = tri_init + 2;

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init + 2;
		trigs[f++] = aux_i;

	aux_i = size(select_vertices);

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init + 2;
		trigs[f++] = tri_init;

		trigs[f++] = aux_i - 1;
		trigs[f++] = tri_init;
		trigs[f++] = 0;

	trigs[f++] = tri_init + 2;
	trigs[f++] = tri_init + 1;
	trigs[f++] = tri_init;

	che * new_off = new che(vertices, n_vertices, trigs, n_trigs);

	delete [] vertices;
	delete [] trigs;

	return new_off;
}


} // namespace gproshan

