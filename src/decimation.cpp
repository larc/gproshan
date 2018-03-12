#include "decimation.h"

decimation::decimation(che * mesh_, const vertex *const & normals, const index_t & levels_)
{
	mesh = mesh_;
	levels = levels_;
	Q = new a_mat[mesh->n_vertices()];

	execute(normals);
}

decimation::~decimation()
{
	if(Q) delete [] Q;
	if(corr) delete [] corr;
}

decimation::operator const corr_t * ()
{
	return corr;
}

void decimation::execute(const vertex *const & normals)
{
	compute_quadrics();

	const size_t n_vertices = mesh->n_vertices();
	index_t * sort_edges = new index_t[mesh->n_edges()];
	vertex_t * error_edges = new vertex_t[mesh->n_edges()];
	vertex * corr_v = new vertex[n_vertices];
	index_t * corr_i = new index_t[n_vertices * P];

	corr_t * corr_aux;
	index_t he, vi;
	vertex a, b, c;

	auto add_he_trigs = [this](vector<index_t> & he_trigs, const corr_t & c)
	{
		const index_t he = c.t * P;
		for_star(t_he, mesh, mesh->vt(he))
			he_trigs.push_back(t_he);
		for_star(t_he, mesh, mesh->vt(next(he)))
			he_trigs.push_back(t_he);
		for_star(t_he, mesh, mesh->vt(prev(he)))
			he_trigs.push_back(t_he);
	};

	order_edges(sort_edges, error_edges);
	corr = mesh->edge_collapse(sort_edges, normals);

	while(--levels)
	{
		#pragma omp parallel for private(he, vi)
		for(index_t v = 0; v < n_vertices; v++)
		{
			corr_v[v] = mesh->corr_vertex(corr[v]);
			he = corr[v].t * P;
			vi = v * P;
			corr_i[vi] = mesh->vt(he);
			corr_i[vi + 1] = mesh->vt(next(he));
			corr_i[vi + 2] = mesh->vt(prev(he));
		}

		order_edges(sort_edges, error_edges);
		corr_aux = mesh->edge_collapse(sort_edges, normals);

		#pragma omp parallel for private(vi, a, b, c)
		for(index_t v = 0; v < n_vertices; v++)
		{
			vi = v * P;
			a = mesh->corr_vertex(corr_aux[corr_i[vi]]);
			b = mesh->corr_vertex(corr_aux[corr_i[vi + 1]]);
			c = mesh->corr_vertex(corr_aux[corr_i[vi + 2]]);
			corr_v[v] = corr[v].alpha[0] * a + corr[v].alpha[1] * b + corr[v].alpha[2] * c;

			vector<index_t> he_trigs;
			add_he_trigs(he_trigs, corr_aux[corr_i[vi]]);
			add_he_trigs(he_trigs, corr_aux[corr_i[vi + 1]]);
			add_he_trigs(he_trigs, corr_aux[corr_i[vi + 2]]);

			corr[v] = mesh->find_corr(corr_v[v], normals[v], he_trigs);
		}

		delete [] corr_aux;
	}

	delete [] sort_edges;
	delete [] error_edges;
	delete [] corr_v;
	delete [] corr_i;
}

void decimation::compute_quadrics()
{
	vertex n;

	#pragma omp parallel for private(n)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		Q[v].resize(4,4);
		Q[v].zeros();
		a_vec p(4);

		for_star(he, mesh, v)
		{
			n = mesh->normal_he(he);
			p(0) = n.x;
			p(1) = n.y;
			p(2) = n.z;
			p(3) = -(n, mesh->gt(v));

			Q[v] += p * p.t();
		}
	}
}

void decimation::order_edges(index_t * const & sort_edges, vertex_t * const & error_edges)
{
	#pragma omp parallel for
	for(int e = 0; e < mesh->n_edges(); e++)
	{
		sort_edges[e] = e;
		error_edges[e] = compute_error(e);
	}

	sort(sort_edges, sort_edges + mesh->n_edges(),
		[&error_edges](const index_t & a, const index_t & b)
		{
			return error_edges[a] < error_edges[b];
		}
		);
}

vertex_t decimation::compute_error(const index_t & e)
{
	vertex ve = create_vertex(e);
	a_vec v(4);

	v(0) = ve.x;
	v(1) = ve.y;
	v(2) = ve.z;
	v(3) = 1;

	return as_scalar(v.t() * (Q[mesh->vt_e(e)] + Q[mesh->vt_e(e, true)]) * v);
}

vertex decimation::create_vertex(const index_t & e)
{
	vertex va, vb;
	va = mesh->gt_e(e);
	vb = mesh->gt_e(e, true);

	return (va + vb) / 2;
}

