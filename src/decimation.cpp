#include "decimation.h"

decimation::decimation(che * mesh_)
{
	mesh = mesh_;
	Q = new mat[mesh->n_vertices()];

	execute();
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

void decimation::execute()
{
	compute_quadrics();

	index_t * sort_edges = new index_t[mesh->n_edges()];
	vertex_t * error_edges = new vertex_t[mesh->n_edges()];

	order_edges(sort_edges, error_edges);
	corr = mesh->edge_collapse(sort_edges);

	delete [] sort_edges;
	delete [] error_edges;
}

void decimation::compute_quadrics()
{
	vertex n;
	
	#pragma omp parallel for private(n)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		Q[v].resize(4,4);
		Q[v].zeros();
		vec p(4);
		
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
	vec v(4);

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

