#include "decimation.h"

decimation::decimation(che * mesh_)
{
	mesh = mesh_;
	Q = new mat[mesh->n_vertices()];
}

decimation::~decimation()
{
	delete [] Q;
}

void decimation::compute_quadrics()
{
	vec p(4);
	vertex n;

	#pragma omp parallel for private(p, n)
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		Q[v].resize(4,4);
		Q[v].zeros();
		
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

priority_queue< err_edge > decimation::order_edges()
{
	priority_queue< err_edge > queue;	
	index_t e;
	err_edge p_edge;
	for(int i = 0; i < mesh->n_edges(); i++)
	{	
		e = mesh->vt_e(i);
		p_edge.edge = e;
		p_edge.error = compute_error(e);
		queue.push( p_edge ) ;
	}
	return queue;
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

