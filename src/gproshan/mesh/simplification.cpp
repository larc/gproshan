#include <gproshan/mesh/simplification.h>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


simplification::simplification(che * mesh_, const index_t & levels_)
{
	mesh = mesh_;
	levels = levels_;
	Q = new a_mat[mesh->n_vertices];

	execute();
}

simplification::~simplification()
{
	if(Q) delete [] Q;
}

void simplification::execute()
{
	compute_quadrics();
}

void simplification::compute_quadrics()
{
	vertex n;

	#pragma omp parallel for private(n)
	for(index_t v = 0; v < mesh->n_vertices; ++v)
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

void simplification::order_edges(index_t * const & sort_edges, real_t * const & error_edges)
{
	#pragma omp parallel for
	for(index_t e = 0; e < mesh->n_edges; ++e)
	{
		sort_edges[e] = e;
		error_edges[e] = compute_error(e);
	}

	sort(sort_edges, sort_edges + mesh->n_edges,
		[&error_edges](const index_t & a, const index_t & b)
		{
			return error_edges[a] < error_edges[b];
		}
		);
}

real_t simplification::compute_error(const index_t & e)
{
	vertex ve = create_vertex(e);
	a_vec v(4);

	v(0) = ve.x;
	v(1) = ve.y;
	v(2) = ve.z;
	v(3) = 1;

	return as_scalar(v.t() * (Q[mesh->vt_e(e)] + Q[mesh->vt_e(e, true)]) * v);
}

vertex simplification::create_vertex(const index_t & e)
{
	vertex va, vb;
	va = mesh->gt_e(e);
	vb = mesh->gt_e(e, true);

	return (va + vb) / 2;
}


} // namespace gproshan
