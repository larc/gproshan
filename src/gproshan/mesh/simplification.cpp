#include <gproshan/mesh/simplification.h>


// geometry processing and shape analysis framework
namespace gproshan {


simplification::simplification(che * mesh_, const index_t levels_)
{
	mesh = mesh_;
	levels = levels_;
	Q = new mat4[mesh->n_vertices];

	execute();
}

simplification::~simplification()
{
	delete [] Q;
}

const mat4 & simplification::operator [] (const index_t i) const
{
	return Q[i];
}

void simplification::execute()
{
	compute_quadrics();
}

void simplification::compute_quadrics()
{
	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
		for(const index_t he: mesh->star(v))
		{
			const vec3 & n = mesh->normal_he(he);
			const vec4 p = (n, -dot(mesh->point(v), n));
			for(int i = 0; i < 4; ++i)
				Q[v][i] += p[i] * p;
		}
}

void simplification::order_edges(index_t * sort_edges, float * error_edges)
{
	#pragma omp parallel for
	for(index_t e = 0; e < mesh->n_edges; ++e)
	{
		sort_edges[e] = e;
		error_edges[e] = compute_error(e);
	}

	std::sort(sort_edges, sort_edges + mesh->n_edges,
		[&error_edges](const index_t a, const index_t b)
		{
			return error_edges[a] < error_edges[b];
		}
		);
}

float simplification::compute_error(const index_t e)
{
	return 0;//as_scalar(v.t() * (Q[mesh->edge_u(e)] + Q[mesh->edge_v(e)]) * v);
}

vertex simplification::create_vertex(const index_t e)
{
	const vertex & va = mesh->vertex_edge_u(e);
	const vertex & vb = mesh->vertex_edge_v(e);

	return (va + vb) / 2;
}


} // namespace gproshan

