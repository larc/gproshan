#include <gproshan/geodesics/dijkstra.h>

#include <cstring>
#include <cmath>


// geometry processing and shape analysis framework
namespace gproshan {


dijkstra::dijkstra(che * mesh, index_t src)
{
	n_vertices = mesh->n_vertices;
	source = src;

	weights = new real_t[n_vertices];
	predecessors = new index_t[n_vertices];

	memset(predecessors, 255, sizeof(real_t)*n_vertices);

	for(index_t i = 0; i < n_vertices; ++i)
		weights[i] = INFINITY;

	run(mesh);
}

dijkstra::~dijkstra()
{
	if(weights)	delete weights;
	if(predecessors) delete predecessors;
}

real_t & dijkstra::operator()(index_t i)
{
	return weights[i];
}

index_t & dijkstra::operator[](index_t i)
{
	return predecessors[i];
}

void dijkstra::print(std::ostream & os)
{
	for(index_t i = 0; i < n_vertices; ++i)
		os << weights[i] << std::endl;
}

void dijkstra::run(che * mesh)
{
	bool * visited = new bool[n_vertices];
	memset(visited, 0, sizeof(bool) * n_vertices);

	visited[source] = true;
	weights[source] = 0;

	real_t min;
	index_t min_i;

	for(index_t i = 0; i < n_vertices; ++i)
	{
		min = INFINITY;
		min_i = NIL;

		#pragma	omp parallel for num_threads(8)
		for(index_t v = 0; v < n_vertices; ++v)
		{
			real_t w;

			if(!visited[v])
			{
				for(const index_t nv: mesh->link(v))
				{
					if(visited[nv])
					{
						w = weights[nv] + norm(mesh->point(nv) - mesh->point(v));

						if(w < weights[v])
							weights[v] = w;
					}
				}

				#pragma omp critical
				if(weights[v] < min)
				{
					min = weights[v];
					min_i = v;
				}
			}
		}

		if(min_i != NIL) visited[min_i] = true;
	}

	delete [] visited;
}


} // namespace gproshan

