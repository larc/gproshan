#include "dijkstra.h"

#include <cstring>
#include <cmath>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


dijkstra::dijkstra(che * shape, index_t src)
{
	n_vertices = shape->n_vertices();
	source = src;

	weights = new distance_t[n_vertices];
	predecessors = new index_t[n_vertices];

	memset(predecessors, 255, sizeof(distance_t)*n_vertices);

	for(index_t i = 0; i < n_vertices; i++)
		weights[i] = INFINITY;

	run(shape);
}

dijkstra::~dijkstra()
{
	if(weights)	delete weights;
	if(predecessors) delete predecessors;
}

distance_t & dijkstra::operator()(index_t i)
{
	return weights[i];
}

index_t & dijkstra::operator[](index_t i)
{
	return predecessors[i];
}

void dijkstra::print(ostream & os)
{
	for(index_t i = 0; i < n_vertices; i++)
		os<<weights[i]<<endl;
}

void dijkstra::run(che * shape)
{
	bool * visited = new bool[n_vertices];
	memset(visited, 0, sizeof(bool)*n_vertices);

	visited[source] = true;
	weights[source] = 0;

	distance_t min;
	index_t min_i;

	for(index_t i = 0; i < n_vertices; i++)
	{
		min = INFINITY;
		min_i = NIL;

		#pragma	omp parallel for num_threads(8)
		for(index_t v = 0; v < n_vertices; v++)
		{
			distance_t w;
			index_t nv;

			if(!visited[v])
			{
				link_t link_he;
				shape->link(link_he, v);

				for(index_t he: link_he)
				{
					nv = shape->vt(next(he));

					if(visited[nv])
					{
						w = weights[nv] + *(shape->get_vertex(nv) - shape->get_vertex(v));

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

	delete visited;
}


} // namespace gproshan

