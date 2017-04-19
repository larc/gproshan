#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "che.h"
#include "include.h"

class dijkstra
{
	private:
		distance_t * weights;
		index_t * predecessors;
		size_t n_vertices;
		index_t source;

	public:
		dijkstra(che * shape, index_t src);
		~dijkstra();
		distance_t & operator()(index_t i);
		index_t & operator[](index_t i);
		void print(ostream & os);

	private:
		void run(che * shape);
};

#endif // DIJKSTRA_H
