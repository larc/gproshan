#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class dijkstra
{
	private:
		real_t * weights;
		index_t * predecessors;
		size_t n_vertices;
		index_t source;

	public:
		dijkstra(che * mesh, index_t src);
		~dijkstra();
		real_t & operator()(index_t i);
		index_t & operator[](index_t i);
		void print(std::ostream & os);

	private:
		void run(che * mesh);
};


} // namespace gproshan

#endif // DIJKSTRA_H

