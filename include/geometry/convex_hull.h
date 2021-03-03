#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include "mesh/vertex.h" 

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


///< 2D Convex Hull 
class convex_hull
{
	private:
		std::vector<vertex> CH;		///< convex hull points clockwise

	public:
		convex_hull(std::vector<vertex> & points);
		convex_hull(vertex * points, const size_t & n_points);
		operator const std::vector<vertex> & ();

	private:
		void andrew_algorithm(vertex * points, const size_t & n_points);
		bool ccw(const vertex & p, const vertex & q, const vertex & r);
};


} // namespace gproshan

#endif // CONVEX_HULL_H

