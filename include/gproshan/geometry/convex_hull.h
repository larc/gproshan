#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <gproshan/geometry/vec.h>

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


using vertex = vec3;


///< 2D Convex Hull
class convex_hull
{
	private:
		std::vector<index_t> CH;		///< convex hull points clockwise

	public:
		convex_hull(const std::vector<vertex> & points);
		convex_hull(const vertex * points, const size_t & n_points);
		operator const std::vector<index_t> & ();

	private:
		void andrew_algorithm(const vertex * points, const size_t & n_points);
		bool ccw(const vertex & p, const vertex & q, const vertex & r);
};


} // namespace gproshan

#endif // CONVEX_HULL_H

