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
		convex_hull(const std::vector<vertex> & points, const real_t precision = 1000);
		convex_hull(const vertex * points, const size_t & n_points, const real_t precision = 1000);
		operator const std::vector<index_t> & () const;

	private:
		void andrew_algorithm(const std::vector<ivec2> & points);
		bool ccw(const ivec2 & p, const ivec2 & q, const ivec2 & r);
};


} // namespace gproshan

#endif // CONVEX_HULL_H

