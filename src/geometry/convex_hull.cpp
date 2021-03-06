#include "geometry/convex_hull.h"


#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan {


convex_hull::convex_hull(std::vector<vertex> & points): convex_hull(points.data(), points.size()) {}

convex_hull::convex_hull(vertex * points, const size_t & n_points)
{
	andrew_algorithm(points, n_points);
}

convex_hull::operator const std::vector<vertex> & ()
{
	return CH;
}

///< Andrew's Convex Hull Algorithm: Competitive Programming 4
void convex_hull::andrew_algorithm(vertex * points, const size_t & n_points)
{
	std::sort(points, points + n_points);

	CH.resize(2 * n_points);

	index_t k = 0;
	for(index_t i = 0; i < n_points; ++i)
	{
		while(k > 1 && !ccw(CH[k - 2], CH[k - 1], points[i])) --k;
		CH[k++] = points[i];
	}
	
	index_t t = k;
	for(index_t i = n_points - 2; i > 0; --i)
	{
		while(k > t && !ccw(CH[k - 2], CH[k - 1], points[i])) --k;
		CH[k++] = points[i];
	}
	
	while(k > t && !ccw(CH[k - 2], CH[k - 1], points[0])) --k;

	CH.resize(k);
}

bool convex_hull::ccw(const vertex & p, const vertex & q, const vertex & r)
{
	return ((q - p) * (r - p)).z > 0;
}


} // namespace gproshan

