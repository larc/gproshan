#include <gproshan/geometry/convex_hull.h>

#include <algorithm>
#include <numeric>


// geometry processing and shape analysis framework
namespace gproshan {


convex_hull::convex_hull(const std::vector<vertex> & points): convex_hull(points.data(), points.size()) {}

convex_hull::convex_hull(const vertex * points, const size_t & n_points)
{
	andrew_algorithm(points, n_points);
}

convex_hull::operator const std::vector<index_t> & ()
{
	return CH;
}

///< Andrew's Convex Hull Algorithm: Competitive Programming 4
void convex_hull::andrew_algorithm(const vertex * points, const size_t & n_points)
{
	std::vector<index_t> idx(n_points);
	std::iota(begin(idx), end(idx), 0);

	std::sort(begin(idx), end(idx),
		[&points](const index_t & i, const index_t & j)
		{
			return points[i] < points[j];
		});

	CH.resize(2 * n_points);

	index_t k = 0;
	for(index_t p = 0; p < n_points; ++p)
	{
		const index_t & i = idx[p];
		while(k > 1 && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[i])) --k;
		CH[k++] = i;
	}

	index_t t = k;
	for(index_t p = n_points - 2; p > 0; --p)
	{
		const index_t & i = idx[p];
		while(k > t && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[i])) --k;
		CH[k++] = i;
	}

	while(k > t && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[idx[0]])) --k;

	CH.resize(k);
}

bool convex_hull::ccw(const vertex & p, const vertex & q, const vertex & r)
{
	// TODO vec2
	return ((q - p) * (r - p)).z() > 0;
}


} // namespace gproshan

