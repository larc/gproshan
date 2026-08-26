#include <gproshan/geometry/convex_hull.h>

#include <algorithm>
#include <numeric>


// geometry processing and shape analysis framework
namespace gproshan {


convex_hull::convex_hull(const std::vector<vertex> & points, const float precision): convex_hull(points.data(), size(points), precision) {}

convex_hull::convex_hull(const vertex * points, const size_t n_points, const float precision)
{
	std::vector<ivec2> points2d(n_points);

	#pragma omp parallel for
	for(index_t i = 0; i < n_points; ++i)
	{
		const vertex & p = points[i] * precision;
		points2d[i] = {int(p.x()), int(p.y())};
	}
	andrew_algorithm(points2d);
}

convex_hull::operator const std::vector<index_t> & () const
{
	return CH;
}

///< Andrew's Convex Hull Algorithm: Competitive Programming 4
void convex_hull::andrew_algorithm(const std::vector<ivec2> & points)
{
	std::vector<index_t> idx(size(points));
	std::iota(begin(idx), end(idx), 0);

	std::ranges::sort(idx, [&points](const index_t i, const index_t j)
					{
						return points[i] < points[j];
					});

	CH.resize(2 * size(points));

	index_t k = 0;
	for(index_t p = 0; p < size(points); ++p)
	{
		const index_t i = idx[p];
		while(k > 1 && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[i])) --k;
		CH[k++] = i;
	}

	index_t t = k;
	for(index_t p = size(points) - 2; p > 0; --p)
	{
		const index_t i = idx[p];
		while(k > t && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[i])) --k;
		CH[k++] = i;
	}

	while(k > t && !ccw(points[CH[k - 2]], points[CH[k - 1]], points[idx[0]])) --k;

	CH.resize(k);
}

bool convex_hull::ccw(const ivec2 & p, const ivec2 & q, const ivec2 & r)
{
	return cross(q - p, r - p) > 0;
}


} // namespace gproshan

