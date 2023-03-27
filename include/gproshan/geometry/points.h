#ifndef POINTS_H
#define POINTS_H

#include <gproshan/geometry/vec.h>


// geometry processing and shape analysis framework
namespace gproshan {


template <class T, size_t N>
std::vector<vec<T, N> > sampling_4points(const size_t & n, const vec<T, N> & a, const vec<T, N> & b, const vec<T, N> & c, const vec<T, N> & d)
{
	std::vector<vec<T, N> > points;
	points.reserve((n + 1) * (n + 1));

	/*
		a --- b
		|     |
		d --- c
	*/

	const real_t alpha = 1.0 / n;
	for(index_t i = 0; i <= n; ++i)
	for(index_t j = 0; j <= n; ++j)
		points.push_back(	j * alpha * (i * alpha * a + (1.0 - i * alpha) * d) + 
					(1.0 - j * alpha) * (i * alpha * b + (1.0 - i * alpha) * c)
					);

	return points;
}


} // namespace gproshan

#endif // POINTS_H

