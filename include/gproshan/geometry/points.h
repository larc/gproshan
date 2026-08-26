#ifndef POINTS_H
#define POINTS_H

#include <gproshan/geometry/vec.h>


// geometry processing and shape analysis framework
namespace gproshan {


template <class T, size_t N>
std::vector<vec<T, N> > sampling_4points(const size_t n, const vec<T, N> & a, const vec<T, N> & b, const vec<T, N> & c, const vec<T, N> & d)
{
	std::vector<vec<T, N> > points;
	points.reserve((n + 1) * (n + 1));

	/*
		a --- b
		|     |
		d --- c
	*/

	const size_t n2 = n * n;
	for(index_t i = 0; i <= n; ++i)
	for(index_t j = 0; j <= n; ++j)
		points.push_back((j * (i * a + (n - i) * d) +
					(n - j) * (i * b + (n - i) * c)) / n2);

	return points;
}


} // namespace gproshan

#endif // POINTS_H

