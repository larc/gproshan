#ifndef MAT_H
#define MAT_H

#include <gproshan/geometry/vec.h>


// geometry processing and shape analysis framework
namespace gproshan {


template<class T, size_t N>
using row = vec<T, N>;

template<class T, size_t N>
class mat
{
	public:
		row<T, N> rows[N] = {};

	public:
		__host_device__
		mat(const std::initializer_list<row<T, N> > & list = {})
		{
			int i = -1;
			for(const row<T, N> & r: list)
				rows[++i] = r;
		}

		__host_device__
		T & operator () (const index_t i, const index_t j)
		{
			return rows[i][j];
		}

		__host_device__
		T operator () (const index_t i, const index_t j) const
		{
			return rows[i][j];
		}

		__host_device__
		row<T, N> & operator [] (const index_t i)
		{
			assert(i < N);
			return rows[i];
		}

		__host_device__
		const row<T, N> & operator [] (const index_t i) const
		{
			assert(i < N);
			return rows[i];
		}

		__host_device__
		mat<T, N> operator * (const mat<T, N> & b) const
		{
			mat<T, N> res;
			mat<T, N> bt = transpose(b);
			for(index_t i = 0; i < N; ++i)
			for(index_t j = 0; j < N; ++j)
				res[i][j] = dot(rows[i], bt[j]);
			return res;
		}

		__host_device__
		vec<T, N> operator * (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = dot(rows[i], v);
			return res;
		}

		__host_device__
		static mat<T, N> identity()
		{
			mat<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i][i] = 1;
			return res;
		}

		__host_device__
		static mat<T, N> transpose(const mat<T, N> & m)
		{
			mat<T, N> res;
			for(index_t i = 0; i < N; ++i)
			for(index_t j = 0; j < N; ++j)
				res[i][j] = m[j][i];
			return res;
		}
};

///< std std::ostream
template<class T, size_t N>
std::ostream & operator << (std::ostream & os, const mat<T, N> & m)
{
	for(index_t i = 0; i < N; ++i)
		os << m[i] << std::endl;
	return os;
}

///< std std::istream
template<class T, size_t N>
std::istream & operator >> (std::istream & is, mat<T, N> & m)
{
	for(index_t i = 0; i < N; ++i)
		is >> m[i];
	return is;
}


using mat2 = mat<real_t, 2>;
using mat3 = mat<real_t, 3>;
using mat4 = mat<real_t, 4>;

mat4 inverse(const mat4 & m);


} // namespace gproshan

#endif // MAT_H

