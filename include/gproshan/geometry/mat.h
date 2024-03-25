#ifndef MAT_H
#define MAT_H

#include <gproshan/geometry/vec.h>

#ifndef __CUDACC__
	#include <armadillo>
#endif // __CUDACC__


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
		mat<T, N> t() const
		{
			return transpose(*this);
		}

		__host_device__
		static mat<T, N> identity()
		{
			mat<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i][i] = 1;
			return res;
		}
};


template<class T, size_t N>
__host_device__
mat<T, N> transpose(const mat<T, N> & m)
{
	mat<T, N> res;
	for(index_t i = 0; i < N; ++i)
	for(index_t j = 0; j < N; ++j)
		res[i][j] = m[j][i];
	return res;
}


template<class T, size_t N>
std::ostream & operator << (std::ostream & os, const mat<T, N> & m)
{
	for(index_t i = 0; i < N; ++i)
		os << m[i] << std::endl;
	return os;
}

template<class T, size_t N>
std::istream & operator >> (std::istream & is, mat<T, N> & m)
{
	for(index_t i = 0; i < N; ++i)
		is >> m[i];
	return is;
}


#ifndef __CUDACC__
template<class T, size_t N>
mat<T, N> inverse(const mat<T, N> & m)
{
	mat<T, N> inv;
	mat<T, N> mt = m.t();

	arma::Mat<T> ainv((T *) &inv, N, N, false, true);
	arma::Mat<T> amt((T *) &mt, N, N, false, true);

	arma::inv(ainv, amt);

	return inv.t();
}
#endif // __CUDACC__


using mat2 = mat<real_t, 2>;
using mat3 = mat<real_t, 3>;
using mat4 = mat<real_t, 4>;


} // namespace gproshan

#endif // MAT_H

