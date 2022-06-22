#ifndef VEC_H
#define VEC_H

#include <gproshan/include.h>

#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>

#ifndef __CUDACC__
	#define __host__
	#define __device__
#endif

// geometry processing and shape analysis framework
namespace gproshan {


/*!
	The vec<T, N> class represents a generic N dimensional vector and implements the vector operations.
*/
template<class T, size_t N>
class vec
{
	public:
		union
		{
			T values[N] = {};
			struct
			{
				T x;
				T y;
				T z;
			};
		};

	public:
		__host__ __device__
		vec(const std::initializer_list<T> & list)
		{
			int i = -1;
			for(const T & v: list)
				values[++i] = v;
		}

		__host__ __device__
		vec(const vec<T, N - 1> & v, const T & val = 0)
		{
			for(index_t i = 0; i < N - 1; ++i)
				values[i] = v[i];
			values[N - 1] = val;
		}

		__host__ __device__
		vec(const vec<T, N + 1> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] = v[i];
		}

		__host__ __device__
		vec(const T & val = 0)
		{
			for(T & v: values)
				v = val;
		}

		__host__ __device__
		T & operator [] (const index_t & i)
		{
			assert(i < N);
			return values[i];
		}

		__host__ __device__
		const T & operator [] (const index_t & i) const
		{
			assert(i < N);
			return values[i];
		}

		///< norm
		__host__ __device__
		T norm() const
		{
			T res = 0;
			for(const T & v: values)
				res += v * v;
			return std::sqrt(res);
		}

		///< length
		__host__ __device__
		T length() const
		{
			return norm();
		}

		///< dot product
		__host__ __device__
		T operator , (const vec<T, N> & v) const
		{
			T res = 0;
			for(index_t i = 0; i < N; ++i)
				res += values[i] * v[i];
			return res;
		}

		///< scalar product
		__host__ __device__
		vec<T, N> operator * (const T & a) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] * a;
			return res;
		}

		///< scalar division
		__host__ __device__
		vec<T, N> operator / (const T & a) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] / a;
			return res;
		}

		///< sum of vectors
		__host__ __device__
		vec<T, N> operator + (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] + v[i];
			return res;
		}

		///< difference of vectors
		__host__ __device__
		vec<T, N> operator - (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] - v[i];
			return res;
		}

		///< negative of vector
		__host__ __device__
		vec<T, N> operator - () const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = - values[i];
			return res;
		}

		///< scalar product self assign
		__host__ __device__
		const vec<T, N> & operator *= (const T & a)
		{
			for(T & v: values)
				v *= a;
			return *this;
		}

		///< scalar division self assign
		__host__ __device__
		const vec<T, N> & operator /= (const T & a)
		{
			for(T & v: values)
				v /= a;
			return *this;
		}

		///< sum of vectors self assign
		__host__ __device__
		const vec<T, N> & operator += (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] += v[i];
			return *this;
		}

		///< difference of vectors self assign
		__host__ __device__
		const vec<T, N> & operator -= (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] -= v[i];
			return *this;
		}

		///< comparison less than
		__host__ __device__
		bool operator < (const vec<T, N> & v) const
		{
			if(x != v.x) return x < v.x;
			if(y != v.y) return y < v.y;
			return z < v.z;
		}

		///< comparison equal than
		__host__ __device__
		bool operator == (const vec<T, N> & v) const
		{
			return x == v.x && y == v.y && z == v.z;
		}

		__host__ __device__
		bool is_zero()
		{
			double eps = std::numeric_limits<double>::epsilon();

			return abs(double(x)) < eps && abs(double(y)) < eps && abs(double(z)) < eps;
		}
};


///< scalar product
template<class T, size_t N>
__host__ __device__
vec<T, N> operator * (const T & a, const vec<T, N> & v)
{
	return v * a;
}

///< cross product
template<class T>
__host__ __device__
vec<T, 3> operator * (const vec<T, 3> & u, const vec<T, 3> & v)
{
	return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

///< cross product
template<class T>
__host__ __device__
vec<T, 3> cross(const vec<T, 3> & u, const vec<T, 3> & v)
{
	return u * v;
}

///< dot product
template<class T, size_t N>
__host__ __device__
T dot(const vec<T, N> & u, const vec<T, N> & v)
{
	return (u, v);
}

///< norm
template<class T, size_t N>
__host__ __device__
T norm(const vec<T, N> & v)
{
	return v.norm();
}

///< length
template<class T, size_t N>
__host__ __device__
T length(const vec<T, N> & v)
{
	return v.length();
}

///< normalize vector: divide by its norm
template<class T, size_t N>
__host__ __device__
vec<T, N> normalize(const vec<T, N> & v)
{
	return v / norm(v);
}

///< std ostream
template<class T, size_t N>
__host__ __device__
std::ostream & operator << (std::ostream & os, const vec<T, N> & v)
{
	for(index_t i = 0; i < N - 1; ++i)
		os << v[i] << " ";
	return os << v[N - 1];
}

///< std istream
template<class T, size_t N>
__host__ __device__
std::istream & operator >> (std::istream & is, vec<T, N> & v)
{
	for(index_t i = 0; i < N; ++i)
		is >> v[i];
	return is;
}


using vec2 = vec<real_t, 2>;
using vec3 = vec<real_t, 3>;
using vec4 = vec<real_t, 4>;

using uvec2 = vec<unsigned int, 2>;
using uvec3 = vec<unsigned int, 3>;
using uvec4 = vec<unsigned int, 4>;

using ivec2 = vec<int, 2>;
using ivec3 = vec<int, 3>;
using ivec4 = vec<int, 4>;


} // namespace gproshan

#endif // VEC_H

