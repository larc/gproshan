#ifndef VEC_H
#define VEC_H

#include <gproshan/include.h>

#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>

#ifdef __CUDACC__
	#define __host_device__ __host__ __device__
#else
	#define __host_device__
#endif // __CUDACC__


// geometry processing and shape analysis framework
namespace gproshan {


/*!
	The vec<T, N> class represents a generic N dimensional vector and implements the vector operations.
*/
template<class T, size_t N>
class vec
{
	public:
		T values[N] = {};

	public:
		__host_device__
		vec(const T & val = 0)
		{
			for(T & v: values)
				v = val;
		}

		__host_device__
		vec(const std::initializer_list<T> & list)
		{
			int i = -1;
			for(const T & v: list)
				values[++i] = v;
		}

		__host_device__
		vec(const vec<T, N - 1> & v, const T & val = 0)
		{
			for(index_t i = 0; i < N - 1; ++i)
				values[i] = v[i];
			values[N - 1] = val;
		}

		__host_device__
		vec(const vec<T, N + 1> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] = v[i];
		}

		__host_device__
		T & operator [] (const index_t & i)
		{
			assert(i < N);
			return values[i];
		}

		__host_device__
		const T & operator [] (const index_t & i) const
		{
			assert(i < N);
			return values[i];
		}

		///< concatenate with comma operator
		__host_device__
		vec<T, N + 1> operator , (const T & a) const
		{
			return {*this, a};
		}

		__host_device__
		T & x()
		{
			assert(N > 0);
			return values[0];
		}

		__host_device__
		const T & x() const
		{
			assert(N > 0);
			return values[0];
		}

		__host_device__
		T & y()
		{
			assert(N > 1);
			return values[1];
		}

		__host_device__
		const T & y() const
		{
			assert(N > 1);
			return values[1];
		}

		__host_device__
		T & z()
		{
			assert(N > 2);
			return values[2];
		}

		__host_device__
		const T & z() const
		{
			assert(N > 2);
			return values[2];
		}


		///< norm
		__host_device__
		T norm() const
		{
			T res = 0;
			for(const T & v: values)
				res += v * v;
			return std::sqrt(res);
		}

		///< length
		__host_device__
		T length() const
		{
			return norm();
		}

		///< element wise product & scalar product
		__host_device__
		vec<T, N> operator * (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] * v[i];
			return res;
		}

		///< scalar division
		__host_device__
		vec<T, N> operator / (const T & a) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] / a;
			return res;
		}

		///< sum of vectors
		__host_device__
		vec<T, N> operator + (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] + v[i];
			return res;
		}

		///< difference of vectors
		__host_device__
		vec<T, N> operator - (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] - v[i];
			return res;
		}

		///< negative of vector
		__host_device__
		vec<T, N> operator - () const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = - values[i];
			return res;
		}

		///< element wise product self assign
		__host_device__
		const vec<T, N> & operator *= (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] *= v[i];
			return *this;
		}

		///< scalar division self assign
		__host_device__
		const vec<T, N> & operator /= (const T & a)
		{
			for(T & v: values)
				v /= a;
			return *this;
		}

		///< sum of vectors self assign
		__host_device__
		const vec<T, N> & operator += (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] += v[i];
			return *this;
		}

		///< difference of vectors self assign
		__host_device__
		const vec<T, N> & operator -= (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] -= v[i];
			return *this;
		}

		///< comparison less than
		__host_device__
		bool operator < (const vec<T, N> & v) const
		{
			for(index_t i = 0; i < N; ++i)
				if(values[i] != v[i])
					return values[i] < v[i];

			return false;
		}

		///< comparison equal than
		__host_device__
		bool operator == (const vec<T, N> & v) const
		{
			for(index_t i = 0; i < N; ++i)
				if(values[i] != v[i])
					return false;

			return true;
		}

		__host_device__
		bool is_zero()
		{
			T eps = std::numeric_limits<T>::epsilon();
			for(index_t i = 0; i < N; ++i)
				if(abs(values[i]) > eps)
					return false;

			return true;
		}
};


///< scalar product
template<class T, size_t N>
__host_device__
vec<T, N> operator * (const T & a, const vec<T, N> & v)
{
	return v * a;
}


///< cross product
template<class T>
__host_device__
vec<T, 3> cross(const vec<T, 3> & u, const vec<T, 3> & v)
{
	const T & ux = u[0];
	const T & uy = u[1];
	const T & uz = u[2];
	const T & vx = v[0];
	const T & vy = v[1];
	const T & vz = v[2];

	return {uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx};
}

///< dot product
template<class T, size_t N>
__host_device__
T dot(const vec<T, N> & u, const vec<T, N> & v)
{
	T res = 0;
	for(index_t i = 0; i < N; ++i)
		res += u[i] * v[i];
	return res;
}

///< norm
template<class T, size_t N>
__host_device__
T norm(const vec<T, N> & v)
{
	return v.norm();
}

///< length
template<class T, size_t N>
__host_device__
T length(const vec<T, N> & v)
{
	return v.length();
}

///< normalize vector: divide by its norm
template<class T, size_t N>
__host_device__
vec<T, N> normalize(const vec<T, N> & v)
{
	return v / norm(v);
}

template<class T, size_t N>
std::ostream & operator << (std::ostream & os, const vec<T, N> & v)
{
	os << v[0];
	for(index_t i = 1; i < N; ++i)
		os << " " << v[i];
	return os;
}

template<class T, size_t N>
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

