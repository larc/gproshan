#ifndef VEC_H
#define VEC_H

#include <gproshan/include.h>

#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>


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
		vec() = default;

		vec(const std::initializer_list<T> & list)
		{
			memcpy(values, list.begin(), sizeof(values));
		}

		vec(const T & val)
		{
			for(T & v: values)
				v = val;
		}

		const vec<T, N> & operator = (const T & val)
		{
			for(T & v: values)
				v = val;
			return *this;
		}

		T & operator [] (const index_t & i)
		{
			assert(i < N);
			return values[i];
		}

		const T & operator [] (const index_t & i) const
		{
			assert(i < N);
			return values[i];
		}

		///< norm
		T norm() const
		{
			T res = 0;
			for(const T & v: values)
				res += v * v;
			return std::sqrt(res);
		}

		///< length
		T length() const
		{
			return norm();
		}

		///< dot product
		T operator , (const vec<T, N> & v) const
		{
			T res = 0;
			for(index_t i = 0; i < N; ++i)
				res += values[i] * v[i];
			return res;
		}

		///< scalar product
		vec<T, N> operator * (const T & a) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] * a;
			return res;
		}

		///< scalar division
		vec<T, N> operator / (const T & a) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] / a;
			return res;
		}

		///< sum of vectors
		vec<T, N> operator + (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] + v[i];
			return res;
		}

		///< difference of vectors
		vec<T, N> operator - (const vec<T, N> & v) const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = values[i] - v[i];
			return res;
		}

		///< negative of vector
		vec<T, N> operator - () const
		{
			vec<T, N> res;
			for(index_t i = 0; i < N; ++i)
				res[i] = - values[i];
			return res;
		}

		///< scalar product self assign
		const vec<T, N> & operator *= (const T & a)
		{
			for(T & v: values)
				v *= a;
			return *this;
		}

		///< scalar division self assign
		const vec<T, N> & operator /= (const T & a)
		{
			for(T & v: values)
				v /= a;
			return *this;
		}

		///< sum of vectors self assign
		const vec<T, N> & operator += (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] += v[i];
			return *this;
		}

		///< difference of vectors self assign
		const vec<T, N> & operator -= (const vec<T, N> & v)
		{
			for(index_t i = 0; i < N; ++i)
				values[i] -= v[i];
			return *this;
		}

		///< comparison less than
		bool operator < (const vec<T, N> & v) const
		{
			if(x != v.x) return x < v.x;
			if(y != v.y) return y < v.y;
			return z < v.z;
		}

		///< comparison equal than
		bool operator == (const vec<T, N> & v) const
		{
			return x == v.x && y == v.y && z == v.z;
		}

		bool is_zero()
		{
			double eps = std::numeric_limits<double>::epsilon();

			return abs(double(x)) < eps && abs(double(y)) < eps && abs(double(z)) < eps;
		}
};


///< scalar product
template<class T, size_t N>
vec<T, N> operator * (const double & a, const vec<T, N> & v)
{
	return v * a;
}

///< cross product
template<class T>
vec<T, 3> operator * (const vec<T, 3> & u, const vec<T, 3> & v)
{
	return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

///< cross product
template<class T>
vec<T, 3> cross(const vec<T, 3> & u, const vec<T, 3> & v)
{
	return u * v;
}

///< dot product
template<class T, size_t N>
T dot(const vec<T, N> & u, const vec<T, N> & v)
{
	return (u, v);
}

///< norm
template<class T, size_t N>
T norm(const vec<T, N> & v)
{
	return v.norm();
}

///< length
template<class T, size_t N>
T length(const vec<T, N> & v)
{
	return v.length();
}

///< normalize vector: divide by its norm
template<class T, size_t N>
vec<T, N> normalize(const vec<T, N> & v)
{
	return v / norm(v);
}

///< std ostream
template<class T, size_t N>
std::ostream & operator << (std::ostream & os, const vec<T, N> & v)
{
	for(index_t i = 0; i < N - 1; ++i)
		os << v[i] << " ";
	return os << v[N - 1];
}

///< std istream
template<class T, size_t N>
std::istream & operator >> (std::istream & is, vec<T, N> & v)
{
	for(index_t i = 0; i < N; ++i)
		is >> v[i];
	return is;
}


} // namespace gproshan

#endif // VEC_H

