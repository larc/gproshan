#ifndef VEC_H
#define VEC_H

#include <gproshan/include.h>

#include <cstring>
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
	protected:
		T values[N] = {};

	public:
		vec() = default;

		vec(const std::initializer_list<T> & list)
		{
			memcpy(values, list.begin(), sizeof(values));
		}

		virtual ~vec() = default;

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

		///< norm or length
		T operator * () const
		{
			T res = 0;
			for(const T & v: values)
				res += v * v;
			return std::sqrt(res);
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
		const vec<T, N> & operator *= (const T & a) const
		{
			for(T & v: values)
				v *= a;
			return *this;
		}

		///< scalar division self assign
		const vec<T, N> & operator /= (const T & a) const
		{
			for(T & v: values)
				v /= a;
			return *this;
		}

		///< sum of vectors self assign
		const vec<T, N> & operator += (const vec<T, N> & v) const
		{
			for(index_t i = 0; i < N; ++i)
				values[i] += v[i];
			return *this;
		}

		///< difference of vectors self assign
		const vec<T, N> & operator -= (const vec<T, N> & v) const
		{
			for(index_t i = 0; i < N; ++i)
				values[i] -= v[i];
			return *this;
		}

		/*

		bool operator < (const vec<T, N> & v) const;
		bool operator == (const vec<T, N> & v) const;

		bool is_zero();*/
};

///< scalar product
template<class T, size_t N>
vec<T, N> operator * (const double & a, const vec<T, N> & v)
{
	return v * a;
}

///< dot product
template<class T, size_t N>
T dot(const vec<T, N> & u, const vec<T, N> & v)
{
	return (u, v);
}

///< norm or length
template<class T, size_t N>
T norm(const vec<T, N> & v)
{
	return *v;
}

///< normalize vector: divide by its norm
template<class T, size_t N>
vec<T, N> normalize(const vec<T, N> & v)
{
	return v / *v;
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

