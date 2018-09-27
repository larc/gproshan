#ifndef VERTEX_CUH
#define VERTEX_CUH

#include <cmath>

#include "../include.h"

struct vertex_cu
{
	real_t x;
	real_t y;
	real_t z;

	__host__ __device__
	vertex_cu(const real_t & x_ = 0, const real_t & y_ = 0, const real_t & z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	__host__ __device__
	~vertex_cu()
	{

	}


	/**
		doc product
	*/
	__host__ __device__
	vertex_cu operator*(const vertex_cu & v) const
	{
		return vertex_cu(y * v.z - z * v.y, -(x * v.z - z * v.x), x * v.y - y * v.x);
	}

	/**
		Norma
	*/
	__host__ __device__
	real_t operator*()
	{
		return sqrt(x * x + y * y + z * z);
	}

	__host__ __device__
	vertex_cu operator/(const real_t & a) const
	{
		return vertex_cu(x / a, y / a, z / a);
	}

	__host__ __device__
	void operator/=(const real_t & v)
	{
		x /= v;
		y /= v;
		z /= v;
	}

	/**
	cross product
	*/
	__host__ __device__
	real_t operator,(const vertex_cu & v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}

	__host__ __device__
	vertex_cu operator+(const vertex_cu & v) const
	{
		return vertex_cu(x+v.x, y+v.y, z+v.z);
	}

	__host__ __device__
	void operator+=(const vertex_cu & v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	__host__ __device__
	vertex_cu operator-(const vertex_cu & v) const
	{
		return vertex_cu(x - v.x, y - v.y, z - v.z);
	}

	__host__ __device__
	void operator-=(const vertex_cu & v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	__host__ __device__
	vertex_cu operator-() const
	{
		return vertex_cu(-x, -y, -z);
	}

	__host__ __device__
	real_t & operator[](const int & i)
	{
		return (&x)[i];
	}
};

__host__ __device__
vertex_cu operator*(const real_t & a, const vertex_cu & v);

#endif

