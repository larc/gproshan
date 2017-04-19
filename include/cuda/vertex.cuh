#ifndef VERTEX_CUH
#define VERTEX_CUH

#include <cmath>

#include "../include.h"

struct vertex_cu
{
	vertex_t x;
	vertex_t y;
	vertex_t z;

	__host__ __device__
	vertex_cu(const vertex_t & x_ = 0, const vertex_t & y_ = 0, const vertex_t & z_ = 0)
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
	vertex_t operator*()
	{
		return sqrt(x * x + y * y + z * z);
	}

	__host__ __device__
	vertex_cu operator/(const vertex_t & a) const
	{
		return vertex_cu(x / a, y / a, z / a);
	}

	__host__ __device__
	void operator/=(const vertex_t & v)
	{
		x /= v;
		y /= v;
		z /= v;
	}

	/**
	cross product
	*/
	__host__ __device__
	vertex_t operator,(const vertex_cu & v) const
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
	vertex_t & operator[](const int & i)
	{
		return (&x)[i];
	}
};

__host__ __device__
vertex_cu operator*(const vertex_t & a, const vertex_cu & v);

#endif

