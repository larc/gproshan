#include "vertex.h"

#include <cmath>

vertex::vertex(const vertex_t & x_, const vertex_t & y_, const vertex_t & z_)
{
	x = x_;
	y = y_;
	z = z_;
}

vertex::~vertex()
{

}

vertex vertex::operator*(const vertex & v) const
{
	return vertex(y * v.z - z * v.y, -(x * v.z - z * v.x), x * v.y - y * v.x);
}

void vertex::operator*=(const vertex_t & a)
{
	x *= a;
	y *= a;
	z *= a;
}

vertex_t vertex::operator*() const
{
	return sqrt(x * x + y * y + z * z);
}

vertex vertex::operator/(const vertex_t & a) const
{
	return (1. / a) * (*this);
}

void vertex::operator/=(const vertex_t & a)
{
	(*this) *= (1. / a);
}

vertex_t vertex::operator,(const vertex & v) const
{
	return x * v.x + y * v.y + z * v.z;
}

vertex vertex::operator+(const vertex & v) const
{
	return vertex(x+v.x, y+v.y, z+v.z);
}

void vertex::operator+=(const vertex & v)
{
	x += v.x;
	y += v.y;
	z += v.z;
}

vertex vertex::operator-(const vertex & v) const
{
	return vertex(x-v.x, y-v.y, z-v.z);
}

void vertex::operator-=(const vertex & v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
}

vertex vertex::operator-() const
{
	return vertex(-x, -y, -z);
}

vertex vertex::unit() const
{
	return *this / **this;
}

vertex_t & vertex::operator[](const index_t & i)
{
	return (&x)[i];
}

vertex operator*(const vertex_t & a, const vertex & v)
{
	return vertex(a * v.x, a * v.y, a * v.z);
}

bool operator < (const vertex & a, const vertex & b)
{
	if(a.x != b.x) return a.x < b.x;
	if(a.y != b.y) return a.y < b.y;
	return a.z < b.z;
}

bool operator == (const vertex & a, const vertex & b)
{
	return a.x == b.x && a.y == b.y && a.z == b.z;
}

ostream & operator<<(ostream & os, const vertex & v)
{
	os << v.x << " " << v.y << " " << v.z;
	return os;
}

istream & operator>>(istream & is, vertex & v)
{
	is >> v.x >> v.y >> v.z;
	return is;
}

