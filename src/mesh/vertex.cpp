#include "mesh/vertex.h"

#include <cmath>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


vertex::vertex(const real_t & x_, const real_t & y_, const real_t & z_)
{
	x = x_;
	y = y_;
	z = z_;
}

real_t & vertex::operator [] (const index_t & i)
{
	return (&x)[i];
}

const real_t & vertex::operator [] (const index_t & i) const
{
	return (&x)[i];
}

vertex vertex::unit() const
{
	return *this / **this;
}

real_t vertex::operator * () const
{
	return sqrt(x * x + y * y + z * z);
}

real_t vertex::operator , (const vertex & v) const
{
	return x * v.x + y * v.y + z * v.z;
}

vertex vertex::operator * (const vertex & v) const
{
	return vertex(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

vertex vertex::operator / (const real_t & a) const
{
	return (1. / a) * (*this);
}

vertex vertex::operator + (const vertex & v) const
{
	return vertex(x + v.x, y + v.y, z + v.z);
}

vertex vertex::operator - (const vertex & v) const
{
	return vertex(x - v.x, y - v.y, z - v.z);
}

vertex vertex::operator - () const
{
	return vertex(-x, -y, -z);
}

void vertex::operator *= (const real_t & a)
{
	x *= a;
	y *= a;
	z *= a;
}

void vertex::operator /= (const real_t & a)
{
	(*this) *= (1. / a);
}

void vertex::operator += (const vertex & v)
{
	x += v.x;
	y += v.y;
	z += v.z;
}

void vertex::operator -= (const vertex & v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
}

bool vertex::operator < (const vertex & v)
{
	if(x != v.x) return x < v.x;
	if(y != v.y) return y < v.y;
	return z < v.z;
}

bool vertex::operator == (const vertex & v)
{
	return x == v.x && y == v.y && z == v.z;
}

bool vertex::is_zero()
{
	real_t eps = std::numeric_limits<real_t>::epsilon();

	return abs(x) < eps && abs(y) < eps && abs(z) < eps;
}

vertex operator * (const real_t & a, const vertex & v)
{
	return vertex(a * v.x, a * v.y, a * v.z);
}

ostream & operator << (ostream & os, const vertex & v)
{
	os << v.x << " " << v.y << " " << v.z;
	return os;
}

istream & operator >> (istream & is, vertex & v)
{
	is >> v.x >> v.y >> v.z;
	return is;
}


} // namespace gproshan

