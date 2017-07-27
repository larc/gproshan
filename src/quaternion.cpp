#include <cmath>
#include <iostream>

using namespace std;

#include "quaternion.h"

quaternion :: quaternion(void)
: s(0.), v(0., 0., 0.)
{
}

quaternion :: quaternion(const quaternion & q)
: s(q.s), v(q.v)
{
}

quaternion :: quaternion(vertex_t s_, vertex_t vi, vertex_t vj, vertex_t vk)
: s(s_), v(vi, vj, vk)
{
}

quaternion :: quaternion(vertex_t s_, const vertex & v_)
: s(s_), v(v_)
{
}

quaternion :: quaternion(vertex_t s_)
: s(s_),v(0., 0., 0.)
{
}

quaternion :: quaternion(const vertex & v_)
: s(0.), v(v_)
{
}

const quaternion & quaternion :: operator=(vertex_t _s)
{
	s = _s;
	v = vertex(0., 0., 0.);

	return *this;
}

const quaternion & quaternion :: operator=(const vertex & _v)
{
	s = 0.;
	v = _v;

	return *this;
}


vertex_t & quaternion::operator[](int index)
{
	return ( &s)[ index ];
}

const vertex_t & quaternion::operator[](int index) const
{
	return ( &s)[ index ];
}

void quaternion::toMatrix(vertex_t Q[4][4]) const
{
	Q[0][0] =	s; Q[0][1] = -v.x; Q[0][2] = -v.y; Q[0][3] = -v.z;
	Q[1][0] = v.x; Q[1][1] =	 s; Q[1][2] = -v.z; Q[1][3] =	v.y;
	Q[2][0] = v.y; Q[2][1] =	v.z; Q[2][2] =	 s; Q[2][3] = -v.x;
	Q[3][0] = v.z; Q[3][1] = -v.y; Q[3][2] =	v.x; Q[3][3] =	 s;
}

vertex_t & quaternion::re(void)
{
	return s;
}

const vertex_t & quaternion::re(void) const
{
	return s;
}

vertex & quaternion::im(void)
{
	return v;
}

const vertex & quaternion::im(void) const
{
	return v;
}

quaternion quaternion::operator+(const quaternion & q) const
{
	return quaternion(s + q.s, v + q.v);
}

quaternion quaternion::operator-(const quaternion & q) const
{
	return quaternion(s - q.s, v - q.v);
}

quaternion quaternion::operator-(void) const
{
	return quaternion(-s, -v);
}

quaternion quaternion::operator*(vertex_t c) const
{
	return quaternion(c * s, c * v);
}

quaternion operator*(vertex_t c, const quaternion & q)
{
	return q * c;
}

quaternion quaternion::operator/(vertex_t c) const
{
	return quaternion(s / c, v / c);
}

void quaternion::operator+=(const quaternion & q)
{
	s += q.s;
	v += q.v;
}

void quaternion::operator+=(vertex_t c)
{
	s += c;
}

void quaternion::operator-=(const quaternion & q)
{
	s -= q.s;
	v -= q.v;
}

void quaternion::operator-=(vertex_t c)
{
	s -= c;
}

void quaternion::operator*=(vertex_t c)
{
	s *= c;
	v *= c;
}

void quaternion::operator/=(vertex_t c)
{
	s /= c;
	v /= c;
}

// Hamilton product
quaternion quaternion::operator*(const quaternion & q) const
{
	const vertex_t & s1(s);
	const vertex_t & s2(q.s);
	const vertex & v1(v);
	const vertex & v2(q.v);

	return quaternion(s1*s2 - (v1,v2), s1*v2 + s2*v1 + (v1*v2));
}

void quaternion::operator*=(const quaternion & q)
{
	*this = (*this * q);
}

quaternion quaternion::conj(void) const
{
	return quaternion(s, -v);
}

quaternion quaternion::inv(void) const
{
	return (this->conj()) / this->norm2();
}

vertex_t quaternion::norm(void) const
{
	return sqrt(norm2());
}

vertex_t quaternion::norm2(void) const
{
	return s * s + (v , v);
}

quaternion quaternion::unit(void) const
{
	return *this / norm();
}

void quaternion::normalize(void)
{
	*this /= norm();
}
	
// spherical-linear interpolation
quaternion slerp(const quaternion & q0, const quaternion & q1, vertex_t t)
{
	// interpolate length
	vertex_t m0 = q0.norm();
	vertex_t m1 = q1.norm();
	vertex_t m = (1-t)*m0 + t*m1;

	// interpolate direction
	quaternion p0 = q0 / m0;
	quaternion p1 = q1 / m1;
	vertex_t theta = acos((p0.conj() * p1).re());
	quaternion p = (sin((1 - t) * theta) * p0 + sin(t * theta) * p1) / sin(theta);

	return m * p;
}

std::ostream & operator<<(std::ostream & os, const quaternion & q)
{
	return os << "(" << q.re() << ", " << q.im() << ")";	
}

