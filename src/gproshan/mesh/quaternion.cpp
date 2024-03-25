#include <gproshan/mesh/quaternion.h>

#include <cmath>
#include <iostream>


// geometry processing and shape analysis framework
namespace gproshan {


quaternion::quaternion(float s_, float vi, float vj, float vk): s(s_), v{vi, vj, vk} {}

quaternion::quaternion(float s_, const vertex & v_): s(s_), v(v_) {}

quaternion::quaternion(const vertex & v_): s(0), v(v_) {}

quaternion::operator const vertex & () const
{
	return v;
}

const quaternion & quaternion::operator = (float _s)
{
	s = _s;
	v = {0, 0, 0};

	return *this;
}

const quaternion & quaternion::operator = (const vertex & _v)
{
	s = 0;
	v = _v;

	return *this;
}


float & quaternion::operator [] (int index)
{
	return v[index];
}

float quaternion::operator [] (int index) const
{
	return v[index];
}

float & quaternion::re()
{
	return s;
}

float quaternion::re() const
{
	return s;
}

vertex & quaternion::im()
{
	return v;
}

const vertex & quaternion::im() const
{
	return v;
}

quaternion quaternion::operator + (const quaternion & q) const
{
	return quaternion(s + q.s, v + q.v);
}

quaternion quaternion::operator - (const quaternion & q) const
{
	return quaternion(s - q.s, v - q.v);
}

quaternion quaternion::operator - () const
{
	return quaternion(-s, -v);
}

quaternion quaternion::operator * (float c) const
{
	return quaternion(c * s, c * v);
}

quaternion operator * (float c, const quaternion & q)
{
	return q * c;
}

quaternion quaternion::operator / (float c) const
{
	return quaternion(s / c, v / c);
}

void quaternion::operator += (const quaternion & q)
{
	s += q.s;
	v += q.v;
}

void quaternion::operator += (float c)
{
	s += c;
}

void quaternion::operator -= (const quaternion & q)
{
	s -= q.s;
	v -= q.v;
}

void quaternion::operator -= (float c)
{
	s -= c;
}

void quaternion::operator *= (float c)
{
	s *= c;
	v *= c;
}

void quaternion::operator /= (float c)
{
	s /= c;
	v /= c;
}

// Hamilton product
quaternion quaternion::operator * (const quaternion & q) const
{
	const float s1(s);
	const float s2(q.s);
	const vertex & v1(v);
	const vertex & v2(q.v);

	return quaternion(s1 * s2 - dot(v1, v2), s1 * v2 + s2 * v1 + cross(v1, v2));
}

void quaternion::operator *= (const quaternion & q)
{
	*this = (*this * q);
}

quaternion quaternion::conj() const
{
	return quaternion(s, -v);
}

quaternion quaternion::inv() const
{
	return (this->conj()) / this->norm2();
}

float quaternion::norm() const
{
	return sqrt(norm2());
}

float quaternion::norm2() const
{
	return s * s + dot(v, v);
}

quaternion quaternion::unit() const
{
	return *this / norm();
}

void quaternion::normalize()
{
	*this /= norm();
}

// spherical-linear interpolation
quaternion slerp(const quaternion & q0, const quaternion & q1, float t)
{
	// interpolate length
	float m0 = q0.norm();
	float m1 = q1.norm();
	float m = (1-t)*m0 + t*m1;

	// interpolate direction
	quaternion p0 = q0 / m0;
	quaternion p1 = q1 / m1;
	float theta = acos((p0.conj() * p1).re());
	quaternion p = (sin((1 - t) * theta) * p0 + sin(t * theta) * p1) / sin(theta);

	return m * p;
}

std::ostream & operator << (std::ostream & os, const quaternion & q)
{
	return os << q.s << " " << q.v;
}

std::istream & operator >> (std::istream & is, quaternion & q)
{
	return is >> q.s >> q.v;
}

} // namespace gproshan

