#include <gproshan/geometry/quaternion.h>

#include <cmath>
#include <iostream>


// geometry processing and shape analysis framework
namespace gproshan {


quaternion::quaternion(float s, const vec3 & v): s(s), v(v) {}

quaternion::quaternion(const vec3 & v): v(v) {}

quaternion::operator const vec3 & () const
{
	return v;
}

float & quaternion::operator [] (int index)
{
	return v[index];
}

float quaternion::operator [] (int index) const
{
	return v[index];
}

quaternion quaternion::operator + (const quaternion & q) const
{
	return {s + q.s, v + q.v};
}

quaternion quaternion::operator - (const quaternion & q) const
{
	return {s - q.s, v - q.v};
}

quaternion quaternion::operator - () const
{
	return {-s, -v};
}

quaternion quaternion::operator * (float c) const
{
	return {c * s, c * v};
}

quaternion quaternion::operator / (float c) const
{
	return {s / c, v / c};
}

// Hamilton product
quaternion quaternion::operator * (const quaternion & q) const
{
	return {s * q.s - dot(v, q.v), s * q.v + q.s * v + cross(v, q.v)};
}

quaternion quaternion::conj() const
{
	return {s, -v};
}

quaternion quaternion::inv() const
{
	return conj() / norm2();
}

float quaternion::norm() const
{
	return sqrt(norm2());
}

float quaternion::norm2() const
{
	return s * s + dot(v, v);
}


float norm(const quaternion & q)
{
	return q.norm();
}

quaternion normalize(const quaternion & q)
{
	return q / norm(q);
}

quaternion operator * (float c, const quaternion & q)
{
	return q * c;
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

