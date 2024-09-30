#ifndef QUATERNION_H
#define QUATERNION_H

#include <gproshan/geometry/vec.h>

#include <ostream>


// geometry processing and shape analysis framework
namespace gproshan {


class quaternion
{
	public:
		float s;
		vec3 v;

	public:
		quaternion(float s = 0, float vi = 0, float vj = 0, float vk = 0);
		quaternion(float s, const vec3 & v);
		quaternion(const vec3 & v);

		operator const vec3 & () const;
		const quaternion & operator = (float s);
		const quaternion & operator = (const vec3 & v);
		float & operator [] (int index);
		float operator [] (int index) const;
		float & re(void);
		float re(void) const;
		vec3 & im(void);
		const vec3 & im(void) const;

		quaternion operator + (const quaternion & q) const;
		quaternion operator - (const quaternion & q) const;
		quaternion operator - (void) const;
		quaternion operator * (float c) const;
		quaternion operator / (float c) const;
		void operator += (const quaternion & q);
		void operator += (float c);
		void operator -= (const quaternion & q);
		void operator -= (float c);
		void operator *= (float c);
		void operator /= (float c);
		quaternion operator * (const quaternion & q) const;
		void operator *= (const quaternion & q);

		quaternion conj() const;
		quaternion inv() const;
		float norm() const;
		float norm2() const;
		quaternion unit() const;
		void normalize();

	friend std::ostream & operator << (std::ostream & os, const quaternion & q);
	friend std::istream & operator >> (std::istream & is, quaternion & q);
};

quaternion operator * (float c, const quaternion & q);


} // namespace gproshan

#endif // QUATERNION_H

