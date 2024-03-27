#ifndef QUATERNION_H
#define QUATERNION_H

#include <gproshan/geometry/vec.h>

#include <ostream>


// geometry processing and shape analysis framework
namespace gproshan {


using vertex = vec3;


class quaternion
{
	public:
		float s;
		vertex v;

	public:
		quaternion(float s = 0, float vi = 0, float vj = 0, float vk = 0);
		quaternion(float s, const vertex & v);
		quaternion(const vertex & v);

		operator const vertex & () const;
		const quaternion & operator = (float s);
		const quaternion & operator = (const vertex & v);
		float & operator [] (int index);
		float operator [] (int index) const;
		float & re(void);
		float re(void) const;
		vertex & im(void);
		const vertex & im(void) const;

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

