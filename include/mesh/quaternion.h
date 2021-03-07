#ifndef QUATERNION_H
#define QUATERNION_H

#include "mesh/vertex.h"

#include <ostream>


// geometry processing and shape analysis framework
namespace gproshan {


class quaternion
{
	protected:
		real_t s;
		vertex v;

	public:
		quaternion(real_t s = 0, real_t vi = 0, real_t vj = 0, real_t vk = 0);
		quaternion(real_t s, const vertex & v);
		quaternion(const vertex & v);

		const quaternion & operator = (real_t s);
		const quaternion & operator = (const vertex & v);
		real_t & operator [] (int index);
		const real_t & operator [] (int index) const;
		void toMatrix(real_t Q[4][4]) const;
		real_t & re(void);
		const real_t & re(void) const;
		vertex & im(void);
		const vertex & im(void) const;

		quaternion operator + (const quaternion & q) const;
		quaternion operator - (const quaternion & q) const;
		quaternion operator - (void) const;
		quaternion operator * (real_t c) const;
		quaternion operator / (real_t c) const;
		void operator += (const quaternion & q);
		void operator += (real_t c);
		void operator -= (const quaternion & q);
		void operator -= (real_t c);
		void operator *= (real_t c);
		void operator /= (real_t c);
		quaternion operator * (const quaternion & q) const;
		void operator *= (const quaternion & q);

		quaternion conj() const;
		quaternion inv() const;
		real_t norm() const;
		real_t norm2() const;
		quaternion unit() const;
		void normalize();

	friend std::ostream & operator << (std::ostream & os, const quaternion & q);
	friend std::istream & operator >> (std::istream & is, quaternion & q);
};

quaternion operator * (real_t c, const quaternion & q);


} // namespace gproshan

#endif // QUATERNION_H

