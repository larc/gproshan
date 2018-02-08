#ifndef VIEWER_QUATERNION_H
#define VIEWER_QUATERNION_H

#include "vertex.h"
#include <ostream>

class quaternion
{
	protected:
		vertex_t s;
		vertex v;

	public:
		quaternion(void);
		quaternion(const quaternion & q);
		quaternion(vertex_t s, vertex_t vi, vertex_t vj, vertex_t vk);
		quaternion(vertex_t s, const vertex & v);
		quaternion(vertex_t s);
		quaternion(const vertex & v);

		const quaternion & operator=(vertex_t s);
		const quaternion & operator=(const vertex & v);
		vertex_t & operator[](int index);
		const vertex_t & operator[](int index) const;
		void toMatrix(vertex_t Q[4][4]) const;
		vertex_t & re(void);
		const vertex_t & re(void) const;
		vertex & im(void);
		const vertex & im(void) const;

		quaternion operator+(const quaternion & q) const;
		quaternion operator-(const quaternion & q) const;
		quaternion operator-(void) const;
		quaternion operator*(vertex_t c) const;
		quaternion operator/(vertex_t c) const;
		void operator+=(const quaternion & q);
		void operator+=(vertex_t c);
		void operator-=(const quaternion & q);
		void operator-=(vertex_t c);
		void operator*=(vertex_t c);
		void operator/=(vertex_t c);
		quaternion operator*(const quaternion & q) const;
		void operator*=(const quaternion & q);

		quaternion conj(void) const;
		quaternion inv(void) const;
		vertex_t norm(void) const;
		vertex_t norm2(void) const;
		quaternion unit(void) const;
		void normalize(void);

};

quaternion operator*(vertex_t c, const quaternion & q);
std::ostream & operator<<(std::ostream & os, const quaternion & q);

#endif

