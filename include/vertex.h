#ifndef VERTEX_H
#define VERTEX_H

#include "include.h"

#include <iostream>

using namespace std;

/**
 * @brief The vertex class is a point in the space 3D.
 */
class vertex
{
	public:
		real_t x;
		real_t y;
		real_t z;

	public:
		vertex(const real_t & x_ = 0, const real_t & y_ = 0, const real_t & z_ = 0);
		~vertex();
		vertex operator*(const vertex & v) const;		//cross product
		void operator*=(const real_t & v);			//scalar produc
		real_t operator*() const;						//norm
		vertex operator/(const real_t & v) const;
		void operator/=(const real_t & v);
		real_t operator,(const vertex & v) const;		//dot product
		vertex operator+(const vertex & v) const;
		void operator+=(const vertex & v);
		vertex operator-(const vertex & v) const;
		void operator-=(const vertex & v);
		vertex operator-() const;

		vertex unit() const;
		real_t & operator[](const index_t & i);
};

vertex operator*(const real_t & a, const vertex & v);
bool operator < (const vertex & a, const vertex & b);
bool operator == (const vertex & a, const vertex & b);
ostream & operator<<(ostream & os, const vertex & v);
istream & operator>>(istream & is, vertex & v);

#endif // VERTEX_H
