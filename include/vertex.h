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
		vertex_t x;
		vertex_t y;
		vertex_t z;

	public:
		vertex(const vertex_t & x_ = 0, const vertex_t & y_ = 0, const vertex_t & z_ = 0);
		~vertex();
		vertex operator*(const vertex & v) const;		//cross product
		void operator*=(const vertex_t & v);			//scalar produc
		vertex_t operator*() const;						//norm
		vertex operator/(const vertex_t & v) const;
		void operator/=(const vertex_t & v);
		vertex_t operator,(const vertex & v) const;		//dot product
		vertex operator+(const vertex & v) const;
		void operator+=(const vertex & v);
		vertex operator-(const vertex & v) const;
		void operator-=(const vertex & v);
		vertex operator-() const;

		vertex unit() const;
		vertex_t & operator[](const index_t & i);
};

vertex operator*(const vertex_t & a, const vertex & v);
bool operator < (const vertex & a, const vertex & b);
bool operator == (const vertex & a, const vertex & b);
ostream & operator<<(ostream & os, const vertex & v);
istream & operator>>(istream & is, vertex & v);

#endif // VERTEX_H
