#ifndef VERTEX_H
#define VERTEX_H

#include <gproshan/include.h>

#include <iostream>


#define glm_vec3(v) glm::vec3((v)[0], (v)[1], (v)[2])


// geometry processing and shape analysis framework
namespace gproshan {


/*!
	The vertex class represents a 3D point and implements 3D vector operations.
*/
class vertex
{
	public:
		real_t x;
		real_t y;
		real_t z;

	public:
		vertex(const real_t & x_ = 0, const real_t & y_ = 0, const real_t & z_ = 0);
		~vertex() = default;

		real_t & operator [] (const index_t & i);
		const real_t & operator [] (const index_t & i) const;

		vertex unit() const;
		real_t operator * () const;						// norm
		real_t operator , (const vertex & v) const;		// dot product

		vertex operator * (const vertex & v) const;		// cross product
		vertex operator * (const real_t & v) const;		// scalar product
		vertex operator / (const real_t & v) const;		// scalar division
		vertex operator + (const vertex & v) const;
		vertex operator - (const vertex & v) const;
		vertex operator - () const;

		void operator *= (const real_t & v);			// scalar produc
		void operator /= (const real_t & v);
		void operator += (const vertex & v);
		void operator -= (const vertex & v);

		bool operator < (const vertex & v) const;
		bool operator == (const vertex & v) const;

		bool is_zero();
};

vertex operator * (const real_t & a, const vertex & v);

std::ostream & operator << (std::ostream & os, const vertex & v);
std::istream & operator >> (std::istream & is, vertex & v);


} // namespace gproshan

#endif // VERTEX_H

