// -----------------------------------------------------------------------------
// libDDG -- Quaternion.h
// -----------------------------------------------------------------------------
//
// Quaternion represents an element of the quaternions, along with all the usual
// vertexs space operations (addition, multiplication by scalars, etc.).  The
// Hamilton product is expressed using the * operator:
//
//	 Quaternion p, q, r;
//	 r = q * p;
//
// and conjugation is expressed using the method Quaternion::conj():
//
//	 Quaternion q;
//	 vertex_t normQSquared = -q.conj()*q;
//
// Individual components can be accessed in several ways: the real and imaginary
// parts can be accessed using the methods Quaternion::re() and Quaternion::im():
//
//	Quaternion q;
//	vertex_t a = q.re();
//	vertex b = q.im();
//
// or by index:
//
//	Quaternion q;
//	vertex_t a  = q[0];
//	vertex_t bi = q[1];
//	vertex_t bj = q[2];
//	vertex_t bk = q[3];
//

#ifndef VIEWER_QUATERNION_H
#define VIEWER_QUATERNION_H

#include "vertex.h"
#include <ostream>

namespace DDG
{
	class Quaternion
	{
		public:
			Quaternion( void );
			// initializes all components to zero

			Quaternion( const Quaternion& q );
			// initializes from existing quaternion

			Quaternion( vertex_t s, vertex_t vi, vertex_t vj, vertex_t vk );
			// initializes with specified real (s) and imaginary (v) components

			Quaternion( vertex_t s, const vertex& v );
			// initializes with specified real (s) and imaginary (v) components

			Quaternion( vertex_t s );
			// initializes purely real quaternion with specified real (s) component (imaginary part is zero)

			Quaternion( const vertex& v );
			// initializes purely imaginary quaternion with specified imaginary (v) component (real part is zero)

			const Quaternion& operator=( vertex_t s );
			// assigns a purely real quaternion with real value s

			const Quaternion& operator=( const vertex& v );
			// assigns a purely real quaternion with imaginary value v

			vertex_t& operator[]( int index );
			// returns reference to the specified component (0-based indexing: r, i, j, k)

			const vertex_t& operator[]( int index ) const;
			// returns const reference to the specified component (0-based indexing: r, i, j, k)

			void toMatrix( vertex_t Q[4][4] ) const;
			// builds 4x4 matrix Q representing (left) quaternion multiplication

			vertex_t& re( void );
			// returns reference to vertex_t part

			const vertex_t& re( void ) const;
			// returns const reference to vertex_t part

			vertex& im( void );
			// returns reference to imaginary part

			const vertex& im( void ) const;
			// returns const reference to imaginary part

			Quaternion operator+( const Quaternion& q ) const;
			// addition

			Quaternion operator-( const Quaternion& q ) const;
			// subtraction

			Quaternion operator-( void ) const;
			// negation

			Quaternion operator*( vertex_t c ) const;
			// right scalar multiplication

			Quaternion operator/( vertex_t c ) const;
			// scalar division

			void operator+=( const Quaternion& q );
			// addition / assignment

			void operator+=( vertex_t c );
			// addition / assignment of pure real

			void operator-=( const Quaternion& q );
			// subtraction / assignment

			void operator-=( vertex_t c );
			// subtraction / assignment of pure real

			void operator*=( vertex_t c );
			// scalar multiplication / assignment

			void operator/=( vertex_t c );
			// scalar division / assignment

			Quaternion operator*( const Quaternion& q ) const;
			// Hamilton product

			void operator*=( const Quaternion& q );
			// Hamilton product / assignment

			Quaternion conj( void ) const;
			// conjugation

			Quaternion inv( void ) const;
			// inverse
			
			vertex_t norm( void ) const;
			// returns Euclidean length

			vertex_t norm2( void ) const;
			// returns Euclidean length squared

			Quaternion unit( void ) const;
			// returns unit quaternion

			void normalize( void );
			// divides by Euclidean length
	
		protected:
			vertex_t s;
			// scalar (vertex_t) part

			vertex v;
			// vertex (imaginary) part
	};
	
	Quaternion operator*( vertex_t c, const Quaternion& q );
	// left scalar multiplication
	
	std::ostream& operator<<( std::ostream& os, const Quaternion& q );
	// prints components
}

#endif

