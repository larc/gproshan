#ifndef FACE_H
#define FACE_H

#include <iostream>

#include "include.h"

#define FSIZE 3

using namespace std;

/**
 * @brief The face class save of indices of vertices.
 */
class face
{
	private:
		size_t vertices[FSIZE];

	public:
		face();
		~face();
		size_t & operator[](size_t i);
		bool isVertex(size_t v, size_t & i);

	public:
		friend ostream & operator<<(ostream & os, face & f);
		friend istream & operator>>(istream & is, face & f);
};

#endif // FACE_H
