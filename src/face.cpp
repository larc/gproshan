#include "face.h"

face::face()
{

}

face::~face()
{

}

size_t & face::operator[](size_t i)
{
	return vertices[i];
}

bool face::isVertex(size_t v, size_t & i)
{
	i = NIL;
	if(vertices[0] == v) i = 0;
	if(vertices[1] == v) i = 1;
	if(vertices[2] == v) i = 2;
	return i != NIL;
}

ostream & operator<<(ostream & os, face & f)
{
	size_t n = FSIZE - 1;
	os<<FSIZE<<" ";
	for(size_t i = 0; i < n; i++)
		os<<f[i]<<" ";
	os<<f[n];
	return os;
}

istream & operator>>(istream & is, face & f)
{
	size_t tmp;
	is>>tmp;
	for(size_t i = 0; i < FSIZE; i++)
		is>>f[i];
	return is;
}
