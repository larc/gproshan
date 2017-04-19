#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <set>
#include <vector>

#include "include.h"
#include "vertex.h"
#include "face.h"

typedef set<size_t> ring_t;
typedef vector<size_t> face_t;

class mesh
{
	public:
		mesh();
		virtual size_t get_nvertices() = 0;
		virtual size_t get_nfaces() = 0;
		virtual vertex & operator()(size_t i) = 0;
		virtual face & operator[](size_t i) = 0;
		virtual ring_t & get_rings(size_t i) = 0;
		virtual face_t & get_faces(size_t i) = 0;
		virtual void save(string file) = 0;
};

#endif // MESH_H
