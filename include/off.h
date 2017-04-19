#ifndef OFF_H
#define OFF_H

#include "mesh.h"

using namespace std;

class off : public mesh
{
	private:
		vertex * vertices;
		face * faces;
		ring_t * rings;
		size_t n_vertices;
		size_t n_faces;
		size_t value;
		face_t * ring_faces;

	public:
		off(string file);
		off(size_t v = 0, size_t f = 0);
		off(vector<vertex> & v_vertices, vector<face> & v_faces);
		~off();
		size_t get_nvertices();
		size_t get_nfaces();
		vertex & operator()(size_t i);
		face & operator[](size_t i);
		ring_t & get_rings(size_t i);
		face_t & get_faces(size_t i);
		void save(string file);
		void generate_grid(size_t s, size_t f, vertex_t r = 1);

	private:
		void set_rings();

	public:
		friend ostream & operator<<(ostream & os, off & o);
		friend istream & operator>>(istream & is, off & o);
};

#endif // OFF_H
