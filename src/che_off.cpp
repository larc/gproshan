#include "che_off.h"

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>


che_off::che_off(const string & file)
{
	init(file);
}

che_off::~che_off()
{

}

void che_off::read_file(const string & file)
{
	string soff;
	size_t n_v, n_f, v;

	ifstream is(file);

	assert(is.good());

	is >> soff;
	is >> n_v >> n_f >> v;
	init(n_v, n_f);

	int r, g, b, a;
	for(index_t i = 0; i < n_vertices_; i++)
	{
		is >> GT[i];
		if(soff[0] == 'C') // COFF file, ignore RGBA
			is >> r >> g >> b >> a;
	}

	index_t he = 0;
	for(index_t i = 0; i < n_faces_; i++)
	{
		is>>v;
		if(!i && v > che::P)
		{
			vertex * tGT = GT; GT = NULL;

			delete_me();
			init(n_v, n_f * (v - che::P + 1));

			GT = tGT;
		}

		for(index_t j = 0; j < v; j++)
			is>>VT[he++];

		// divide face
		if(v > che::P)
		{
			VT[he++] = VT[he - v - 1];
			VT[he++] = VT[he - v];

			i += (v - che::P);
		}
	}

	is.close();
}

void che_off::write_file(const che * mesh, const string & file)
{
	ofstream os(file + ".off");

	os << "OFF" << endl;
	os << mesh->n_vertices() << " " << mesh->n_faces() << " 0" << endl;

	for(size_t v = 0; v < mesh->n_vertices(); v++)
		os << mesh->gt(v) << endl;

	for(index_t he = 0; he < mesh->n_half_edges(); he++)
	{
		if(!(he % che::P)) os << che::P;
		os << " " << mesh->vt(he);
		if(he % che::P == che::P - 1) os << endl;
	}

	os.close();
}

