#include "mesh/che_off.h"

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_off::che_off(const string & file)
{
	init(file);
}

che_off::che_off(const che_off & mesh): che(mesh)
{
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
	
	if(soff[0] == 'N')
		VN = new vertex[n_vertices_];

	real_t alpha;	// color
	for(index_t i = 0; i < n_vertices_; i++)
	{
		is >> GT[i];
		if(soff[0] == 'C' || soff[1] == 'C')
			is >> VC[i] >> alpha;
		if(soff[0] == 'N')
			is >> VN[i];
	}
	
	if(soff[0] == 'C' || soff[1] == 'C')
	{
		#pragma omp parallel for
		for(index_t i = 0; i < n_vertices_; i++)
			VC[i] /= 255;
	}

	index_t he = 0;
	for(index_t i = 0; i < n_faces_; i++)
	{
		is >> v;
		if(!i && v > che::mtrig)
		{
			vertex * tGT = GT; GT = nullptr;

			delete_me();
			init(n_v, n_f * (v - che::mtrig + 1));

			GT = tGT;
		}

		for(index_t j = 0; j < v; j++)
			is >> VT[he++];

		// divide face
		if(v == che::mquad)
		{
			VT[he] = VT[he - v];		he++;
			VT[he] = VT[he - che::mtrig];	he++;

			i++;
		}
	}

	is.close();
}

void che_off::write_file(const che * mesh, const string & file, const che_off::type & off, const bool & pointcloud)
{
	ofstream os(file + ".off");

	os << off << endl;
	os << mesh->n_vertices() << " " << (pointcloud ? 0 : mesh->n_faces()) << " 0" << endl;
	
	for(size_t v = 0; v < mesh->n_vertices(); v++)
	{
		os << mesh->gt(v);

		if(off == NOFF) os << " " << mesh->normal(v);	// NOFF file

		os << endl;
	}
	
	if(!pointcloud)
		for(index_t he = 0; he < mesh->n_half_edges(); )
		{
			os << che::mtrig;
			for(index_t i = 0; i < che::mtrig; i++)
				os << " " << mesh->vt(he++);
			os << endl;
		}

	os.close();
}

ostream & operator << (ostream & os, const che_off::type & off)
{
	switch(off)
	{
		case che_off::OFF	: os << "OFF";		break;
		case che_off::NOFF	: os << "NOFF";		break;
		case che_off::COFF	: os << "COFF";		break;
		case che_off::CNOFF	: os << "CNOFF";	break;
	}

	return os;
}


} // namespace gproshan

