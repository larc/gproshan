#include "mesh/che_off.h"

#include <cstring>
#include <cassert>
#include <cstdio>
#include <fstream>

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
	char soff[32];
	size_t nv, nf, ne;
	
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	fgets(soff, sizeof(soff), fp);
	fscanf(fp, "%lu %lu %lu", &nv, &nf, &ne);
	
	alloc(nv, nf);
	
	float x, y, z, r, g, b, a;
	for(index_t v = 0; v < n_vertices; ++v)
	{
		fscanf(fp, "%f %f %f", &x, &y, &z);
		GT[v] = { x, y, z };

		if(soff[0] == 'C' || soff[1] == 'C')
		{
			fscanf(fp, "%f %f %f %f", &r, &g, &b, &a);
			VC[v] = { r, g, b };
		}

		if(soff[0] == 'N')
		{
			fscanf(fp, "%f %f %f", &x, &y, &z);
			VN[v] = { x, y, z };
		}
	}
	
	if(soff[0] == 'C' || soff[1] == 'C')
	{
		#pragma omp parallel for
		for(index_t i = 0; i < n_vertices; ++i)
			VC[i] /= 255;
	}

	index_t he = 0;
	for(index_t i = 0; i < n_faces; ++i)
	{
		fscanf(fp, "%lu", &ne);
		if(!i && ne > che::mtrig)
		{
			vertex * tGT = GT; GT = nullptr;

			free();
			alloc(nv, nf * (ne - che::mtrig + 1));

			GT = tGT;
		}

		for(index_t j = 0; j < ne; ++j)
			fscanf(fp, "%u", VT + he++);

		// divide face
		if(ne == che::mquad)
		{
			VT[he] = VT[he - ne];			++he;
			VT[he] = VT[he - che::mtrig];	++he;

			++i;
		}
	}

	fclose(fp);
}

void che_off::write_file(const che * mesh, const string & file, const che_off::type & off, const bool & pointcloud)
{
	ofstream os(file + ".off");

	os << off << endl;
	os << mesh->n_vertices << " " << (pointcloud ? 0 : mesh->n_faces) << " 0" << endl;
	
	for(size_t v = 0; v < mesh->n_vertices; ++v)
	{
		os << mesh->gt(v);

		if(off == NOFF) os << " " << mesh->normal(v);	// NOFF file

		os << endl;
	}
	
	if(!pointcloud)
		for(index_t he = 0; he < mesh->n_half_edges; )
		{
			os << che::mtrig;
			for(index_t i = 0; i < che::mtrig; ++i)
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

