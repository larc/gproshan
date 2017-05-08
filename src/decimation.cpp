#include "decimation.h"


double compute_error(const  vertex_t & a, const  vertex_t & b)
{

}

decimation::decimation( che * mesh_)
{
	mesh = mesh_;
	Q = new mat[mesh->n_vertices()];
}

void decimation::compute_quadrics()
{
	for(int i = 0; i < mesh->n_vertices(); i++)
	{
	/*	vec normal(4);
		normal[]
		= mesh->gt[i];*/
	}
}

decimation::~decimation()
{
}


double decimation::eval(vertex_t v)
{
	double x, y, z;
}

