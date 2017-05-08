#include "decimation.h"

decimation::decimation()
{
	var_size = 10;
	for (int i = 0; i < var_size; i++)
		Q[i] = 0;
}


double compute_error(const  vertex_t & a, const  vertex_t & b)
{

}

decimation::decimation( vertex normal, vertex p)
{
	double nx, ny, nz, d;
	var_size = 10;
	nx = normal[0];
	ny = normal[1];
	nz = normal[2];
	d = -normal[0]* p[0] - normal[1]*p[1] - normal[2]*p[2];

	Q[0] = nx*nx;
	Q[1] = nx*ny;
	Q[2] = nx*nz;
	Q[3] = nx*d;
	Q[4] = nx*nx;
	Q[5] = ny*ny;
	Q[6] = ny*nz;
	Q[7] = ny*d;
	Q[8] = nz*nz;
	Q[9] = nz*d;
	Q[10] = d*d;

}

decimation::~decimation()
{
}


double decimation::eval(vertex v)
{
	double x, y, z;
	x = v[0];
	y = v[1];
	z = v[2];

	return (x*x*Q[0] + x*y*Q[1] + x*z*Q[2] + x*Q[3] +
			x*y*Q[1] + y*y*Q[5] + y*z*Q[6] + y*Q[7] +
			x*z*Q[2] + y*z*Q[6] + z*z*Q[8] + z*Q[9] +
			x*Q[3] + y*Q[7] + z*Q[9] + Q[10]);
}

