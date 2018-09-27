#include "vertex.cuh"

__host__ __device__
vertex_cu operator*(const real_t & a, const vertex_cu & v)
{
	return vertex_cu(a * v.x, a * v.y, a * v.z);
}

