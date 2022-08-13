#ifndef CHE_CUH
#define CHE_CUH

#include <gproshan/mesh/che.h>


#define cu_for_star(he, mesh, v) for(index_t stop = mesh->EVT[v], he = mesh->EVT[v]; he != NIL; he = (he = mesh->OT[cu_prev(he)]) != stop ? he : NIL)


// geometry processing and shape analysis framework
namespace gproshan {


__host__ __device__
index_t cu_trig(index_t he);

__host__ __device__
index_t cu_next(index_t he);

__host__ __device__
index_t cu_prev(index_t he);


struct CHE
{
	size_t n_vertices = 0;
	size_t n_faces = 0;
	size_t n_half_edges = 0;

	vertex * GT	= nullptr;
	vertex * VN	= nullptr;
	che::rgb_t * VC	= nullptr;
	index_t * VT	= nullptr;
	index_t * OT	= nullptr;
	index_t * EVT	= nullptr;

	CHE() = default;
	CHE(const che * mesh);
};


void cuda_create_CHE(CHE * h_che, CHE *& dd_che, CHE *& d_che, const bool & normal = false, const bool & color = false);

void cuda_free_CHE(CHE *& dd_che, CHE *& d_che);


} // namespace gproshan

#endif // CHE_CUH

