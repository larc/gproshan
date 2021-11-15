#ifndef CHE_CUH
#define CHE_CUH

#include "mesh/vertex.cuh"
#include "mesh/che.h"


#define cu_for_star(he, mesh, v) for(index_t stop = mesh->EVT[v], he = mesh->EVT[v]; he != NIL; he = (he = mesh->OT[cu_prev(he)]) != stop ? he : NIL)


// geometry processing and shape analysis framework
namespace gproshan {


__host__ __device__
index_t cu_trig(index_t he);

__host__ __device__
index_t cu_next(index_t he);

__host__ __device__
index_t cu_prev(index_t he);

void cuda_create_CHE(CHE * h_che, CHE *& dd_che, CHE *& d_che, const bool & normal = false, const bool & color = false);

void cuda_free_CHE(CHE *& dd_che, CHE *& d_che);


} // namespace gproshan

#endif // CHE_CUH

