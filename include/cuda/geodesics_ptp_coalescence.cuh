#ifndef GEODESICS_PTP_COALESCENCE_CUH
#define GEODESICS_PTP_COALESCENCE_CUH

#include "che.cuh"
#include "geodesics_ptp.cuh"
#include "geodesics_ptp.h"


// geometry processing and shape analysis framework
namespace gproshan {


index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const std::vector<index_t> & sources, const toplesets_t & inv, distance_t * d_error, index_t * h_clusters = nullptr, index_t ** d_clusters = nullptr);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t end, index_t start = 0);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t end, index_t start = 0);


} // namespace gproshan

#endif // GEODESICS_PTP_COALESCENCE_CUH

