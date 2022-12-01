#ifndef GEODESICS_PTP_COALESCENCE_CUH
#define GEODESICS_PTP_COALESCENCE_CUH

#include <gproshan/mesh/che.cuh>
#include <gproshan/geodesics/geodesics_ptp.cuh>
#include <gproshan/geodesics/geodesics_ptp.h>


// geometry processing and shape analysis framework
namespace gproshan {


index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const std::vector<index_t> & sources, const toplesets_t & inv, real_t * d_error, index_t * h_clusters = nullptr, index_t ** d_clusters = nullptr);


} // namespace gproshan

#endif // GEODESICS_PTP_COALESCENCE_CUH

