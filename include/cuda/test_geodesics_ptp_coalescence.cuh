#ifndef TEST_GEODESICS_PTP_COALESCENCE_CUH
#define TEST_GEODESICS_PTP_COALESCENCE_CUH

#include "che.cuh"


// geometry processing and shape analysis framework
namespace gproshan {


/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
std::vector<std::pair<index_t, real_t> > iter_error_run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const std::vector<index_t> & sources, const std::vector<index_t> & limits, const index_t * inv, const real_t * exact_dist, real_t * d_error);


} // namespace gproshan

#endif // TEST_GEODESICS_PTP_COALESCENCE_CUH

