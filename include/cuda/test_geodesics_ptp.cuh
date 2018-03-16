#ifndef TEST_GEODESICS_PTP_CUH
#define TEST_GEODESICS_PTP_CUH

#include "che.cuh"

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
distance_t * iter_error_run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, const distance_t * exact_dist);

#endif // TEST_GEODESICS_PTP_CUH

