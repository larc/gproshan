#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include "che.h"

index_t iterations(const vector<index_t> &);
index_t start_v(const index_t & i, const vector<index_t> & limits);
index_t end_v(const index_t & i, const vector<index_t> & limits);

distance_t * parallel_toplesets_propagation_coalescence_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, double & time_ptp, index_t * clusters = NULL);

distance_t * parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, double & time_ptp, index_t * clusters = NULL);

distance_t * parallel_toplesets_propagation_cpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, index_t * clusters = NULL);

distance_t farthest_point_sampling_ptp_gpu(che * mesh, vector<index_t> & samples, double & time_fps, size_t n, distance_t radio = 0);

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he);

void normalize_ptp(distance_t * dist, const size_t & n);

#endif // GEODESICS_PTP_H

