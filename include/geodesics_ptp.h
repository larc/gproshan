#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include "che.h"

#define PTP_TOL 1e-3

index_t iterations(const vector<index_t> &);
index_t start_v(const index_t & i, const vector<index_t> & limits);
index_t end_v(const index_t & i, const vector<index_t> & limits);

double parallel_toplesets_propagation_coalescence_gpu(distance_t * dist, che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, index_t * clusters = NULL);

double parallel_toplesets_propagation_gpu(distance_t * dist, che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, index_t * clusters = NULL);

void parallel_toplesets_propagation_cpu(distance_t *& dist, che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, index_t * clusters = NULL);

distance_t farthest_point_sampling_ptp_gpu(che * mesh, vector<index_t> & samples, double & time_fps, size_t n, distance_t radio = 0);

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he);

void normalize_ptp(distance_t * dist, const size_t & n);

#endif // GEODESICS_PTP_H

