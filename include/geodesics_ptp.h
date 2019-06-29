#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include "che.h"

#define PTP_TOL 1e-3


struct ptp_out_t
{
	distance_t * dist;
	index_t * clusters;
	
	ptp_out_t(distance_t *const & d, index_t *const & c = NULL);
};

struct toplesets_t
{
	const vector<index_t> & limits;
	const index_t *const & index;
};

che * ptp_coalescence(index_t * & inv, che * mesh, const toplesets_t & toplesets);

double parallel_toplesets_propagation_coalescence_gpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets, const bool & set_inf = 1);

double parallel_toplesets_propagation_gpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_coalescence_cpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const vector<index_t> & sources, const toplesets_t & toplesets);

distance_t farthest_point_sampling_ptp_gpu(che * mesh, vector<index_t> & samples, double & time_fps, size_t n, distance_t radio = 0);

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he);

void normalize_ptp(distance_t * dist, const size_t & n);

#endif // GEODESICS_PTP_H

