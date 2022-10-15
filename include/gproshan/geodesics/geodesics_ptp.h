#ifndef GEODESICS_PTP_H
#define GEODESICS_PTP_H

/**
	Parallel Toplesets Propagation algorithm
	topleset: (Topological Level Set)
*/

#include <gproshan/mesh/che.h>

#define PTP_TOL 1e-4


// geometry processing and shape analysis framework
namespace gproshan {


struct ptp_out_t
{
	real_t * dist;
	index_t * clusters;

	ptp_out_t(real_t *const & d, index_t *const & c = nullptr);
};

struct toplesets_t
{
	const std::vector<index_t> & limits;
	const index_t *const & index;
};

che * ptp_coalescence(index_t * & inv, const che * mesh, const toplesets_t & toplesets);

double parallel_toplesets_propagation_coalescence_gpu(const ptp_out_t & ptp_out, const che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets, const bool & set_inf = 1);

double parallel_toplesets_propagation_gpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_coalescence_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

void parallel_toplesets_propagation_cpu(const ptp_out_t & ptp_out, che * mesh, const std::vector<index_t> & sources, const toplesets_t & toplesets);

real_t farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, double & time_fps, size_t n, real_t radio = 0);

real_t update_step(che * mesh, const real_t * dist, const index_t & he);

void normalize_ptp(real_t * dist, const size_t & n);


} // namespace gproshan

#endif // GEODESICS_PTP_H

