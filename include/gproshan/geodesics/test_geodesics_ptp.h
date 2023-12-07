#ifndef TEST_GEODESICS_PTP_H
#define TEST_GEODESICS_PTP_H

#include <gproshan/geodesics/geodesics.h>
#include <gproshan/geodesics/geodesics_ptp.h>


// geometry processing and shape analysis framework
namespace gproshan {


/// Execute performance and accuracy test for ptp algorithm on cpu and gpu.
void main_test_geodesics_ptp(const int nargs, const char ** args);

double test_fast_marching(real_t & error, const real_t * exact, che * mesh, const std::vector<index_t> & source, const int n_test);

double test_ptp_cpu(real_t & error, const real_t * exact, che * mesh, const std::vector<index_t> & source, const toplesets_t & toplesets, const int n_test);

double test_heat_method_cholmod(real_t & error, double & stime, const real_t * exact, che * mesh, const std::vector<index_t> & source, const int n_test);


#ifdef GPROSHAN_CUDA

double test_ptp_gpu(real_t & error, const real_t * exact, che * mesh, const std::vector<index_t> & source, const toplesets_t & toplesets, const int n_test);

double test_heat_method_gpu(real_t & error, double & stime, const real_t * exact, che * mesh, const std::vector<index_t> & source, const int n_test);

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
std::vector<std::pair<index_t, real_t> > iter_error_parallel_toplesets_propagation_coalescence_gpu(che * mesh, const std::vector<index_t> & sources, const std::vector<index_t> & limits, const index_t * sorted_index, const real_t * exact_dist, double & time_ptp);

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
std::vector<std::pair<index_t, real_t> > iter_error_parallel_toplesets_propagation_gpu(che * mesh, const std::vector<index_t> & sources, const std::vector<index_t> & limits, const index_t * sorted_index, const real_t * exact_dist, double & time_ptp);

/// Return an array with the time per iteration.
double * times_farthest_point_sampling_ptp_gpu(che * mesh, std::vector<index_t> & samples, size_t n, real_t radio = 0);

/// Return an array with the time per iteration.
double * times_farthest_point_sampling_ptp_coalescence_gpu(che * mesh, std::vector<index_t> & samples, size_t n, real_t radio = 0);

#endif // GPROSHAN_CUDA


/// Exact geodesics computed using MeshLP https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP,
/// Geodesics code: http://code.google.com/p/geodesic/
real_t * load_exact_geodesics(const std::string & file, const size_t n);

real_t compute_error(const real_t * dist, const real_t * exact, const size_t n, const size_t s);


} // namespace gproshan

#endif // TEST_GEODESICS_PTP_H

