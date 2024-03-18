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

#endif // GPROSHAN_CUDA


/// Exact geodesics computed using MeshLP https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP,
/// Geodesics code: http://code.google.com/p/geodesic/
real_t * load_exact_geodesics(const std::string & file, const size_t n);

real_t compute_error(const real_t * dist, const real_t * exact, const size_t n, const size_t s);


} // namespace gproshan

#endif // TEST_GEODESICS_PTP_H

