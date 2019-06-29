#include "geodesics.h"
#include "geodesics_ptp.h"

/// Execute performance and accuracy test for ptp algorithm on cpu and gpu.
void main_test_geodesics_ptp(const int & nargs, const char ** args);

double test_fast_marching(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test);

double test_ptp_gpu(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const toplesets_t & toplesets, const int & n_test);

double test_ptp_cpu(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const toplesets_t & toplesets, const int & n_test);

double test_heat_method_cholmod(distance_t & error, double & stime, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test);

double test_heat_method_gpu(distance_t & error, double & stime, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test);

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
vector<pair<index_t, distance_t> > iter_error_parallel_toplesets_propagation_coalescence_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, const distance_t * exact_dist, double & time_ptp);

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
vector<pair<index_t, distance_t> > iter_error_parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, const distance_t * exact_dist, double & time_ptp);

/// Return an array with the time per iteration.
double * times_farthest_point_sampling_ptp_gpu(che * mesh, vector<index_t> & samples, size_t n, distance_t radio = 0);

/// Return an array with the time per iteration.
double * times_farthest_point_sampling_ptp_coalescence_gpu(che * mesh, vector<index_t> & samples, size_t n, distance_t radio = 0);

/// Exact geodesics computed using MeshLP https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP,
/// Geodesics code: http://code.google.com/p/geodesic/
distance_t * load_exact_geodesics(const string & file, const size_t & n);

distance_t compute_error(const distance_t * dist, const distance_t * exact, const size_t & n, const size_t & s);

