#include "geodesics.h"

/// Execute performance and accuracy test for ptp algorithm on cpu and gpu.
void main_test_geodesics_ptp(const int & nargs, const char ** args);

/// Return an array with the error per iteration.
/// Starting to store (position 0) errors after number of toplesets.
distance_t * iter_error_parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, 
const distance_t * exact_dist, float & time_ptp);

/// Exact geodesics computed using MeshLP https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP,
/// Geodesics code: http://code.google.com/p/geodesic/
distance_t * load_exact_geodesics(const string & file, const size_t & n);

