#include "geodesics.h"

/// Execute performance and accuracy test for ptp algorithm on cpu and gpu.
void main_test_geodesics_ptp(int nargs, const char ** args);

/// Exact geodesics computed using MeshLP https://github.com/areslp/matlab/tree/master/MeshLP/MeshLP,
/// Geodesics code: http://code.google.com/p/geodesic/
distance_t * load_exact_geodesics(const string & file, const size_t & n);

