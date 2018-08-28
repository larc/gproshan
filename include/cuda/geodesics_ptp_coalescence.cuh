#ifndef GEODESICS_PTP_COALESCENCE_CUH
#define GEODESICS_PTP_COALESCENCE_CUH

#include "che.cuh"

#define NT 32
#define NB(x) (x + NT - 1) / NT

index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * inv, index_t * h_clusters = NULL, index_t ** d_clusters = NULL);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t end, index_t start = 0);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t end, index_t start = 0);

#endif // GEODESICS_PTP_COALESCENCE_CUH

