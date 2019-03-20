#ifndef GEODESICS_PTP_CUH
#define GEODESICS_PTP_CUH

#include "che.cuh"

#define NT 64
#define NB(x) (x + NT - 1) / NT

index_t run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, distance_t * d_error, index_t * h_clusters = NULL, index_t ** d_clusters = NULL);

__forceinline__ __device__
distance_t cu_update_step(CHE * mesh, const distance_t * dist, const index_t & he);

__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t * sorted, index_t end, index_t start = 0);

__global__
void relax_ptp(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * sorted, index_t end, index_t start = 0);

__global__
void relative_error(distance_t * error, const distance_t * new_dist, const distance_t * old_dist, const index_t start, const index_t end, const index_t * sorted = NULL);


struct is_ok
{
	__host__ __device__
	bool operator()(const distance_t & val) const;
};

#endif // GEODESICS_PTP_CUH

