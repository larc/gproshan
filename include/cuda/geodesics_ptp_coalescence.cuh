#ifndef GEODESICS_PTP_COALESCENCE_CUH
#define GEODESICS_PTP_COALESCENCE_CUH

#include "che.cuh"
#include "geodesics_ptp.cuh"

index_t run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, distance_t * h_dist, distance_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * inv, distance_t * d_error, index_t * h_clusters = NULL, index_t ** d_clusters = NULL);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t end, index_t start = 0);

__global__
void relax_ptp_coalescence(CHE * mesh, distance_t * new_dist, distance_t * old_dist, index_t end, index_t start = 0);

__global__
void relative_error(distance_t * error, distance_t * new_dist, distance_t * old_dist, index_t n);

struct is_ok
{
	__host__ __device__
	bool operator()(const distance_t & val) const
	{
		return val < 1e-5;
	}
};

#endif // GEODESICS_PTP_COALESCENCE_CUH

