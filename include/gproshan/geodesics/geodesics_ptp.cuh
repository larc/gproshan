#ifndef GEODESICS_PTP_CUH
#define GEODESICS_PTP_CUH

#include <gproshan/mesh/che.cuh>


#define NT 64
#define NB(x) (x + NT - 1) / NT


// geometry processing and shape analysis framework
namespace gproshan {


index_t run_ptp_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const std::vector<index_t> & sources, const std::vector<index_t> & limits, const index_t * h_sorted, index_t * d_sorted, real_t * d_error, index_t * h_clusters = nullptr, index_t ** d_clusters = nullptr);

__forceinline__ __device__
real_t cu_update_step(CHE * mesh, const real_t * dist, const index_t & he);

__global__
void relax_ptp(CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, index_t * sorted, index_t end, index_t start = 0);

__global__
void relax_ptp(CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * sorted, index_t end, index_t start = 0);

__global__
void relative_error(real_t * error, const real_t * new_dist, const real_t * old_dist, const index_t start, const index_t end, const index_t * sorted = nullptr);


struct is_ok
{
	__host__ __device__
	bool operator()(const real_t & val) const;
};


} // namespace gproshan

#endif // GEODESICS_PTP_CUH

