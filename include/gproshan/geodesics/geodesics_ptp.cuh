#ifndef GEODESICS_PTP_CUH
#define GEODESICS_PTP_CUH

#include <gproshan/mesh/che.cuh>
#include <gproshan/geodesics/geodesics_ptp.h>


#define NT 64
#define NB(x) (x + NT - 1) / NT


// geometry processing and shape analysis framework
namespace gproshan {


index_t run_ptp_gpu(const CHE * d_mesh, const std::vector<index_t> & sources, const index_t & n_vertices,
					real_t * h_dist, real_t ** d_dist, const toplesets_t & inv, real_t * d_error,
					index_t * h_clusters = nullptr, index_t ** d_clusters = nullptr, index_t * d_sorted = nullptr);

__global__
void relax_ptp(const CHE * mesh, real_t * new_dist, real_t * old_dist, index_t * new_clusters, index_t * old_clusters, const index_t start, const index_t end, const index_t * sorted = nullptr);

__global__
void relative_error(real_t * error, const real_t * new_dist, const real_t * old_dist, const index_t start, const index_t end, const index_t * sorted = nullptr);


struct is_ok
{
	__host__ __device__
	bool operator()(const real_t & val) const;
};


} // namespace gproshan

#endif // GEODESICS_PTP_CUH

