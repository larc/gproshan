#ifndef GEODESICS_H
#define GEODESICS_H

#include "che.h"
#include "include_arma.h"


// geometry processing and shape analysis framework
namespace gproshan {


/*!
	Compute the geodesics distances on a mesh from a source or multi-source. This class implements
	the Fast Marching algorithm without deal with obtuse triangles. Also, if the options PTP_CPU or
	PTP_GPU are enables, compute the geodesics distances executing the Parallel Toplesets Propagation
	algorithm.
*/
class geodesics
{
	public:
		enum algorithm {	FM,				///< Execute Fast Marching algorithm
						#ifdef GPROSHAN_CUDA
							PTP_GPU,		///< Execute Parallel Toplesets Propagation GPU algorithm
							HEAT_FLOW_GPU,	///< Execute Heat Flow - cusparse (GPU)
						#endif // GPROSHAN_CUDA
							PTP_CPU,		///< Execute Parallel Toplesets Propagation CPU algorithm
							HEAT_FLOW		///< Execute Heat Flow - cholmod (CPU)
						};
		
		struct params
		{
			algorithm alg		= FM;					///< specific the algorithm to execute.
			size_t n_iter		= 0;					///< maximum number of iterations.
			real_t radio		= INFINITY;				///< execute until the specific radio.
			real_t * dist_alloc	= nullptr;				///< external dist allocation
			bool cluster		= false;				///< to cluster vertices to closest source.
			bool (*fun) (const index_t & u) = nullptr;	///< fun is executed inside FM loop
		};

	public:
		index_t * clusters;			///< Clustering vertices to closest source.

	private:
		real_t * dist;				///< Results of computation geodesic distances.
		index_t * sorted_index;		///< Sort vertices by topological level or geodesic distance.
		size_t n_sorted;			///< Number of vertices sorted by their geodesics distance.
		bool free_dist;
		
		const size_t & n_vertices;	///< Number of vertices, const reference

	public:
		geodesics(	che * mesh,								///< input triangular mesh.
					const std::vector<index_t> & sources,	///< source vertices.
					const params & p = {FM, 0, INFINITY, nullptr, false, nullptr}
					);

		virtual ~geodesics();
		const real_t & operator[](const index_t & i) const;
		const index_t & operator()(const index_t & i) const;
		const real_t & radio() const;
		const index_t & farthest() const;
		const size_t & n_sorted_index() const;
		void copy_sorted_index(index_t * indexes, const size_t & n) const;
		void normalize();

	private:
		void execute(che * mesh, const std::vector<index_t> & sources, const size_t & n_iter, const real_t & radio, const algorithm & alg);
		void run_fastmarching(che * mesh, const std::vector<index_t> & sources, const size_t & n_iter, const real_t & radio);
		void run_parallel_toplesets_propagation_cpu(che * mesh, const std::vector<index_t> & sources, const size_t & n_iter, const real_t & radio);
		void run_heat_flow(che * mesh, const std::vector<index_t> & sources);
		
#ifdef GPROSHAN_CUDA
		void run_parallel_toplesets_propagation_gpu(che * mesh, const std::vector<index_t> & sources, const size_t & n_iter, const real_t & radio);
		void run_heat_flow_gpu(che * mesh, const std::vector<index_t> & sources);
#endif // GPROSHAN_CUDA

		real_t update(index_t & d, che * mesh, const index_t & he, vertex & vx);
		real_t planar_update(index_t & d, a_mat & X, index_t * x, vertex & vx);
};


} // namespace gproshan

#endif //GEODESICS_H

