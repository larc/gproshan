#ifndef GEODESICS_H
#define GEODESICS_H

#include <gproshan/mesh/che.h>

#include <cmath>
#include <functional>


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
	using fm_function_t = std::function<bool (const index_t)>;

	public:
		enum algorithm {	FM,				///< Execute Fast Marching algorithm
							PTP_CPU,		///< Execute Parallel Toplesets Propagation CPU algorithm
							HEAT_METHOD,	///< Execute Heat Method - cholmod (CPU)
						#ifdef GPROSHAN_CUDA
							PTP_GPU,		///< Execute Parallel Toplesets Propagation GPU algorithm
							HEAT_METHOD_GPU	///< Execute Heat Method - cusparse (GPU)
						#endif // GPROSHAN_CUDA
						};

		struct params
		{
			algorithm alg		= FM;					///< specific the algorithm to execute.
			size_t n_iter		= 0;					///< maximum number of iterations.
			float radio		= INFINITY;				///< execute until the specific radio.
			float * dist_alloc	= nullptr;				///< external dist allocation
			bool cluster		= false;				///< to cluster vertices to closest source.
			fm_function_t fun	= nullptr;				///< fun is executed inside FM loop
		};

	public:
		index_t * clusters;			///< Clustering vertices to closest source.

	private:
		float * dist;				///< Results of computation geodesic distances.
		index_t * sorted_index;		///< Sort vertices by topological level or geodesic distance.
		size_t n_sorted;			///< Number of vertices sorted by their geodesics distance.
		bool free_dist;

		const size_t n_vertices;	///< Number of vertices, const reference

	public:
		geodesics(	che * mesh,								///< input triangular mesh.
					const std::vector<index_t> & sources,	///< source vertices.
					const params & p = {FM, 0, INFINITY, nullptr, false, nullptr}
					);

		virtual ~geodesics();
		operator const float * () const;
		float operator[](const index_t i) const;
		index_t operator()(const index_t i) const;
		float radio() const;
		index_t farthest() const;
		size_t n_sorted_index() const;
		void copy_sorted_index(index_t * indexes, const size_t n) const;
		void normalize();

	private:
		void execute(che * mesh, const std::vector<index_t> & sources, const params & p);
		void run_fastmarching(che * mesh, const std::vector<index_t> & sources, const size_t n_iter, const float radio, const fm_function_t & fun);
		void run_parallel_toplesets_propagation_cpu(che * mesh, const std::vector<index_t> & sources);
		void run_heat_method(che * mesh, const std::vector<index_t> & sources);

#ifdef GPROSHAN_CUDA
		void run_parallel_toplesets_propagation_gpu(che * mesh, const std::vector<index_t> & sources);
		void run_heat_method_gpu(che * mesh, const std::vector<index_t> & sources);
#endif // GPROSHAN_CUDA

};


} // namespace gproshan

#endif //GEODESICS_H

