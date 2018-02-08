#ifndef GEODESICS_H
#define GEODESICS_H

#include "include.h"
#include "che.h"

#include <armadillo>

class geodesics
{
	public:
		enum option_t {FM, PTP_CPU, PTP_GPU};

	public:
		distance_t * distances;
		index_t * clusters;

	private:
		index_t * sorted_index;
		size_t n_vertices;
		size_t n_sorted;

	public:
		geodesics(che * mesh, const vector<index_t> & _sources, const option_t & opt = FM, const size_t & n_iter = 0, const distance_t & radio = INFINITY);
		virtual ~geodesics();
		const index_t & operator[](const index_t & i) const;
		const index_t & farthest() const;
		const size_t & n_sorted_index() const;
		void copy_sorted_index(index_t * indexes, const size_t & n) const;
		void normalize();

	private:
		void execute(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio, const option_t & opt);
		void run_fastmarching(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio);
		void run_parallel_toplesets_propagation_cpu(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio);
		void run_parallel_toplesets_propagation_gpu(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio);

		distance_t update(index_t & d, che * mesh, const index_t & he, vertex & vx);
		distance_t planar_update(index_t & d, arma::mat & X, index_t * x, vertex & vx);
};

#endif //GEODESICS_H

