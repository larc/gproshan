#ifndef GEODESICS_H
#define GEODESICS_H

#include "include.h"
#include "che.h"

#include <armadillo>

class geodesics
{
	public:
		distance_t * distances;
		index_t * clusters;

	private:
		index_t * color;
		index_t * sorted_index;
		size_t n_vertices;
		size_t n_sorted;
		size_t n_iter;

	public:
		geodesics(che * mesh, const vector<index_t> & _sources, const size_t & _n_iter = 0, const distance_t & radio = INFINITY, const bool & parallel = 0);
		virtual ~geodesics();
		const index_t & operator[](const index_t & i) const;
		const index_t & farthest() const;
		const size_t & n_sorted_index() const;
		void copy_sorted_index(index_t * indexes, const size_t & n) const;
		void normalize();

	private:
		void execute(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio, const bool & parallel);
		void run_fastmarching(che * mesh, const vector<index_t> & sources, const size_t & n_iter, const distance_t & radio);
		
		distance_t update(index_t & d, che * mesh, const index_t & he, vertex & vx);
		distance_t planar_update(index_t & d, arma::mat & X, index_t * x, vertex & vx);
};

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he);

#endif //GEODESICS_H

