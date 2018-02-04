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
		index_t * paths;
		index_t * color;
		size_t n_vertices;
		bool spherical;
		bool normalized;
		distance_t radio;
		index_t * sort_index;
		size_t n_pesos;
		size_t n_iter;
		const vector<index_t> & sources;

	public:
		geodesics(che * mesh, const vector<index_t> & _sources, size_t _n_iter = 0,
				 distance_t _radio = INFINITY, bool normalized = 0, bool parallel = 0, bool _spherical = 0);
		virtual ~geodesics();
		index_t operator[](index_t i);
		index_t farthest();
		size_t get_n_radio();
		void get_sort_indexes(index_t * indexes, size_t K);
		void execute(che * mesh);
		void path_to(vector<vertex> & v_path, che * mesh, const index_t & v);
		void path_to(vector<vertex> & v_path, che * mesh, const vertex & v, const index_t & he);
		void normalize();

	private:
		void run(che * mesh);
		distance_t update(index_t & d, che * mesh, const index_t & he, vertex & vx);
		distance_t planar_update(index_t & d, arma::mat & X, index_t * x, vertex & vx);
};

distance_t update_step(che * mesh, const distance_t * dist, const index_t & he);
distance_t * fast_geodesics(che * mesh, index_t * source, length_t source_size, const vector<index_t> & limites, const index_t * sorted);

#endif //GEODESICS_H

