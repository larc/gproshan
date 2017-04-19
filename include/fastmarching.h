#ifndef FASTMARCHING_H
#define FASTMARCHING_H

#include "include.h"
#include "off.h"

#include <cmath>
#include <vector>
#include <armadillo>

struct path
{
	face f;
	vertex v;
	size_t i;

	path(face fx_ = face(), vertex vx_ = vertex(), size_t vi_ = NIL);
};

class fastmarching
{
	public:
		distance_t * distances;
		index_t * cluster;

	private:
		path * paths;
		size_t * color;
		size_t n_vertices;
		bool spherical;
		distance_t radio;
		size_t * sort_index;
		size_t n_pesos;

	public:
		fastmarching(bool _spherical = 0);
		fastmarching(off & model, vector<index_t> & source_, distance_t r = INFINITY, bool normalized = 0, bool _spherical = 0);
		~fastmarching();
		size_t get_npesos();
		void print(ostream & os);
		void print_radio(ostream & os);
		distance_t operator[](size_t i); //distance's vertex i
		size_t operator()(size_t i); //distance's vertex i-th less distance
		void print_path(ostream & os);
		void print_path(off & model, size_t x, ostream & os, size_t source);
		void print_path(off & model, path p, ostream & os, distance_t d, size_t source);

	private:
		void calcular_FM(off & model, vector<index_t> & source);
		distance_t update(path p, path & r, off & model);
		distance_t planar_update(path & r, arma::mat & X, size_t ix0, size_t ix1, bool debug = 0);
		distance_t spherical_update(path & r, arma::mat & X, size_t ix0, size_t ix1);
		void normalize();
};

#endif // FASTMARCHING_H

