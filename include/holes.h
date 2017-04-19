#ifndef HOLES_H
#define HOLES_H

#include <vector>
#include <map>
#include <Eigen/Eigen>

#include "keypoints.h"
#include "off.h"

using namespace Eigen;

typedef double curv_t;

class holes
{
	private:
		off shape_holes;
		size_t n_vertices;
		size_t n_faces;
		off * mesh_holes;
		pair<curv_t, curv_t> * curv_holes; //[SUM +][SUM -]
		curv_t * tangentes;
		bool * holes_vertices;
		bool * holes_faces;
		bool * deleted_faces;
		size_t * neighbors_bound;
		size_t * mapea;
		vector<size_t> index_holes;

	public:
		holes(keypoints & kps, off & shape);
		~holes();
		void save_off(string file);
		void print_holes(ostream & os);
		size_t getNroHoles();

	private:
		void calcular_holes(keypoints & kps, off & shape);
		void calcular_bounds(off & shape);
		void find_holes();
		void repair_hole(off & shape, size_t index);
		void add_mesh(size_t index, vector<vertex> & vertices, vector<face> & faces);
		void divide_faces(size_t k, size_t nf, vector<vertex> & vertices, vector<face> & faces);
		void divide_faces(size_t k, vector<vertex> & vertices, map<vertex, size_t> & m_vertices, vector<face> & faces);
		void divide_faces(size_t k, size_t nf, vector<vertex> & vertices, map<vertex, size_t> & m_vertices, vector<face> & faces);
		void biharmonic_interp_2(MatrixXd & P, MatrixXd & H);
		void analysis_curvature_2(size_t index, MatrixXd & P, MatrixXd & H, off & p_mesh, off & h_mesh);
		curv_t analysis_curvature_2(MatrixXd & alpha, MatrixXd & P, MatrixXd & H, off & mesh);
		MatrixXd pca(MatrixXd & p, MatrixXd & m);
		void sub_mesh(off & shape, vector<vertex> & vertices, vector<face> & faces, size_t k = 1);
		off * sub_mesh(off & shape, size_t k = 1);
};

#endif // HOLES_H
