#ifndef D_MESH_APPS_H
#define D_MESH_APPS_H

#include "include.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

void dictionary_learning_process(che * mesh, vector<index_t> & points, const size_t & freq,  size_t & rt, const size_t & m, size_t & M, const distance_t & f, const index_t & pf, const bool & op_dict = true);

void mesh_denoising(che * mesh, vector<index_t> & points, const size_t & freq, size_t & rt, const size_t & m, size_t & M, const distance_t & f, const  bool & learn);

void mesh_super_resolution(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

void mesh_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

void mesh_iterative_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

size_t sort_first_valid_vertex(index_t * indexes, const size_t & size, const vector<patches_map_t> & patches_map);

void plot_atoms(phi_function_t phi, params_t params, const distance_t & radio, const mat & A, string file = "atoms.gpi");

void plot_phi(phi_function_t phi, params_t params, distance_t radio, string file = "phi.gpi");

void plot_phi(phi_function_t phi, params_t params, const vertex_t & radio, const size_t & freq, const size_t & rt, string file = "phi.gpi");

#endif

