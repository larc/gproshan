#ifndef D_MESH_APPS_H
#define D_MESH_APPS_H

#include "include.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

void dictionary_learning_process(che * mesh, vector<index_t> & points, const size_t & freq, size_t & rt, const size_t & m, size_t & M, const distance_t & f, const index_t & pf, const bool & op_dict = true);

void mesh_denoising(che * mesh, vector<index_t> & points, const size_t & freq, size_t & rt, const size_t & m, size_t & M, const distance_t & f, const bool & learn);

void mesh_super_resolution(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

void mesh_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

void mesh_iterative_inpaiting(che * mesh, vector<index_t> & points, size_t freq, size_t rt, size_t m, size_t M, double f, const bool & learn);

size_t sort_first_valid_vertex(index_t * indexes, const size_t & size, const vector<patches_map_t> & patches_map);

} // mdict

#endif // D_MESH_APPS_H

