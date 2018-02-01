#ifndef D_DICT_LEARNING_H
#define D_DICT_LEARNING_H

#include "include.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

void OMP(vec & alpha, vec & x, mat & D, size_t L);

void KSVD(mat & D, mat & X, size_t L);

void OMP_patch(mat & alpha, const mat & A, const index_t & i, patch_t & p, const size_t & L);

void OMP_all_patches_ksvt(mat & alpha, mat & A, vector<patch_t> & patches, size_t M, size_t L);

void KSVDT(mat & A, vector<patch_t> & patches, size_t M, size_t L);

} // mdict

#endif // D_DICT_LEARNING_H

