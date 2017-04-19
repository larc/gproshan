#ifndef D_OMP_H
#define D_OMP_H

#include "include.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

void OMP(vec & alpha, vec & x, mat & D, size_t L);

void KSVD(mat & D, mat & X, size_t L);

void OMP_patch(mat & alpha, const mat & A, const index_t & i, patch & p, const size_t & L);

void OMP_all_patches_ksvt(mat & alpha, mat & A, vector<patch> & patches, size_t M, size_t L);

void KSVDT(mat & A, vector<patch> & patches, size_t M, size_t L);

#endif

