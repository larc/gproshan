#ifndef MDICT_H
#define MDICT_H

#include "include.h"
#include "patch.h"
#include "d_mesh.h"

#include "include_arma.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


a_vec OMP(const a_vec & x, const a_mat & D, const size_t & L);

void KSVD(a_mat & D, const a_mat & X, const size_t & L, size_t k);

void OMP_patch(a_mat & alpha, const a_mat & A, const index_t & i, patch & p, const size_t & L);

void OMP_all_patches_ksvt(a_mat & alpha, a_mat & A, std::vector<patch> & patches, size_t M, size_t L);

void KSVDT(a_mat & A, std::vector<patch> & patches, size_t M, size_t L);

[[deprecated]]
void OMP_patch(a_mat & alpha, const a_mat & A, const index_t & i, patch_t & p, const size_t & L);

[[deprecated]]
void OMP_all_patches_ksvt(a_mat & alpha, a_mat & A, std::vector<patch_t> & patches, size_t M, size_t L);

[[deprecated]]
void KSVDT(a_mat & A, std::vector<patch_t> & patches, size_t M, size_t L);


} // namespace gproshan::mdict

#endif // MDICT_H

