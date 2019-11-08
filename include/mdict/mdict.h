#ifndef MDICT_H
#define MDICT_H

#include "include.h"
#include "patch.h"
#include "d_mesh.h"

#include "include_arma.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


// SPARSE

struct locval_t
{
	arma::uword i, j;
	real_t val;

	bool operator < (const locval_t & lc);
};

void OMP(vector<locval_t> & alpha, const a_vec & x, const index_t & i, const a_mat & D, const size_t & L);

a_sp_mat OMP_all(vector<locval_t> & locval, const a_mat & X, const a_mat & D, const size_t & L);

void sp_KSVD(a_mat & D, const a_mat & X, const size_t & L, size_t k);


// DENSE

tuple<a_vec, arma::uvec> _OMP(const a_vec & x, const a_mat & D, const size_t & L);

a_vec OMP(const a_vec & x, const a_mat & D, const size_t & L);

a_mat OMP_all(const a_mat & X, const a_mat & D, const size_t & L);

void KSVD(a_mat & D, const a_mat & X, const size_t & L, size_t k);


// MESH DENSE

a_vec OMP(const a_mat & A, const patch & p, const size_t & L);

a_mat OMP_all(const vector<patch> & patches, const a_mat & A, const size_t & L);

void KSVD(a_mat & A, const vector<patch> & patches, const size_t & L, size_t k);


} // namespace gproshan::mdict

#endif // MDICT_H

