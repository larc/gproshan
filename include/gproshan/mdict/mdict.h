#ifndef MDICT_H
#define MDICT_H

#include <gproshan/include.h>
#include <gproshan/mdict/patch.h>
#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


// SPARSE

struct locval_t
{
	arma::uword i, j;
	real_t val;
};


void OMP(std::vector<locval_t> & alpha, const arma::fvec & x, const index_t i, const arma::fmat & D, const size_t L);

arma::sp_fmat OMP_all(std::vector<locval_t> & locval, const arma::fmat & X, const arma::fmat & D, const size_t L);

void sp_KSVD(arma::fmat & D, const arma::fmat & X, const size_t L, size_t k);


// DENSE

std::tuple<arma::fvec, arma::uvec> _OMP(const arma::fvec & x, const arma::fmat & D, const size_t L);

arma::fvec OMP(const arma::fvec & x, const arma::fmat & D, const size_t L);
arma::fvec OMP(const arma::fvec & x, const arma::fmat & D, const size_t L, const arma::uchar_vec & mask);

arma::fmat OMP_all(const arma::fmat & X, const arma::fmat & D, const size_t L);

void KSVD(arma::fmat & D, const arma::fmat & X, const size_t L, size_t k);


// MESH DENSE

arma::fvec OMP(const patch & p, const arma::fmat & A, const size_t L);
arma::fvec OMP(const patch & p, const arma::fmat & A, const size_t L);

arma::fmat OMP_all(const std::vector<patch> & patches, const arma::fmat & A, const size_t L);
arma::fmat OMP_all(const std::vector<patch> & patches, basis * phi_basis, const arma::fmat & A, const size_t L);


void KSVD(arma::fmat & A, const std::vector<patch> & patches, const size_t L, size_t k);


// MESH SPARSE

void OMP(std::vector<locval_t> & alpha, const patch & p, const index_t i, const arma::fmat & A, const size_t L);

arma::sp_fmat OMP_all(std::vector<locval_t> & locval, const std::vector<patch> & patches, const arma::fmat & A, const size_t L);

void sp_KSVD(arma::fmat & A, const std::vector<patch> & patches, const size_t L, size_t k);


} // namespace gproshan::mdict

#endif // MDICT_H

