#include "d_dict_learning.h"
#include "sampling.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


a_vec OMP(const a_vec & x, const a_mat & D, const size_t & L)
{
	a_vec alpha(D.n_cols, arma::fill::zeros);
	arma::uvec selected_atoms(L, arma::fill::zeros);

	real_t sigma = 0.001;
	real_t threshold = norm(x) * sigma;

	a_mat DD;
	a_vec aa, r = x;
	
	index_t l = 0;
	while(norm(r) > threshold && l < L)
	{
		selected_atoms(l) = index_max(abs(D.t() * r));
		
		DD = D.cols(selected_atoms.head(l + 1));
		aa = pinv(DD) * x;
		r = x - DD * aa;

		l++;
	}

	alpha.elem(selected_atoms.head(l)) = aa;
	return alpha;
}

void KSVD(a_mat & D, a_mat & X, size_t L)
{
//	size_t n = X.n_rows;
	size_t m = D.n_cols;
	size_t M = X.n_cols;

	a_mat alpha(m, M);

	size_t iter = L;
	while(iter--)
	{
		#pragma omp parallel for
		for(index_t i = 0; i < M; i++)
			alpha.col(i) = OMP(X.col(i), D, L);

		#pragma omp parallel for
		for(index_t j = 0; j < m; j++)
		{
			arma::uvec omega = find(abs(alpha.row(j)) > 0);
			if(omega.n_elem)
			{
				a_mat E = X - D * alpha - D.col(j) * alpha.row(j);
				E = E.cols(omega);
				a_mat U;
				a_vec s;
				a_mat V;
				svd(U, s, V, E);
				D.col(j) = U.col(0);
				a_rowvec a = alpha.row(j);
				a.elem(omega) = s(0) * V.col(0);
				alpha.row(j) = a;
			}
			//else
			//	D.col(j).randu();
		}
	}
}

void OMP_patch(a_mat & alpha, const a_mat & A, const index_t & i, patch & p, const size_t & L)
{
	alpha.col(i) = OMP(p.xyz.row(2).t(), p.phi * A, L);
}

void OMP_all_patches_ksvt(a_mat & alpha, a_mat & A, vector<patch> & patches, size_t M, size_t L)
{
	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
		OMP_patch(alpha, A, i, patches[i], L);
}

void KSVDT(a_mat & A, vector<patch> & patches, size_t M, size_t L)
{
	size_t K = A.n_rows;
	size_t m = A.n_cols;

	a_mat alpha(m, M);

	size_t iter = L;
	while(iter--)
	{
		OMP_all_patches_ksvt(alpha, A, patches, M, L);

		#pragma omp parallel for
		for(index_t j = 0; j < m; j++)
		{
			arma::uvec omega = find(abs(alpha.row(j)) > 0);

			a_mat sum(K, K, arma::fill::zeros);
			a_vec sum_error(K, arma::fill::zeros);

			for(arma::uword o: omega)
			{
				sum += alpha(j, o) * patches[o].phi.t() * patches[o].phi;

				a_mat D = patches[o].phi * A;
				D.col(j).zeros();
				a_mat e = patches[o].xyz.row(2).t() - D * alpha.col(o);

				sum_error += alpha(j, o) * patches[o].phi.t() * e;
			}
			if(omega.size())
			{
				a_vec X;
				solve(X, sum, sum_error);
				A.col(j) = X;
			}
		}
	}
}

void OMP_patch(a_mat & alpha, const a_mat & A, const index_t & i, patch_t & p, const size_t & L)
{
	alpha.col(i) = OMP(p.xyz.row(2).t(), p.phi * A, L);
}

void OMP_all_patches_ksvt(a_mat & alpha, a_mat & A, vector<patch_t> & patches, size_t M, size_t L)
{
	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
		if(patches[i].valid_xyz())
			OMP_patch(alpha, A, i, patches[i], L);
}

void KSVDT(a_mat & A, vector<patch_t> & patches, size_t M, size_t L)
{
	size_t K = A.n_rows;
	size_t m = A.n_cols;

	a_mat alpha(m, M);

	size_t iter = L;
	while(iter--)
	{
		OMP_all_patches_ksvt(alpha, A, patches, M, L);

		#pragma omp parallel for
		for(index_t j = 0; j < m; j++)
		{
			arma::uvec omega = find(abs(alpha.row(j)) > 0);

			a_mat sum(K, K, arma::fill::zeros);
			a_vec sum_error(K, arma::fill::zeros);

			for(arma::uword o: omega)
			{
				sum += alpha(j, o) * patches[o].phi.t() * patches[o].phi;

				a_mat D = patches[o].phi * A;
				D.col(j).zeros();
				a_mat e = patches[o].xyz.row(2).t() - D * alpha.col(o);

				sum_error += alpha(j, o) * patches[o].phi.t() * e;
			}
			if(omega.size())
			{
				a_vec X;
				solve(X, sum, sum_error);
				A.col(j) = X;
			}
		}
	}
}


} // namespace gproshan::mdict

