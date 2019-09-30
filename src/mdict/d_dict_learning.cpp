#include "d_dict_learning.h"
#include "sampling.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>



// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void OMP(a_vec & alpha, a_vec & x, a_mat & D, size_t L)
{
//	size_t n = x.n_elem;
	size_t m = D.n_cols;

	arma::uword max_i;

	alpha.zeros(m);
	arma::uvec selected_atoms(L, arma::fill::zeros);
	a_vec aa;

	double sigma = 0.001;
	double threshold = norm(x) * sigma;

	a_vec r = x;
	for(index_t l = 0; norm(r) > threshold && l < L; l++)
	{
		a_vec Dtr = abs(D.t() * r);

		Dtr.max(max_i);
		selected_atoms(l) = max_i;

		a_mat DD = D.cols(selected_atoms.head(l + 1));
		aa = pinv(DD) * x;
		r = x - DD * aa;
	}

	alpha.elem(selected_atoms) = aa;
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
		{
			a_vec a;
			a_vec x = X.col(i);
			OMP(a, x, D, L);
			alpha.col(i) = a;
		}

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
	a_vec a;
	a_vec x = p.xyz.row(2).t();
	a_mat D = p.phi * A;

	OMP(a, x, D, L);
	alpha.col(i) = a;
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

/// DEPRECATED
void OMP_patch(a_mat & alpha, const a_mat & A, const index_t & i, patch_t & p, const size_t & L)
{
	a_vec a;
	a_vec x = p.xyz.row(2).t();
	a_mat D = p.phi * A;

	OMP(a, x, D, L);
	alpha.col(i) = a;
}

/// DEPRECATED
void OMP_all_patches_ksvt(a_mat & alpha, a_mat & A, vector<patch_t> & patches, size_t M, size_t L)
{
	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
		if(patches[i].valid_xyz())
			OMP_patch(alpha, A, i, patches[i], L);
}

/// DEPRECATED
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

