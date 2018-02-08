#include "d_dict_learning.h"
#include "sampling.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>

// mesh dictionary learning and sparse coding namespace
namespace mdict {

void OMP(vec & alpha, vec & x, mat & D, size_t L)
{
	size_t n = x.n_elem;
	size_t m = D.n_cols;

	uword max_i;

	alpha = zeros(m);
	uvec selected_atoms(L, fill::zeros);
	vec aa;

	vec r = x;
	for(index_t l = 0; l < L; l++)
	{
		vec Dtr = abs(D.t() * r);

		Dtr.max(max_i);
		selected_atoms(l) = max_i;

		mat DD = D.cols(selected_atoms.head(l + 1));
		aa = pinv(DD) * x;
		r = x - DD * aa;
	}

	alpha.elem(selected_atoms) = aa;
}

void KSVD(mat & D, mat & X, size_t L)
{
	size_t n = X.n_rows;
	size_t m = D.n_cols;
	size_t M = X.n_cols;

	mat alpha(m, M);

	size_t iter = L;
	while(iter--)
	{
		#pragma omp parallel for
		for(index_t i = 0; i < M; i++)
		{
			vec a;
			vec x = X.col(i);
			OMP(a, x, D, L);
			alpha.col(i) = a;
		}

		#pragma omp parallel for
		for(index_t j = 0; j < m; j++)
		{
			uvec omega = find(abs(alpha.row(j)) > 0);
			if(omega.n_elem)
			{
				mat E = X - D * alpha - D.col(j) * alpha.row(j);
				E = E.cols(omega);
				mat U;
				vec s;
				mat V;
				svd(U, s, V, E);
				D.col(j) = U.col(0);
				rowvec a = alpha.row(j);
				a.elem(omega) = s(0) * V.col(0);
				alpha.row(j) = a;
			}
			//else
			//	D.col(j).randu();
		}
	}
}

void OMP_patch(mat & alpha, const mat & A, const index_t & i, patch_t & p, const size_t & L)
{
	vec a;
	vec x = p.xyz.row(2).t();
	mat D = p.phi * A;

	OMP(a, x, D, L);
	alpha.col(i) = a;
}

void OMP_all_patches_ksvt(mat & alpha, mat & A, vector<patch_t> & patches, size_t M, size_t L)
{
	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
		if(patches[i].valid_xyz())
			OMP_patch(alpha, A, i, patches[i], L);
}

void KSVDT(mat & A, vector<patch_t> & patches, size_t M, size_t L)
{
	size_t K = A.n_rows;
	size_t m = A.n_cols;

	mat alpha(m, M);

	size_t iter = L;
	while(iter--)
	{
		OMP_all_patches_ksvt(alpha, A, patches, M, L);

		#pragma omp parallel for
		for(index_t j = 0; j < m; j++)
		{
			uvec omega = find(abs(alpha.row(j)) > 0);

			mat sum(K, K, fill::zeros);
			vec sum_error(K, fill::zeros);

			for(uword o: omega)
			{
				sum += alpha(j, o) * patches[o].phi.t() * patches[o].phi;

				mat D = patches[o].phi * A;
				D.col(j).zeros();
				mat e = patches[o].xyz.row(2).t() - D * alpha.col(o);

				sum_error += alpha(j, o) * patches[o].phi.t() * e;
			}
			if(omega.size())
			{
				vec X;
				solve(X, sum, sum_error);
				A.col(j) = X;
			}
		}
	}
}

} // mdict

