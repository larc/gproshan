#include "mdict.h"

#include "sampling.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {

const real_t sigma = 0.001;


// SPARSE

bool locval_t::operator < (const locval_t & lc)
{
	if(i == lc.i) return j < lc.j;
	return i < lc.i;
}

std::ostream & operator << (std::ostream & os, const locval_t & lc)
{
	return os << '(' << lc.i << ',' << lc.j << ") = " << lc.val; 
}

void OMP(vector<locval_t> & alpha, const a_vec & x, const index_t & i, const a_mat & D, const size_t & L)
{
	a_vec aa;
	arma::uvec selected_atoms;

	tie(aa, selected_atoms) = _OMP(x, D, L);

	for(index_t k = 0; k < selected_atoms.size(); k++)
	{
		#pragma omp critical
		alpha.push_back({selected_atoms(k), i, aa(k)});
	}
}

a_sp_mat OMP_all(vector<locval_t> & locval, const a_mat & X, const a_mat & D, const size_t & L)
{
	locval.clear();

	#pragma omp parallel for
	for(index_t i = 0; i < X.n_cols; i++)
		OMP(locval, X.col(i), i, D, L);

	arma::umat DI(2, locval.size());
	a_vec DV(locval.size());

	#pragma omp parallel for
	for(index_t k = 0; k < locval.size(); k++)
	{
		DI(0, k) = locval[k].i; // row
		DI(1, k) = locval[k].j; // column
		DV(k) = locval[k].val;
	}

	return a_sp_mat(DI, DV, D.n_cols, X.n_cols);
}

void sp_KSVD(a_mat & D, const a_mat & X, const size_t & L, size_t k)
{
	size_t m = D.n_cols;
	size_t M = X.n_cols;

	arma::uvec omega(M);
	a_mat R, E, U, V;
	a_vec s;

	vector<locval_t> locval;
	vector<size_t> rows(m + 1);
	index_t r;

	while(k--)
	{
		a_sp_mat alpha = OMP_all(locval, X, D, L);
		
		sort(locval.begin(), locval.end());

		rows[r = 0] = 0;
		for(index_t k = 1; k < locval.size(); k++)
			if(locval[k].i != locval[k - 1].i)
				rows[++r] = k;
		
		rows[++r] = locval.size();
		
		R = X - D * alpha;

		#pragma omp parallel for firstprivate(omega, E, U, V, s)
		for(index_t j = 0; j < m; j++)
		{
			for(index_t r = rows[j]; r < rows[j + 1]; r++)
				omega(r - rows[j]) = locval[r].j;

			if(rows[j + 1] - rows[j])
			{
				E = R + D.col(j) * alpha.row(j);
				E = E.cols(omega.head(rows[j + 1] - rows[j]));
				svd(U, s, V, E);
				D.col(j) = U.col(0);
			}
		}

		locval.clear();
	}
}


// DENSE

tuple<a_vec, arma::uvec> _OMP(const a_vec & x, const a_mat & D, const size_t & L)
{
	arma::uvec selected_atoms(L);
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

	return {aa, selected_atoms.head(l)};
}

a_vec OMP(const a_vec & x, const a_mat & D, const size_t & L)
{
	a_vec alpha(D.n_cols, arma::fill::zeros);
	a_vec aa;
	arma::uvec selected_atoms;

	tie(aa, selected_atoms) = _OMP(x, D, L);
	alpha.elem(selected_atoms) = aa;
	
	return alpha;
}

a_mat OMP_all(const a_mat & X, const a_mat & D, const size_t & L)
{
	a_mat alpha(D.n_cols, X.n_cols);

	#pragma omp parallel for
	for(index_t i = 0; i < X.n_cols; i++)
		alpha.col(i) = OMP(X.col(i), D, L);

	return alpha;
}

void KSVD(a_mat & D, const a_mat & X, const size_t & L, size_t k)
{
	arma::uvec omega;
	a_mat alpha, R, E, U, V;
	a_vec s;

	while(k--)
	{
		alpha = OMP_all(X, D, L);

		R = X - D * alpha;

		#pragma omp parallel for private(omega, E, U, V, s)
		for(index_t j = 0; j < D.n_cols; j++)
		{
			omega = find(abs(alpha.row(j)) > 0);
			if(omega.n_elem)
			{
				E = R + D.col(j) * alpha.row(j);
				E = E.cols(omega);
				svd(U, s, V, E);
				D.col(j) = U.col(0);
			}
		}
	}
}


// MESH DENSE 

a_vec OMP(const a_mat & A, const patch & p, const size_t & L)
{
	return OMP(p.xyz.row(2).t(), p.phi * A, L);
}

a_mat OMP_all(const vector<patch> & patches, const a_mat & A, const size_t & L)
{
	a_mat alpha(A.n_cols, patches.size());

	#pragma omp parallel for
	for(index_t i = 0; i < patches.size(); i++)
		alpha.col(i) = OMP(A, patches[i], L);
	
	return alpha;
}

void KSVD(a_mat & A, const vector<patch> & patches, const size_t & L, size_t k)
{
	size_t K = A.n_rows;

	a_mat alpha;

	while(k--)
	{
		alpha = OMP_all(patches, A, L);

		#pragma omp parallel for
		for(index_t j = 0; j < A.n_cols; j++)
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

