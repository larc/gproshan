#include "mdict.h"

#include "sampling.h"

#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {

const real_t sigma = 0.01;


// SPARSE

bool operator < (const locval_t & a, const locval_t & b)
{
	return (a.i == b.i) ? a.j < b.j : a.i < b.i;
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
	vector<size_t> rows;
	rows.reserve(m + 1);

	while(k--)
	{
		a_sp_mat alpha = OMP_all(locval, X, D, L);
		
		sort(locval.begin(), locval.end());

		rows.push_back(0);
		for(index_t k = 1; k < locval.size(); k++)
			if(locval[k].i != locval[k - 1].i)
				rows.push_back(k);
		
		rows.push_back(locval.size());
		
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

arma::uword max_index(const a_vec & V,const arma::uchar_vec & mask)
{

	arma::uvec indices = arma::sort_index( V , "desscend");

	for(size_t i=0; i< V.size(); i++)
		if(mask[indices [i]]) return indices[i];
	
}

tuple<a_vec, arma::uvec> _OMP(const a_vec & x, const a_mat & D, const size_t & L, const arma::uchar_vec & mask)
{

	arma::uvec selected_atoms(L);
	real_t threshold = norm(x) * sigma;	

	a_mat DD;
	a_vec aa, r = x;	
	
	index_t l = 0;
	while(norm(r) > threshold && l < L)
	{
	//	gproshan_debug_var(D.t() * r);
		selected_atoms(l) = max_index(abs(D.t() * r), mask);
		
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

a_vec OMP(const a_vec & x, const a_mat & D, const size_t & L, const arma::uchar_vec & mask)
{
	a_vec alpha(D.n_cols, arma::fill::zeros);
	a_vec aa;
	arma::uvec selected_atoms;

	tie(aa, selected_atoms) = _OMP(x, D, L, mask);
	alpha.elem(selected_atoms) = aa;

	return alpha;
}

a_mat OMP_all(const a_mat & X, const a_mat & D, const size_t & L)
{
	a_mat alpha(D.n_cols, X.n_cols);
	#pragma omp parallel for
	for(index_t i = 0; i < X.n_cols; i++)
	{
		alpha.col(i) = OMP(X.col(i), D, L);
	}
		
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

a_vec OMP(const patch & p, const a_mat & A, const size_t & L)
{
	
	return OMP(p.xyz.row(2).t(), p.phi * A, L);
}
a_vec OMP(const patch & p, basis * phi_basis, const a_mat & A, const size_t & L)
{
	arma::uchar_vec mask;
	mask.zeros(A.n_cols);
	for(int i = 0; i < A.n_cols; i++)
		if(  phi_basis->get_frequency(i) >= 0*p.avg_dist) //2.5* if it ismin
		{
		//	gproshan_debug_var(phi_basis->get_frequency(i));
		//	gproshan_debug_var(p.avg_dist);
			 mask(i) = 1;
		}
	return OMP(p.xyz.row(2).t(), p.phi * A, L, mask);
}

a_mat OMP_all(const vector<patch> & patches, basis * phi_basis, const a_mat & A, const size_t & L)
{
	a_mat alpha(A.n_cols, patches.size());

	#pragma omp parallel for
	for(index_t i = 0; i < patches.size(); i++)
		alpha.col(i) = OMP(patches[i],phi_basis, A, L);
	
	return alpha;
}

a_mat OMP_all(const vector<patch> & patches, const a_mat & A, const size_t & L)
{
	a_mat alpha(A.n_cols, patches.size());

	#pragma omp parallel for
	for(index_t i = 0; i < patches.size(); i++)
		alpha.col(i) = OMP(patches[i], A, L);
	
	return alpha;
}

void KSVD(a_mat & A, const vector<patch> & patches, const size_t & L, size_t k)
{
	size_t K = A.n_rows;

	arma::uvec omega;

	a_mat new_A = A;
	a_mat alpha, D, sum, sum_error;
	a_vec a, e;
	
	real_t aj;

	while(k--)
	{
		alpha = OMP_all(patches, A, L);

		#pragma omp parallel for private(omega, a, aj, D, e, sum, sum_error)
		for(index_t j = 0; j < A.n_cols; j++)
		{
			//Taking all alphas that uses atom j
			arma::uvec omega = find(abs(alpha.row(j)) > 0);

			sum.zeros(K, K);
			sum_error.zeros(K);

			for(arma::uword & i: omega)
			{
				a = alpha.col(i);
				a(j) = 0;

				D = patches[i].phi * A; // fetch the discrete dictionary for the patch i
				e = patches[i].xyz.row(2).t() - D * a; // getting the rec error for the patch i 
				aj = as_scalar(e.t() * D.col(j) / (D.col(j).t() * D.col(j)));

				sum += aj * aj * patches[i].phi.t() * patches[i].phi;
				sum_error += aj * patches[i].phi.t() * e;
				//concat e patches[i].phi.t() * e;
				//apply svd to update the atom
				
			}

			if(omega.size())
				new_A.col(j) = solve(sum, sum_error);
		}

		A = new_A;
	}
}


// MESH SPARSE

void OMP(vector<locval_t> & alpha, const patch & p, const index_t & i, const a_mat & A, const size_t & L)
{
	OMP(alpha, p.xyz.row(2).t(), i, p.phi * A, L);
}

a_sp_mat OMP_all(vector<locval_t> & locval, const vector<patch> & patches, const a_mat & A, const size_t & L)
{
	locval.clear();

	#pragma omp parallel for
	for(index_t i = 0; i < patches.size(); i++)
		OMP(locval, patches[i], i, A, L);

	arma::umat DI(2, locval.size());
	a_vec DV(locval.size());

	#pragma omp parallel for
	for(index_t k = 0; k < locval.size(); k++)
	{
		DI(0, k) = locval[k].i; // row
		DI(1, k) = locval[k].j; // column
		DV(k) = locval[k].val;
	}

	return a_sp_mat(DI, DV, A.n_cols, patches.size());
}

void sp_KSVD(a_mat & A, const vector<patch> & patches, const size_t & L, size_t k)
{
	size_t K = A.n_rows;

	a_mat new_A = A;
	a_mat D, sum, sum_error;
	a_vec a, e;
	
	real_t aj;
	
	vector<locval_t> locval;
	vector<size_t> rows;
	rows.reserve(A.n_cols + 1);

	while(k--)
	{
		a_sp_mat alpha = OMP_all(locval, patches, A, L);
		
		sort(locval.begin(), locval.end());

		rows.push_back(0);
		for(index_t k = 1; k < locval.size(); k++)
			if(locval[k].i != locval[k - 1].i)
				rows.push_back(k);
		
		rows.push_back(locval.size());

		#pragma omp parallel for private(a, aj, D, e, sum, sum_error)
		for(index_t j = 0; j < A.n_cols; j++)
		{
			sum.zeros(K, K);
			sum_error.zeros(K);

			for(index_t r = rows[j]; r < rows[j + 1]; r++)
			{
				const index_t & i = locval[r].j;

				a = alpha.col(i);
				a(j) = 0;

				D = patches[i].phi * A;
				e = patches[i].xyz.row(2).t() - D * a;
				aj = as_scalar(e.t() * D.col(j) / (D.col(j).t() * D.col(j)));

				sum += aj * aj * patches[i].phi.t() * patches[i].phi;
				sum_error += aj * patches[i].phi.t() * e;
			}

			if(rows[j + 1] - rows[j])
				new_A.col(j) = solve(sum, sum_error);
		}

		A = new_A;
	}
}


} // namespace gproshan::mdict

