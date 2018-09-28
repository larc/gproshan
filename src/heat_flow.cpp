#include "heat_flow.h"

#include "laplacian.h"

#include <cassert>

cholmod_common context;

distance_t * heat_flow(che * mesh, const vector<index_t> & sources)
{
	if(!sources.size()) return 0;
	
	// build impulse signal
	a_mat u0(mesh->n_vertices(), 1, arma::fill::zeros);
	for(auto & v: sources) u0(v) = 1;
	
	// step
	real_t dt = mesh->mean_edge();
	dt *= dt;

	a_sp_mat L, A;
	laplacian(mesh, L, A);
	
	// make L positive-definite
	L += 1.0e-8 * A;

	// heat flow for short interval
	A += dt * L;
	a_mat u(mesh->n_vertices(), 1);
	
	float time;
	TIC(time)
	cholmod_l_start(&context);
	#ifdef SINGLE_P
		context.dtype = CHOLMOD_SINGLE;
	#endif
	
	solve_positive_definite(u, A, u0);
	//assert(spsolve(u, A, u0));

	// extract geodesics
	distance_t * distances = new distance_t[mesh->n_vertices()];
	
	a_mat div(mesh->n_vertices(), 1);
	compute_divergence(mesh, u, div);

	a_mat phi(distances, mesh->n_vertices(), 1, false);

	solve_positive_definite(phi, L, div);
	//assert(spsolve(phi, L, div));
	
	real_t min_val = phi.min();
	phi.for_each([&min_val](a_mat::elem_type & val) { val -= min_val; });

	debug(phi.max())

	cholmod_l_gpu_stats(&context);
	
	cholmod_l_finish(&context);
	TOC(time)
	debug(time)

	return distances;
}

void compute_divergence(che * mesh, const a_mat & u, a_mat & div)
{
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		real_t & sum = div(v);

		sum = 0;
		for_star(he, mesh, v)
			sum += (
					mesh->normal_he(he) * ( mesh->gt_vt(prev(he)) - mesh->gt_vt(next(he)) ) ,
					- mesh->gradient_he(he, u.memptr()) 
					);
	}
}

void solve_positive_definite(a_mat & x, const a_sp_mat & A, const a_mat & b)
{
	cholmod_sparse * cA = arma_2_cholmod(A); cA->stype = 1;
	cholmod_dense * cb = arma_2_cholmod(b);

	cholmod_factor * L = cholmod_l_analyze(cA, &context);
	cholmod_l_factorize(cA, L, &context);
	cholmod_dense * cx = cholmod_l_solve(CHOLMOD_A, L, cb, &context);

	assert(x.n_rows == b.n_rows);
	memcpy(x.memptr(), cx->x, x.n_rows * sizeof(real_t));

	cholmod_l_free_factor(&L, &context);
	cholmod_l_free_sparse(&cA, &context);
	cholmod_l_free_dense(&cb, &context);
}

cholmod_dense * arma_2_cholmod(const a_mat & D)
{
	cholmod_dense * cD = cholmod_l_allocate_dense(D.n_rows, D.n_cols, D.n_rows, CHOLMOD_REAL, &context);
	memcpy(cD->x, D.memptr(), D.n_elem * sizeof(real_t));

	return cD;
}

cholmod_sparse * arma_2_cholmod(const a_sp_mat & S)
{
	assert(sizeof(arma::uword) == sizeof(SuiteSparse_long));
	
	cholmod_sparse * cS = cholmod_l_allocate_sparse(S.n_rows, S.n_cols, S.n_nonzero, 1, 1, 0, CHOLMOD_REAL, &context);
	
	memcpy(cS->p, S.col_ptrs, (S.n_cols + 1) * sizeof(arma::uword));
	memcpy(cS->i, S.row_indices, S.n_nonzero * sizeof(arma::uword));
	memcpy(cS->x, S.values, S.n_nonzero * sizeof(real_t));
	
	return cS;
}

