#include "heat_flow.h"

#include "laplacian.h"

#include <cassert>

distance_t * heat_flow(che * mesh, const vector<index_t> & sources)
{
	if(!sources.size()) return 0;
	
	// build impulse signal
	a_mat u0(mesh->n_vertices(), 1, arma::fill::zeros);
	for(auto & v: sources) u0(v) = 1;
	
	// step
	real_t dt = mesh->mean_edge();
	dt *= dt;

	debug(dt)

	a_sp_mat L, A;
	laplacian(mesh, L, A);
	
	// make L positive-definite
	L += 1.0e-8 * A;

	// heat flow for short interval
	A += dt * L;
	a_mat u;
	
	debug(L.col_ptrs)
	assert(spsolve(u, A, u0));

	// extract geodesics
	distance_t * distances = new distance_t[mesh->n_vertices()];
	
	a_mat div;
	compute_divergence(mesh, u, div);

	a_mat phi(distances, mesh->n_vertices(), 1, false);
	assert(spsolve(phi, L, div));
	
	real_t mp = phi.min();
	phi.for_each([&mp](a_mat::elem_type & val) { val -= mp; });

	return distances;
}

cholmod_dense * arma_2_cholmod(const a_mat & m)
{
	cholmod_dense * cm;

	#ifdef SINGLE_P
	#else
//		cm = cholmod_l_allocate_dense( m, n, d, CHOLMOD_REAL, context );
	#endif

	memcpy(cm, m.memptr(), m.size());

	return cm;
}

void compute_divergence(che * mesh, const a_mat & u, a_mat & div)
{
	div.resize(mesh->n_vertices(), 1);
	
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		real_t & sum = div(v);

		sum = 0;
		for_star(he, mesh, v)
			sum += ( - mesh->gradient_he(he, u.memptr()) , mesh->gt_vt(he) - mesh->gt_vt(next(he)) );
	}
}

