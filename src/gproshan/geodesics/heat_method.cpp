#include <gproshan/geodesics/heat_method.h>

#include <gproshan/util.h>
#include <gproshan/laplacian/laplacian.h>

#include <cassert>
#include <cholmod.h>


// geometry processing and shape analysis framework
namespace gproshan {


/// cholmod Keenan implementation
/// base on the code https://github.com/larc/dgpdec-course/tree/master/Geodesics
double solve_positive_definite(arma::mat & x, const arma::sp_mat & A, const arma::mat & b, cholmod_common * context);

double heat_method(float * dist, const che * mesh, const std::vector<index_t> & sources, const heat_method_opt & opt)
{
	if(!size(sources)) return 0;

	// build impulse signal
	arma::vec u0(mesh->n_vertices, arma::fill::zeros);
	for(auto & v: sources) u0(v) = 1;

	// step
	float dt = mesh->mean_edge();
	dt *= dt;

	arma::sp_mat L, A;
	laplacian(mesh, L, A);

	// make L positive-definite
	L += 1e-8 * A;

	// heat flow for short interval
	A += dt * L;


	double solve_time = 0;
	arma::vec u(mesh->n_vertices);

	cholmod_common context;
	switch(opt)
	{
		case HEAT_ARMA:
			if(!spsolve(u, A, u0)) gproshan_error(arma no solution);
			break;
		case HEAT_CHOLMOD:
			cholmod_l_start(&context);
			solve_time += solve_positive_definite(u, A, u0, &context);
			break;
	#ifdef GPROSHAN_CUDA
		case HEAT_CUDA:
			solve_time += solve_positive_definite_gpu(u, A, u0);
			break;
	#endif // GPROSHAN_CUDA
	}


	// extract geodesics

	arma::vec div = compute_divergence(mesh, u);
	arma::vec phi(mesh->n_vertices);

	switch(opt)
	{
		case HEAT_ARMA:
			if(!spsolve(phi, L, div)) gproshan_error(arma no solution);
			break;
		case HEAT_CHOLMOD:
			solve_time += solve_positive_definite(phi, L, div, &context);
			cholmod_l_finish(&context);
			break;
	#ifdef GPROSHAN_CUDA
		case HEAT_CUDA:
			solve_time += solve_positive_definite_gpu(phi, L, div);
			break;
	#endif // GPROSHAN_CUDA
	}

	if(phi.size() == mesh->n_vertices)
	{
		phi -= phi.min();
		phi *= 0.5;

		#pragma omp parallel for
		for(index_t v = 0; v < mesh->n_vertices; ++v)
			dist[v] = phi(v);
	}

	return solve_time;
}

arma::vec compute_divergence(const che * mesh, const arma::vec & u)
{
	arma::vec div(mesh->n_vertices);

	std::vector<float> f(mesh->n_vertices);

	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		float sum = 0;
		for(const index_t he: mesh->star(v))
		{
			const vertex & n = mesh->normal_he(he);
			const vertex & v = mesh->vertex_he(he_prev(he)) - mesh->vertex_he(he_next(he));

			const dvec3 nhe = {n.x(), n.y(), n.z()};
			const dvec3 vhe = {v.x(), v.y(), v.z()};
			const dvec3 ghe = mesh->gradient_he(he, u.memptr());
			sum += dot(cross(nhe, vhe), -ghe);
		}

		div(v) = sum;
	}

	return div;
}


cholmod_dense * arma_2_cholmod(const arma::mat & m, cholmod_common * context);
cholmod_sparse * arma_2_cholmod(const arma::sp_mat & m, cholmod_common * context);

double solve_positive_definite(arma::mat & x, const arma::sp_mat & A, const arma::mat & b, cholmod_common * context)
{
	cholmod_sparse * cA = arma_2_cholmod(A, context);
	cA->stype = 1;

	cholmod_dense * cb = arma_2_cholmod(b, context);

	cholmod_factor * L = cholmod_l_analyze(cA, context);
	cholmod_l_factorize(cA, L, context);

	/* fill ratio
	gproshan_debug_var(L->xsize);
	gproshan_debug_var(cA->nzmax);
	gproshan_debug_var(L->xsize / cA->nzmax);
	*/

	double solve_time;
	TIC(solve_time)
	cholmod_dense * cx = cholmod_l_solve(CHOLMOD_A, L, cb, context);
	TOC(solve_time)

	assert(x.n_rows == b.n_rows);
	memcpy(x.memptr(), cx->x, x.n_rows * sizeof(double));

	cholmod_l_free_factor(&L, context);
	cholmod_l_free_sparse(&cA, context);
	cholmod_l_free_dense(&cb, context);

	return solve_time;
}

cholmod_dense * arma_2_cholmod(const arma::mat & D, cholmod_common * context)
{
	cholmod_dense * cD = cholmod_l_allocate_dense(D.n_rows, D.n_cols, D.n_rows, CHOLMOD_REAL, context);
	memcpy(cD->x, D.memptr(), D.n_elem * sizeof(double));

	return cD;
}

cholmod_sparse * arma_2_cholmod(const arma::sp_mat & S, cholmod_common * context)
{
	assert(sizeof(arma::uword) == sizeof(SuiteSparse_long));

	cholmod_sparse * cS = cholmod_l_allocate_sparse(S.n_rows, S.n_cols, S.n_nonzero, 1, 1, 0, CHOLMOD_REAL, context);

	memcpy(cS->p, S.col_ptrs, (S.n_cols + 1) * sizeof(arma::uword));
	memcpy(cS->i, S.row_indices, S.n_nonzero * sizeof(arma::uword));
	memcpy(cS->x, S.values, S.n_nonzero * sizeof(double));

	return cS;
}


#ifdef GPROSHAN_CUDA

double solve_positive_definite_gpu(arma::mat & x, const arma::sp_mat & A, const arma::mat & b)
{
	int * hA_col_ptrs = new int[A.n_cols + 1];
	int * hA_row_indices = new int[A.n_nonzero];

	#pragma omp parallel for
	for(index_t i = 0; i <= A.n_cols; ++i)
		hA_col_ptrs[i] = A.col_ptrs[i];

	#pragma omp parallel for
	for(index_t i = 0; i < A.n_nonzero; ++i)
		hA_row_indices[i] = A.row_indices[i];

	double solve_time = solve_positive_definite_cusolver(A.n_rows, A.n_nonzero, A.values, hA_col_ptrs, hA_row_indices, b.memptr(), x.memptr());

	delete [] hA_col_ptrs;
	delete [] hA_row_indices;

	return solve_time;
}

#endif // GPROSHAN_CUDA


} // namespace gproshan

