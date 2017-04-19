#include "fairing_spectral.h"
#include "laplacian.h"

fairing_spectral::fairing_spectral(size_t k_): fairing()
{
	k = k_;
}

fairing_spectral::~fairing_spectral()
{

}

void fairing_spectral::compute(che * shape)
{
	sp_mat L, A;
	double start_omp = omp_get_wtime();
	laplacian(shape, L, A);
	double time = omp_get_wtime() - start_omp;
	cout<<"time laplacian: "<<time<<endl;

	positions = new vertex[shape->n_vertices()];

	mat X(shape->n_vertices(), 3);
	for(index_t v = 0; v < shape->n_vertices(); v++)
	{
		X(v, 0) = shape->gt(v).x;
		X(v, 1) = shape->gt(v).y;
		X(v, 2) = shape->gt(v).z;
	}

	vec eigval;
	mat eigvec;
	start_omp = omp_get_wtime();
	eigs_sym(eigval, eigvec, L, k, "sm");
	time = omp_get_wtime() - start_omp;
	cout<<"time eigs: "<<time<<endl;

	X = eigvec * eigvec.t() * X;

	for(index_t v = 0; v < shape->n_vertices(); v++)
	{
		positions[v][0] = X(v, 0);
		positions[v][1] = X(v, 1);
		positions[v][2] = X(v, 2);
	}
}
