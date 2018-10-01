#include "include_arma.h"
#include <cusolverSp.h>

int solve_positive_definite_gpu(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx)
{
	cudaDeviceReset();

	// device sparse matrix A to device (CSC format)
	int * dA_col_ptrs, * dA_row_indices;
	real_t * dA_values;
	
	cudaMalloc(&dA_col_ptrs, (m + 1) * sizeof(int));
	cudaMemcpy(dA_col_ptrs, hA_col_ptrs, (m + 1) * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&dA_row_indices, nnz * sizeof(int));
	cudaMemcpy(dA_row_indices, hA_row_indices, nnz * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&dA_values, nnz * sizeof(real_t));
	cudaMemcpy(dA_values, hA_values, nnz * sizeof(real_t), cudaMemcpyHostToDevice); 
	
	// vector b to device
	real_t * db;
	cudaMalloc(&db, nnz * sizeof(real_t));
	cudaMemcpy(db, hb, nnz * sizeof(real_t), cudaMemcpyHostToDevice);

	// vector x to device
	real_t * dx;
	cudaMalloc(&dx, m * sizeof(real_t));

	// solve Ax = b

	int singularity;
	
	cusolverSpHandle_t handle_cusolver;
	cusolverSpCreate(&handle_cusolver);

	cusparseMatDescr_t descr = 0;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
	
#ifdef SINGLE_P
	cusolverSpScsrlsvchol(handle_cusolver, m, nnz, descr, dA_values, dA_col_ptrs, dA_row_indices, db, 0, 0, dx, &singularity);
#else
	cusolverSpDcsrlsvchol(handle_cusolver, m, nnz, descr, dA_values, dA_col_ptrs, dA_row_indices, db, 0, 0, dx, &singularity);
#endif
	
	cusparseDestroyMatDescr(descr);
	cusolverSpDestroy(handle_cusolver);
	
	// copy dx to host x
	cudaMemcpy(hx, dx, m * sizeof(real_t), cudaMemcpyDeviceToHost);

	// free device memory
	cudaFree(dA_col_ptrs);
	cudaFree(dA_row_indices);
	cudaFree(dA_values);
	cudaFree(db);
	cudaFree(dx);
	
	return singularity;
}

