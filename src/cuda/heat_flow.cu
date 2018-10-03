#include "include_arma.h"

#include <cassert>

#include <omp.h>
#include <cusolverSp.h>

double solve_positive_definite_gpu(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx)
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
	
	// aux vector y to device
	real_t * dy;
	cudaMalloc(&dy, m * sizeof(real_t));
	
	cusparseHandle_t handle;
	cusparseCreate(&handle);

	// SOLVE Ax = b
	double solve_time;

	cusparseMatDescr_t descr_M = 0;
	cusparseMatDescr_t descr_L = 0;
	
	csric02Info_t info_M = 0;
	csrsv2Info_t info_L = 0;
	csrsv2Info_t info_Lt = 0;
	
	int buffer_size_M;
	int buffer_size_L;
	int buffer_size_Lt;
	int buffer_size;
	
	void * buffer = 0;

	int structural_zero;
	int numerical_zero;

	const real_t alpha = 1.;
	const cusparseSolvePolicy_t policy_M  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
	const cusparseSolvePolicy_t policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
	const cusparseSolvePolicy_t policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
	const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
	const cusparseOperation_t trans_Lt = CUSPARSE_OPERATION_TRANSPOSE;

	cusparseCreateMatDescr(&descr_M);
	cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);

	cusparseCreateMatDescr(&descr_L);
	cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
	cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);

	cusparseCreateCsric02Info(&info_M);
	cusparseCreateCsrsv2Info(&info_L);
	cusparseCreateCsrsv2Info(&info_Lt);

#ifdef SINGLE_P
#else
	cusparseDcsric02_bufferSize(handle, m, nnz, descr_M, dA_values, dA_col_ptrs, dA_row_indices, info_M, &buffer_size_M);
	cusparseDcsrsv2_bufferSize(handle, trans_L, m, nnz, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_L, &buffer_size_L);
	cusparseDcsrsv2_bufferSize(handle, trans_Lt, m, nnz, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_Lt, &buffer_size_Lt);
#endif

	buffer_size = max(buffer_size_M, max(buffer_size_L, buffer_size_Lt));
	printf("%d\n", buffer_size);
	cudaMalloc(&buffer, buffer_size);

#ifdef SINGLE_P
#else
	cusparseDcsric02_analysis(handle, m, nnz, descr_M, dA_values, dA_col_ptrs, dA_row_indices, info_M, policy_M, buffer);
#endif
	if(CUSPARSE_STATUS_ZERO_PIVOT == cusparseXcsric02_zeroPivot(handle, info_M, &structural_zero))
		printf("A(%d,%d) is missing\n", structural_zero, structural_zero);

#ifdef SINGLE_P
#else
	cusparseDcsrsv2_analysis(handle, trans_L, m, nnz, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_L, policy_L, buffer);
	cusparseDcsrsv2_analysis(handle, trans_Lt, m, nnz, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_Lt, policy_Lt, buffer);

	cusparseDcsric02(handle, m, nnz, descr_M, dA_values, dA_col_ptrs, dA_row_indices, info_M, policy_M, buffer);
#endif
	if(CUSPARSE_STATUS_ZERO_PIVOT == cusparseXcsric02_zeroPivot(handle, info_M, &numerical_zero))
		printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);

	TIC(solve_time)

#ifdef SINGLE_P
#else
	assert(CUSPARSE_STATUS_SUCCESS == cusparseDcsrsv2_solve(handle, trans_L, m, nnz, &alpha, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_L, db, dy, policy_L, buffer));
	assert(CUSPARSE_STATUS_SUCCESS == cusparseDcsrsv2_solve(handle, trans_Lt, m, nnz, &alpha, descr_L, dA_values, dA_col_ptrs, dA_row_indices, info_Lt, dy, dx, policy_Lt, buffer));
#endif
	
	// copy sol x to host
	cudaMemcpy(hx, dx, m * sizeof(real_t), cudaMemcpyDeviceToHost);

	TOC(solve_time)

	// FREE
	cudaFree(buffer);
	cusparseDestroyMatDescr(descr_M);
	cusparseDestroyMatDescr(descr_L);
	cusparseDestroyCsric02Info(info_M);
	cusparseDestroyCsrsv2Info(info_L);
	cusparseDestroyCsrsv2Info(info_Lt);
	cusparseDestroy(handle);

	cudaFree(dA_col_ptrs);
	cudaFree(dA_row_indices);
	cudaFree(dA_values);
	cudaFree(db);
	cudaFree(dx);
	cudaFree(dx);
	
	return solve_time;
}

