#include "include_arma.h"

#include <cassert>

#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>


// geometry processing and shape analysis framework
namespace gproshan {


struct cu_spAxb
{
	int * A_col_ptrs, * A_row_indices;
	real_t * A_values, * x, * b;
	
	cu_spAxb(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx)
	{
		cudaMalloc(&A_col_ptrs, (m + 1) * sizeof(int));
		cudaMemcpy(A_col_ptrs, hA_col_ptrs, (m + 1) * sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc(&A_row_indices, nnz * sizeof(int));
		cudaMemcpy(A_row_indices, hA_row_indices, nnz * sizeof(int), cudaMemcpyHostToDevice);

		cudaMalloc(&A_values, nnz * sizeof(real_t));
		cudaMemcpy(A_values, hA_values, nnz * sizeof(real_t), cudaMemcpyHostToDevice); 
		
		cudaMalloc(&b, nnz * sizeof(real_t));
		cudaMemcpy(b, hb, nnz * sizeof(real_t), cudaMemcpyHostToDevice);

		cudaMalloc(&x, m * sizeof(real_t));
	}

	~cu_spAxb()
	{
		cudaFree(A_col_ptrs);
		cudaFree(A_row_indices);
		cudaFree(A_values);
		cudaFree(b);
		cudaFree(x);
	}
};

double solve_positive_definite_cusolver(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx, const bool host)
{
	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	
	// solve Ax = b

	int singularity;
	
	cusolverSpHandle_t handle_cusolver;
	cusolverSpCreate(&handle_cusolver);

	cusparseMatDescr_t descr = 0;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
	
	if(host)
	{
		#ifdef SINGLE_P
			cusolverSpScsrlsvcholHost(handle_cusolver, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, hb, 0, 0, hx, &singularity);
		#else
			cusolverSpDcsrlsvcholHost(handle_cusolver, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, hb, 0, 0, hx, &singularity);
		#endif
	}
	else
	{
		// allocate A, x, b into device
		cu_spAxb data(m, nnz, hA_values, hA_col_ptrs, hA_row_indices, hb, hx);

		cusolverStatus_t status;
		#ifdef SINGLE_P
			status = cusolverSpScsrlsvchol(handle_cusolver, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, data.b, 0, 0, data.x, &singularity);
		#else
			status = cusolverSpDcsrlsvchol(handle_cusolver, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, data.b, 0, 0, data.x, &singularity);
		#endif

		if(status == CUSOLVER_STATUS_SUCCESS)
			cudaMemcpy(hx, data.x, m * sizeof(real_t), cudaMemcpyDeviceToHost);
		else
			memset(hx, 0, m * sizeof(real_t));
	}

//	printf("%d\n", singularity != -1);

	cusparseDestroyMatDescr(descr);
	cusolverSpDestroy(handle_cusolver);

	// end Ax = b

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return (double) time / 1000;
}

double solve_positive_definite_cusparse(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx)
{
	cudaDeviceReset();
	
	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// allocate A, x, b into device
	cu_spAxb data(m, nnz, hA_values, hA_col_ptrs, hA_row_indices, hb, hx);
	
	// aux vector y to device
	real_t * dy;
	cudaMalloc(&dy, m * sizeof(real_t));
	
	cusparseHandle_t handle;
	cusparseCreate(&handle);

	// SOLVE Ax = b
	
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
		cusparseScsric02_bufferSize(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, &buffer_size_M);
		cusparseScsrsv2_bufferSize(handle, trans_L, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, &buffer_size_L);
		cusparseScsrsv2_bufferSize(handle, trans_Lt, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, &buffer_size_Lt);
	#else
		cusparseDcsric02_bufferSize(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, &buffer_size_M);
		cusparseDcsrsv2_bufferSize(handle, trans_L, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, &buffer_size_L);
		cusparseDcsrsv2_bufferSize(handle, trans_Lt, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, &buffer_size_Lt);
	#endif

	buffer_size = max(buffer_size_M, max(buffer_size_L, buffer_size_Lt));
	cudaMalloc(&buffer, buffer_size);

	#ifdef SINGLE_P
		cusparseScsric02_analysis(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, policy_M, buffer);
	#else
		cusparseDcsric02_analysis(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, policy_M, buffer);
	#endif
	if(CUSPARSE_STATUS_ZERO_PIVOT == cusparseXcsric02_zeroPivot(handle, info_M, &structural_zero))
		printf("A(%d,%d) is missing\n", structural_zero, structural_zero);

	#ifdef SINGLE_P
		cusparseScsrsv2_analysis(handle, trans_L, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, policy_L, buffer);
		cusparseScsrsv2_analysis(handle, trans_Lt, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, policy_Lt, buffer);

		cusparseScsric02(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, policy_M, buffer);
	#else
		cusparseDcsrsv2_analysis(handle, trans_L, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, policy_L, buffer);
		cusparseDcsrsv2_analysis(handle, trans_Lt, m, nnz, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, policy_Lt, buffer);

		cusparseDcsric02(handle, m, nnz, descr_M, data.A_values, data.A_col_ptrs, data.A_row_indices, info_M, policy_M, buffer);
	#endif
	if(CUSPARSE_STATUS_ZERO_PIVOT == cusparseXcsric02_zeroPivot(handle, info_M, &numerical_zero))
		printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);


	// SOLVE
	cudaEventRecord(start, 0);
	
	#ifdef SINGLE_P
		cusparseScsrsv2_solve(handle, trans_L, m, nnz, &alpha, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, data.b, dy, policy_L, buffer);
		cusparseScsrsv2_solve(handle, trans_Lt, m, nnz, &alpha, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, dy, data.x, policy_Lt, buffer);
	#else
		cusparseDcsrsv2_solve(handle, trans_L, m, nnz, &alpha, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_L, data.b, dy, policy_L, buffer);
		cusparseDcsrsv2_solve(handle, trans_Lt, m, nnz, &alpha, descr_L, data.A_values, data.A_col_ptrs, data.A_row_indices, info_Lt, dy, data.x, policy_Lt, buffer);
	#endif
	
	// copy sol x to host
	cudaMemcpy(hx, data.x, m * sizeof(real_t), cudaMemcpyDeviceToHost);
	
	// END SOLVE
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

	
	// FREE
	cudaFree(buffer);
	cusparseDestroyMatDescr(descr_M);
	cusparseDestroyMatDescr(descr_L);
	cusparseDestroyCsric02Info(info_M);
	cusparseDestroyCsrsv2Info(info_L);
	cusparseDestroyCsrsv2Info(info_Lt);
	cusparseDestroy(handle);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return (double) time / 1000;
}

double solve_positive_definite_cusolver_preview(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx, const bool host)
{
	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// SOLVE Ax = b

	cusolverSpHandle_t cusolver_handle = NULL;
	cusparseHandle_t cusparse_handle = NULL;
//	cudaStream_t stream = NULL;

	cusparseMatDescr_t descr = NULL;

	size_t size_iternal = 0;
	size_t size_chol = 0;

	void * buffer = NULL;

	int singularity;

	cusolverSpCreate(&cusolver_handle);
	cusparseCreate(&cusparse_handle);
/*
	cudaStreamCreate(&stream);
	cusolverSpSetStream(cusolver_handle, stream);
	cusparseSetStream(cusparse_handle, stream);
*/	
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	
	if(host)
	{
		csrcholInfoHost_t info;
		cusolverSpCreateCsrcholInfoHost(&info);

		cusolverSpXcsrcholAnalysisHost(cusolver_handle, m, nnz, descr, hA_col_ptrs, hA_row_indices, info);

		#ifdef SINGLE_P
			cusolverSpScsrcholBufferInfoHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, &size_iternal, &size_chol);
		#else
			cusolverSpDcsrcholBufferInfoHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, &size_iternal, &size_chol);
		#endif
	
		buffer = new char[size_chol];
		
		#ifdef SINGLE_P
			cusolverSpScsrcholFactorHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, buffer);
			cusolverSpScsrcholZeroPivotHost(cusolver_handle, info, 0, &singularity);
		#else
			cusolverSpDcsrcholFactorHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, buffer);
			cusolverSpDcsrcholZeroPivotHost(cusolver_handle, info, 0, &singularity);
		#endif
		assert(singularity == -1);

		// SOLVE
		cudaEventRecord(start, 0);
		
		#ifdef SINGLE_P
			cusolverSpScsrcholSolveHost(cusolver_handle, m, hb, hx, info, buffer);
		#else
			cusolverSpDcsrcholSolveHost(cusolver_handle, m, hb, hx, info, buffer);
		#endif
		
		// END SOLVE
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time, start, stop);

		// FREE
		delete [] (char*) buffer;
		cusolverSpDestroyCsrcholInfoHost(info);
	}
	else
	{
		cu_spAxb data(m, nnz, hA_values, hA_col_ptrs, hA_row_indices, hb, hx);

		csrcholInfo_t info;
		cusolverSpCreateCsrcholInfo(&info);

		cusolverSpXcsrcholAnalysis(cusolver_handle, m, nnz, descr, data.A_col_ptrs, data.A_row_indices, info);

		#ifdef SINGLE_P
			cusolverSpScsrcholBufferInfo(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, &size_iternal, &size_chol);
		#else
			cusolverSpDcsrcholBufferInfo(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, &size_iternal, &size_chol);
		#endif
	
		cudaMalloc(&buffer, size_chol);

		#ifdef SINGLE_P
			cusolverSpScsrcholFactor(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, buffer);
			cusolverSpScsrcholZeroPivot(cusolver_handle, info, 0, &singularity);
		#else
			cusolverSpDcsrcholFactor(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, buffer);
			cusolverSpDcsrcholZeroPivot(cusolver_handle, info, 0, &singularity);
		#endif

//		assert(singularity == -1);

		// SOLVE
		cudaEventRecord(start, 0);
		
		#ifdef SINGLE_P
			cusolverSpScsrcholSolve(cusolver_handle, m, data.b, data.x, info, buffer);
		#else
			cusolverSpDcsrcholSolve(cusolver_handle, m, data.b, data.x, info, buffer);
		#endif
		
		// END SOLVE
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time, start, stop);

		cudaMemcpy(hx, data.x, m * sizeof(real_t), cudaMemcpyDeviceToHost);

		// FREE
		cudaFree(buffer);
		cusolverSpDestroyCsrcholInfo(info);
	}

//	cudaStreamDestroy(stream);
	cusparseDestroyMatDescr(descr);
	cusparseDestroy(cusparse_handle);
	cusolverSpDestroy(cusolver_handle);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	return (double) time / 1000;
}


} // namespace gproshan

