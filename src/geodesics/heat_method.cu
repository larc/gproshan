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

	cu_spAxb(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb)
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
		#ifdef GPROSHAN_FLOAT
			cusolverSpScsrlsvcholHost(handle_cusolver, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, hb, 0, 0, hx, &singularity);
		#else
			cusolverSpDcsrlsvcholHost(handle_cusolver, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, hb, 0, 0, hx, &singularity);
		#endif
	}
	else
	{
		// allocate A, x, b into device
		cu_spAxb data(m, nnz, hA_values, hA_col_ptrs, hA_row_indices, hb);

		cusolverStatus_t status;
		#ifdef GPROSHAN_FLOAT
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

double solve_positive_definite_cusolver_preview(const int m, const int nnz, const real_t * hA_values, const int * hA_col_ptrs, const int * hA_row_indices, const real_t * hb, real_t * hx, const bool host)
{
	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// SOLVE Ax = b

	cusolverSpHandle_t cusolver_handle = nullptr;
	cusparseHandle_t cusparse_handle = nullptr;
//	cudaStream_t stream = nullptr;

	cusparseMatDescr_t descr = nullptr;

	size_t size_iternal = 0;
	size_t size_chol = 0;

	void * buffer = nullptr;

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

		#ifdef GPROSHAN_FLOAT
			cusolverSpScsrcholBufferInfoHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, &size_iternal, &size_chol);
		#else
			cusolverSpDcsrcholBufferInfoHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, &size_iternal, &size_chol);
		#endif

		buffer = new char[size_chol];

		#ifdef GPROSHAN_FLOAT
			cusolverSpScsrcholFactorHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, buffer);
			cusolverSpScsrcholZeroPivotHost(cusolver_handle, info, 0, &singularity);
		#else
			cusolverSpDcsrcholFactorHost(cusolver_handle, m, nnz, descr, hA_values, hA_col_ptrs, hA_row_indices, info, buffer);
			cusolverSpDcsrcholZeroPivotHost(cusolver_handle, info, 0, &singularity);
		#endif
		assert(singularity == -1);

		// SOLVE
		cudaEventRecord(start, 0);

		#ifdef GPROSHAN_FLOAT
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
		cu_spAxb data(m, nnz, hA_values, hA_col_ptrs, hA_row_indices, hb);

		csrcholInfo_t info;
		cusolverSpCreateCsrcholInfo(&info);

		cusolverSpXcsrcholAnalysis(cusolver_handle, m, nnz, descr, data.A_col_ptrs, data.A_row_indices, info);

		#ifdef GPROSHAN_FLOAT
			cusolverSpScsrcholBufferInfo(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, &size_iternal, &size_chol);
		#else
			cusolverSpDcsrcholBufferInfo(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, &size_iternal, &size_chol);
		#endif

		cudaMalloc(&buffer, size_chol);

		#ifdef GPROSHAN_FLOAT
			cusolverSpScsrcholFactor(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, buffer);
			cusolverSpScsrcholZeroPivot(cusolver_handle, info, 0, &singularity);
		#else
			cusolverSpDcsrcholFactor(cusolver_handle, m, nnz, descr, data.A_values, data.A_col_ptrs, data.A_row_indices, info, buffer);
			cusolverSpDcsrcholZeroPivot(cusolver_handle, info, 0, &singularity);
		#endif

//		assert(singularity == -1);

		// SOLVE
		cudaEventRecord(start, 0);

		#ifdef GPROSHAN_FLOAT
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

