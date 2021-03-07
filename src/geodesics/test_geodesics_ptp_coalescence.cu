#include "geodesics/test_geodesics_ptp_coalescence.cuh"

#include "geodesics/geodesics_ptp_coalescence.cuh"
#include "geodesics/geodesics_ptp.h"
#include "geodesics/test_geodesics_ptp.h"

#include "mesh/che_off.h"

#include <fstream>
#include <cublas_v2.h>

#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


vector<pair<index_t, real_t> > iter_error_parallel_toplesets_propagation_coalescence_gpu(che * mesh, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * sorted_index, const real_t * exact_dist, double & time_ptp)
{
	// sort data by levels, must be improve the coalescence

	vertex * V = new vertex[mesh->n_vertices];
	index_t * F = new index_t[mesh->n_faces * che::mtrig];
	index_t * inv = new index_t[mesh->n_vertices];
	real_t * exact_dist_sorted = new real_t[mesh->n_vertices];

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		V[i] = mesh->gt(sorted_index[i]);
		inv[sorted_index[i]] = i;
		exact_dist_sorted[i] = exact_dist[sorted_index[i]];
	}

	#pragma omp parallel for
	for(index_t he = 0; he < mesh->n_half_edges; ++he)
		F[he] = inv[mesh->vt(he)];

	mesh = new che(V, mesh->n_vertices, F, mesh->n_faces);

	delete [] V;
	delete [] F;

	// ------------------------------------------------------

	cudaDeviceReset();

	float time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// BEGIN PTP

	CHE * h_mesh = new CHE(mesh);
	CHE * dd_mesh, * d_mesh;
	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	real_t * h_dist = new real_t[h_mesh->n_vertices];

	real_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(real_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * h_mesh->n_vertices);

	real_t * d_error;
	cudaMalloc(&d_error, sizeof(real_t) * h_mesh->n_vertices);

	vector<pair<index_t, real_t> > iter_error = iter_error_run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, sources, limits, inv, exact_dist_sorted, d_error);

	delete [] h_dist;
	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cuda_free_CHE(dd_mesh, d_mesh);

	// END PTP

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time_ptp = time / 1000;

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	delete mesh;
	delete [] inv;

	return iter_error;
}

/// Return an array of time in seconds.
double * times_farthest_point_sampling_ptp_coalescence_gpu(che * mesh, vector<index_t> & samples, size_t n, real_t radio)
{
	cudaDeviceReset();

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// BEGIN FPS PTP

	vertex * V = new vertex[mesh->n_vertices];
	index_t * F = new index_t[mesh->n_faces * che::mtrig];
	index_t * inv = new index_t[mesh->n_vertices];


	real_t * h_dist = new real_t[mesh->n_vertices];

	real_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(real_t) * mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(real_t) * mesh->n_vertices);

	real_t * d_error;
	cudaMalloc(&d_error, sizeof(real_t) * mesh->n_vertices);

	vector<index_t> limits;
	index_t * toplesets = new index_t[mesh->n_vertices];
	index_t * sorted_index = new index_t[mesh->n_vertices];

	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= mesh->n_vertices) n = mesh->n_vertices >> 1;

	double * times = new double[n + 1];

	n -= samples.size();
	samples.reserve(n);

	float time_fps;
	index_t d;
	int f;
	real_t max_dist = INFINITY;
	while(n-- && max_dist > radio)
	{
		cudaEventRecord(start, 0);

		limits.clear();
		mesh->compute_toplesets(toplesets, sorted_index, limits, samples);

		// sort data by levels, must be improve the coalescence

		#pragma omp parallel for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
		{
			V[i] = mesh->gt(sorted_index[i]);
			inv[sorted_index[i]] = i;
		}

		#pragma omp parallel for
		for(index_t he = 0; he < mesh->n_half_edges; ++he)
			F[he] = inv[mesh->vt(he)];

		che * tmp_mesh = new che(V, mesh->n_vertices, F, mesh->n_faces);

		CHE * h_mesh = new CHE(tmp_mesh);
		CHE * dd_mesh, * d_mesh;
		cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

		// exec algorithm
		d = run_ptp_coalescence_gpu(d_mesh, h_mesh->n_vertices, h_dist, d_dist, samples, {limits, inv}, d_error);

		// free memory
		cuda_free_CHE(dd_mesh, d_mesh);
		delete tmp_mesh;

		// 1 indexing
		#ifdef SINGLE_P
			cublasIsamax(handle, mesh->n_vertices, d_dist[d], 1, &f);
		#else
			cublasIdamax(handle, mesh->n_vertices, d_dist[d], 1, &f);
		#endif

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time_fps, start, stop);

		times[samples.size()] = time_fps / 1000;

		if(radio > 0 || !n)
			cudaMemcpy(&max_dist, d_dist[d] + f - 1, sizeof(real_t), cudaMemcpyDeviceToHost);

		samples.push_back(sorted_index[f - 1]);
	}

	cublasDestroy(handle);

	delete [] V;
	delete [] F;
	delete [] inv;
	delete [] h_dist;
	delete [] toplesets;
	delete [] sorted_index;

	cudaFree(d_error);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);

	// END FPS PTP

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return times;
}

vector<pair<index_t, real_t> > iter_error_run_ptp_coalescence_gpu(CHE * d_mesh, const index_t & n_vertices, real_t * h_dist, real_t ** d_dist, const vector<index_t> & sources, const vector<index_t> & limits, const index_t * inv, const real_t * exact_dist, real_t * d_error)
{
	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices; ++v)
		h_dist[v] = INFINITY;

	for(index_t i = 0; i < sources.size(); ++i)
		h_dist[inv[sources[i]]] = 0;

	cudaMemcpy(d_dist[0], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(real_t) * n_vertices, cudaMemcpyHostToDevice);

	vector<pair<index_t, real_t> > iter_error;
	iter_error.reserve(limits.size());

	ofstream os("band");

	index_t d = 0;
	index_t start, end, n_cond;
	index_t i = 1, j = 2;
	index_t n_iter = 0;

	while(i < j)
	{
		++n_iter;
		start = limits[i];
		end = limits[j];
		n_cond = limits[i + 1] - start;

		relax_ptp_coalescence <<< NB(end - start), NT >>> (d_mesh, d_dist[!d], d_dist[d], end, start);
		// print band info
		os << n_iter << " " << i << " " << j << " " << end - start << endl;

		// begin calculating iteration error
		cudaMemcpy(h_dist, d_dist[!d], sizeof(real_t) * n_vertices, cudaMemcpyDeviceToHost);
		if(j == limits.size() - 1)
			iter_error.push_back(make_pair(n_iter, compute_error(h_dist, exact_dist, n_vertices, sources.size())));
		// end

		relative_error <<< NB(n_cond), NT >>> (d_error, d_dist[!d], d_dist[d], start, start + n_cond);
		cudaDeviceSynchronize();

		if(n_cond == thrust::count_if(thrust::device, d_error + start, d_error + start + n_cond, is_ok()))
			++i;
		if(j < limits.size() - 1) ++j;

		d = !d;
	}

	os.close();

	return iter_error;
}


} // namespace gproshan

