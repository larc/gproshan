#include "fastmarching.cuh"

#include <cstdio>
#include <fstream>
#include <cublas_v2.h>

#define NT 32
#define NB(x) (x + NT - 1) / NT

inline index_t iterations(const vector<index_t> & limites)
{
	index_t max_i = 1;
	for(index_t i = 2; i < limites.size(); i++)
		if(limites[i] - limites[i - 1] >= limites[max_i] - limites[max_i - 1])
			max_i = i;
	
//	return max_i + limites.size() + 1;
	return limites.size() << 1;
}

inline index_t start_v(const index_t & i, const vector<index_t> & limites)
{
	return limites[i >> 1];
}

inline index_t end_v(const index_t & i, const vector<index_t> & limites)
{
	index_t di = i - (i >> 1) - 2;
	di = i - (di >> 1);
	return limites[i < limites.size() ? i : limites.size() - 1];
}

__host__ __device__
distance_t update_step(CHE * mesh, const distance_t * dist, const index_t & he)
{
	index_t x[3];
	x[0] = mesh->VT[cu_next(he)];
	x[1] = mesh->VT[cu_prev(he)];
	x[2] = mesh->VT[he];

	vertex_cu X[2];
	X[0] = mesh->GT[x[0]] - mesh->GT[x[2]];
	X[1] = mesh->GT[x[1]] - mesh->GT[x[2]];

	distance_t t[2];
	t[0] = dist[x[0]];
	t[1] = dist[x[1]];

	distance_t q[2][2];
	q[0][0] = (X[0], X[0]);
	q[0][1] = (X[0], X[1]);
	q[1][0] = (X[1], X[0]);
	q[1][1] = (X[1], X[1]);
	
	distance_t det = q[0][0] * q[1][1] - q[0][1] * q[1][0];
	distance_t Q[2][2];
	Q[0][0] = q[1][1] / det;
	Q[0][1] = -q[0][1] / det;
	Q[1][0] = -q[1][0] / det;
	Q[1][1] = q[0][0] / det;


	distance_t delta = t[0] * (Q[0][0] + Q[1][0]) + t[1] * (Q[0][1] + Q[1][1]);

	distance_t dis = delta * delta - (Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1]) * (t[0]*t[0]*Q[0][0] + t[0]*t[1]*(Q[1][0] + Q[0][1]) + t[1]*t[1]*Q[1][1] - 1);
	
	distance_t p;

	if(dis >= 0)
	{
		p = delta + sqrt(dis);
		p /= Q[0][0] + Q[0][1] + Q[1][0] + Q[1][1];
	}

	distance_t tp[2];
	tp[0] = t[0] - p;
	tp[1] = t[1] - p;

	vertex_cu n(tp[0] * (X[0][0]*Q[0][0] + X[1][0]*Q[1][0]) + tp[1] * (X[0][0]*Q[0][1] + X[1][0]*Q[1][1]),
			 tp[0] * (X[0][1]*Q[0][0] + X[1][1]*Q[1][0]) + tp[1] * (X[0][1]*Q[0][1] + X[1][1]*Q[1][1]),
			 tp[0] * (X[0][2]*Q[0][0] + X[1][2]*Q[1][0]) + tp[1] * (X[0][2]*Q[0][1] + X[1][2]*Q[1][1]) );

	distance_t cond[2];
	cond[0] = (X[0] , n);
	cond[1] = (X[1] , n);

	distance_t c[2];
	c[0] = cond[0] * Q[0][0] + cond[1] * Q[0][1];
	c[1] = cond[0] * Q[1][0] + cond[1] * Q[1][1];

	if(t[0] == INFINITY || t[1] == INFINITY || dis < 0 || c[0] >= 0 || c[1] >= 0)
	{
		distance_t dp[2];
		dp[0] = dist[x[0]] + *X[0];
		dp[1] = dist[x[1]] + *X[1];

		p = dp[dp[1] < dp[0]];
	}

	return p;
}

distance_t * cpu_fastmarching(CHE * mesh, index_t * source, length_t source_size, vector<index_t> & limites, index_t * sorted, index_t * clusters)
{
	distance_t * dist[2];
	dist[0] = new distance_t[mesh->n_vertices];
	dist[1] = new distance_t[mesh->n_vertices];

	for(index_t v = 0; v < mesh->n_vertices; v++)
		dist[0][v] = dist[1][v] = INFINITY;

	for(index_t i = 0; i < source_size; i++)
	{
		dist[0][source[i]] = dist[1][source[i]] = 0;
		if(clusters) clusters[source[i]] = i + 1;
	}

	index_t v, d = 1;
	index_t iter = iterations(limites);
	index_t start, end;

	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limites);
		end = end_v(i, limites);

		#pragma omp parallel for private(v)
		for(index_t vi = start; vi < end; vi++)
		{
			v = sorted[vi];

			dist[!d][v] = dist[d][v];
			
			distance_t p;
		
			cu_for_star(he, mesh, v)
			{		
				p = update_step(mesh, dist[d], he);
		
				if(p < dist[!d][v])
				{
					dist[!d][v] = p;
					if(clusters)
						clusters[v] = clusters[mesh->VT[cu_prev(he)]] != NIL ? clusters[mesh->VT[cu_prev(he)]] : clusters[mesh->VT[cu_next(he)]];
				}
			}
		}

		d = !d;
	}

	delete [] dist[d];
	return dist[!d];
}

__global__
void fastmarching_relax(CHE * d_mesh, distance_t * d_dist, distance_t * d_prevdist, index_t end, index_t * d_clusters = NULL, index_t * d_prevclusters = NULL, index_t start = 0, index_t * sorted = NULL, index_t iter = NIL, index_t * arrive = NULL, index_t * prevarrive = NULL)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;
	
	if(v < end)
	{
		v = sorted ? sorted[v] : v;

		if(v < d_mesh->n_vertices)
		{
			d_dist[v] = d_prevdist[v];
			if(d_clusters) d_clusters[v] = d_prevclusters[v];
			if(arrive) arrive[v] = prevarrive[v];

			distance_t d;

			cu_for_star(he, d_mesh, v)
			{
				d = update_step(d_mesh, d_prevdist, he);
		
				if(arrive && arrive[v] == NIL && 
					prevarrive[d_mesh->VT[cu_next(he)]] != NIL && prevarrive[d_mesh->VT[cu_prev(he)]] != NIL &&
					prevarrive[d_mesh->VT[cu_next(he)]] < iter && prevarrive[d_mesh->VT[cu_prev(he)]] < iter)
					arrive[v] = iter;

				if(d < d_dist[v])
				{
					d_dist[v] = d;
					if(d_clusters)
						d_clusters[v] = d_prevdist[d_mesh->VT[cu_prev(he)]] < d_prevdist[d_mesh->VT[cu_next(he)]] ? d_prevclusters[d_mesh->VT[cu_prev(he)]] : d_prevclusters[d_mesh->VT[cu_next(he)]];
				}
			}
		}
	}
}

__global__
void fastmarching_relax(CHE * mesh, distance_t * new_dist, distance_t * prev_dist, index_t end, index_t * sorted, index_t start = 0)
{
	index_t v = blockDim.x * blockIdx.x + threadIdx.x + start;
	
	if(v < end)
	{
		v = sorted ? sorted[v] : v;

		if(v < mesh->n_vertices)
		{
			new_dist[v] = prev_dist[v];

			distance_t d;
			cu_for_star(he, mesh, v)
			{
				d = update_step(mesh, prev_dist, he);
				if(d < new_dist[v]) new_dist[v] = d;
			}
		}
	}
}

index_t cuda_fastmarching(CHE * d_mesh, const length_t & n, distance_t * h_dist, distance_t ** d_dist, const index_t *const source, const length_t & source_size, const vector<index_t> & limites, const index_t *const h_sorted, index_t * d_sorted)
{
	for(index_t i = 0; i < n; i++)
		h_dist[i] = INFINITY;

	for(index_t i = 0; i < source_size; i++)
		h_dist[source[i]] = 0;
	
	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sorted, h_sorted, sizeof(index_t) * n, cudaMemcpyHostToDevice);

	index_t d = 1;
	index_t start, end;
	
	index_t iter = iterations(limites);
	for(index_t i = 2; i < iter; i++)
	{
		start = start_v(i, limites);
		end = end_v(i, limites);
		fastmarching_relax<<<NB(end - start), NT>>>(d_mesh, d_dist[!d], d_dist[d], end, d_sorted, start);
		cudaDeviceSynchronize();
		d = !d;
	}
	
	return d;
}

distance_t * cuda_fastmarching(CHE * d_mesh, const length_t & n, const index_t *const source, const length_t & source_size, const vector<index_t> & limites, const index_t *const h_sorted, const bool & device)
{
	distance_t * h_dist = new distance_t[n];
	
	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * n);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * n);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * n);
	
	index_t d = cuda_fastmarching(d_mesh, n, h_dist, d_dist, source, source_size, limites, h_sorted, d_sorted);
	
	if(device)
	{
		cudaFree(d_dist[!d]);
		delete [] h_dist;
	}
	else
	{
		cudaMemcpy(h_dist, d_dist[d], sizeof(distance_t) * n, cudaMemcpyDeviceToHost);
		cudaFree(d_dist[0]);
		cudaFree(d_dist[1]);
	}
	cudaFree(d_sorted);

	return device ? d_dist[d] : h_dist;
}

distance_t * cuda_fastmarching(CHE * h_mesh, CHE * d_mesh, index_t * source, length_t source_size, vector<index_t> & limites, index_t * h_sorted, index_t * h_clusters, distance_t * real_dist)
{
	FILE * f_iter_error;
	if(real_dist) f_iter_error = fopen("iter_error", "w");


	distance_t * h_dist = new distance_t[h_mesh->n_vertices];

	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * h_mesh->n_vertices);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * h_mesh->n_vertices);

	index_t * d_clusters[2] = {NULL, NULL};
	if(h_clusters) cudaMalloc(&d_clusters[0], sizeof(index_t) * h_mesh->n_vertices);	
	if(h_clusters) cudaMalloc(&d_clusters[1], sizeof(index_t) * h_mesh->n_vertices);	

	index_t * h_arrive = NULL; //new index_t[h_mesh->n_vertices];
	index_t * d_arrive[2] = {NULL, NULL};
	if(h_arrive) cudaMalloc(&d_arrive[0], sizeof(index_t) * h_mesh->n_vertices);	
	if(h_arrive) cudaMalloc(&d_arrive[1], sizeof(index_t) * h_mesh->n_vertices);	
	
	
	if(h_arrive)
	{
		memset(h_arrive, 255, sizeof(index_t) * h_mesh->n_vertices);
		for(index_t i = 0; i < limites[1]; i++)
			h_arrive[h_sorted[i]] = 0;
		for(index_t i = limites[1]; i < limites[2]; i++)
			h_arrive[h_sorted[i]] = 1;
	}

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * h_mesh->n_vertices);

	for(index_t i = 0; i < h_mesh->n_vertices; i++)
		h_dist[i] = INFINITY;

	for(index_t i = 0; i < source_size; i++)
	{
		h_dist[source[i]] = 0;
		if(h_clusters) h_clusters[source[i]] = i + 1;
	}

	cudaMemcpy(d_dist[0], h_dist, sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dist[1], h_dist, sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_sorted, h_sorted, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_clusters) cudaMemcpy(d_clusters[0], h_clusters, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_clusters) cudaMemcpy(d_clusters[1], h_clusters, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_arrive) cudaMemcpy(d_arrive[0], h_arrive, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);
	if(h_arrive) cudaMemcpy(d_arrive[1], h_arrive, sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyHostToDevice);

	index_t d = 1;	
	index_t start, end, i = 2;

//	bool stop = false;

	index_t iter = real_dist ? limites.size() << 1 : iterations(limites);

	for(i = 2; i < iter; i++)
//	while(i && !stop)
	{
		start = start_v(i, limites);
		end = end_v(i, limites);

		fastmarching_relax<<<NB(end - start), NT>>>(d_mesh, d_dist[!d], d_dist[d], end, d_clusters[!d], d_clusters[d], start, d_sorted, i, d_arrive[!d], d_arrive[d]);
		
		if(real_dist && i >= limites.size())
		{
			cudaMemcpy(h_dist, d_dist[!d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
			distance_t error = 0;
			
			#pragma omp parallel for reduction(+: error)
			for(index_t v = 1; v < h_mesh->n_vertices; v++)
				error += abs(h_dist[v] - real_dist[v]) / real_dist[v];
			
			error /= h_mesh->n_vertices - source_size;
			if(error < INFINITY) fprintf(f_iter_error, "%d %.10f\n", i, error);
		}
		else
			cudaDeviceSynchronize();

		d = !d;
		
	/*
		if(h_arrive)
		{
			cudaMemcpy(h_arrive, d_arrive[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
			#pragma omp parallel for
			for(index_t s = limites[limites.size() - 2]; s < h_mesh->n_vertices; s++)
				if(h_arrive[h_sorted[s]] != NIL)
				{
					#pragma omp atomic
					stop |= true;
				}
		}
	
		i++;
	*/
	}

	if(real_dist) fclose(f_iter_error);

	cudaMemcpy(h_dist, d_dist[d], sizeof(distance_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
	if(h_clusters) cudaMemcpy(h_clusters, d_clusters[d], sizeof(index_t) * h_mesh->n_vertices, cudaMemcpyDeviceToHost);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);
	if(h_clusters)
	{
		cudaFree(d_clusters[0]);
		cudaFree(d_clusters[1]);
	}
	return h_dist;
}

inline index_t farthest(distance_t * d, size_t n)
{
	index_t f = 0;
	
	#pragma omp parallel for
	for(index_t v = 0; v < n; v++)
		#pragma omp critical
		if(d[v] < INFINITY && d[f] < d[v])
			f = v;
	
	return f;
}

distance_t farthest_point_sampling_gpu(vector<index_t> & points, float & time, che * mesh, size_t n, distance_t radio)
{
	CHE * h_mesh;
	CHE * dd_mesh;
	CHE * d_mesh;

	h_mesh = new CHE(mesh);

	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	ofstream os((PATH_TEST + "fastmarching/" + mesh->name() + ".fps").c_str());
	
	size_t n_v = mesh->n_vertices();

	index_t * rings = new index_t[n_v];
	index_t * h_sorted = new index_t[n_v];
	
	distance_t * h_dist = new distance_t[n_v];
	
	distance_t * d_dist[2];
	cudaMalloc(&d_dist[0], sizeof(distance_t) * n_v);
	cudaMalloc(&d_dist[1], sizeof(distance_t) * n_v);

	index_t * d_sorted;
	cudaMalloc(&d_sorted, sizeof(index_t) * n_v);

	// ---------------------------------------------------------------------------------------------
	cudaEvent_t start;
	cudaEvent_t stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cublasHandle_t handle;
	cublasCreate(&handle);

	if(n >= mesh->n_vertices())
		n = mesh->n_vertices() / 2;

	n -= points.size();
	points.reserve(n);
	
	time = 0;
	float time_aux;

	index_t d;
	int f;
	distance_t max_dis = INFINITY;
	while(n-- && max_dis > radio)
	{
		cudaEventRecord(start, 0);
		
		vector<index_t> limites;
		mesh->sort_by_rings(rings, h_sorted, limites, points);
		
		d = cuda_fastmarching(d_mesh, n_v, h_dist, d_dist, points.data(), points.size(), limites, h_sorted, d_sorted);
				
		// 1 indexing T_T
		#ifdef SINGLE_P
		cublasIsamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#else
		cublasIdamax(handle, mesh->n_vertices(), d_dist[d], 1, &f);
		#endif
		
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&time_aux, start, stop);

		time_aux /= 1000;
		time += time_aux;
		
		os << points.size() << " " << time_aux << endl;
		
		if(radio > 0 || !n)
			cudaMemcpy(&max_dis, d_dist[d] + f - 1, sizeof(distance_t), cudaMemcpyDeviceToHost);
		points.push_back(f - 1);
	}
	

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cublasDestroy(handle);
	// ---------------------------------------------------------------------------------------------
	
	os.close();

	delete [] rings;
	delete [] h_sorted;
	delete [] h_dist;

	cuda_free_CHE(dd_mesh, d_mesh);
	cudaFree(d_dist[0]);
	cudaFree(d_dist[1]);
	cudaFree(d_sorted);

	return max_dis;
}

distance_t * parallel_fastmarching(che * mesh, index_t * source, length_t source_size, float & time_g, vector<index_t> & limites, index_t * sorted, bool normalize, index_t * clusters, bool GPU, distance_t * real_dist)
{
	CHE * h_mesh;
	CHE * dd_mesh;
	CHE * d_mesh;

	h_mesh = new CHE(mesh);

	cuda_create_CHE(h_mesh, dd_mesh, d_mesh);

	cudaEvent_t start;
	cudaEvent_t stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	distance_t * distances;
	cudaEventRecord(start, 0);
	
	if(GPU)
		distances = cuda_fastmarching(h_mesh, d_mesh, source, source_size, limites, sorted, clusters, real_dist);
//		distances = cuda_fastmarching(d_mesh, mesh->n_vertices(), source, source_size, limites, sorted, false);
	else
		distances = cpu_fastmarching(h_mesh, source, source_size, limites, sorted, clusters);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_g, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	time_g /= 1000;

	if(normalize)
	{
		distance_t max_d = 0;
		
		#pragma omp parallel for reduction(max: max_d)
		for(index_t v = 0; v < h_mesh->n_vertices; v++)
			if(distances[v] < INFINITY)
				max_d = max(distances[v], max_d);
		
		#pragma omp parallel for shared(max_d)
		for(index_t v = 0; v < h_mesh->n_vertices; v++)
			distances[v] /= max_d;
	}


	//free_CHE(h_mesh);
	cuda_free_CHE(dd_mesh, d_mesh);

	return distances;
}

