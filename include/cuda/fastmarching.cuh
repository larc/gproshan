#ifndef FASTMARCHING_CUH
#define FASTMARCHING_CUH

#include "che.cuh"

distance_t * cpu_fastmarching(CHE * mesh, index_t * source, length_t source_size, vector<index_t> & limites, index_t * sorted, index_t * clusters = NULL);

distance_t * cuda_fastmarching(CHE * h_che, CHE * d_che, index_t * source, length_t source_size, vector<index_t> & limites, index_t * h_sorted, index_t * h_clusters = NULL, distance_t * real_dist = NULL);

distance_t * cuda_fastmarching(CHE * h_che, CHE * d_che, index_t * source, length_t source_size, index_t iter);

distance_t * parallel_fastmarching(che * mesh, index_t * source, length_t source_size, float & time_g, index_t iter, bool cuda, bool normalize, index_t * clusters = NULL, bool GPU = true, distance_t * real_dist = NULL);

#endif

