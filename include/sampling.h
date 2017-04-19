/*! \file */

#ifndef SAMPLING_H
#define SAMPLING_H

#include "geodesics.h"

#include <vector>
#include <string>

distance_t * parallel_fastmarching(che * mesh, index_t * source, length_t source_size, float & time_g, vector<index_t> & limites, index_t * sorted, bool normalize = false, index_t * clusters = NULL, bool GPU = true, distance_t * real_dist = NULL);

distance_t * parallel_fastmarching(const char * filename, int * source, length_t source_size, float & time, int iter, bool cuda = 1, bool normalize = true);

distance_t parallel_farthest_point_sampling(vector<index_t> & points, che * shape, size_t n, distance_t radio = 0);

index_t ** sampling_shape(vector<index_t> & points, size_t *& sizes, vertex *& normals, che * shape, size_t n_points, distance_t radio);

void sampling_shape(int nargs, char ** args);

void sampling_shape(const char * name);

bool load_sampling(vector<index_t> & points, distance_t & radio, che * mesh, size_t M);

distance_t farthest_point_sampling_gpu(vector<index_t> & points, float & time, che * mesh, size_t n, distance_t radio = 0);

#endif // SAMPLING_H

