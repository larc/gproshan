#include "util.h"

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan {


void copy_real_t_array(float * destination, const float * source, const size_t & n_elem)
{
	memcpy(destination, source, n_elem * sizeof(float));
}

void copy_real_t_array(float * destination, const double * source, const size_t & n_elem)
{
	#pragma omp parallel for
	for(index_t i = 0; i < n_elem; ++i)
		destination[i] = source[i];
}

void copy_real_t_array(double * destination, const float * source, const size_t & n_elem)
{
	#pragma omp parallel for
	for(index_t i = 0; i < n_elem; ++i)
		destination[i] = source[i];
}

void copy_real_t_array(double * destination, const double * source, const size_t & n_elem)
{
	memcpy(destination, source, n_elem * sizeof(double));
}


} // namespace gproshan

