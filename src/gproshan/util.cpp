#include <gproshan/util.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan {


partitions::partitions(index_t * s): sorted(s)
{
	splits.push_back(0);
}

void partitions::add(const index_t & size)
{
	return splits.push_back(size + splits.back());
}

size_t partitions::size(const index_t & i) const
{
	return splits[i + 1] - splits[i];
}

partitions::part partitions::operator () (const index_t & i) const
{
	assert(i > 0 && i < size(splits));
	return {splits[i], splits[i + 1], sorted};
}


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

