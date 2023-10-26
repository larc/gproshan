#include <gproshan/util.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan {


partitions::part partitions::operator () (const index_t & i) const
{
	assert(i > 0 && i < splits.size());
	return {splits[i - 1], splits[i], sorted};
}

partitions::part::iterator partitions::part::begin() const
{
	return {_begin, sorted};
}

partitions::part::iterator partitions::part::end() const
{
	return {_end, sorted};
}

partitions::part::iterator & partitions::part::iterator::operator ++ ()
{
	++i;
	return *this;
}

bool partitions::part::iterator::operator != (const iterator & it) const
{
	return i != it.i;
}

const index_t & partitions::part::iterator::operator * ()
{
	return sorted ? sorted[i] : i;
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

