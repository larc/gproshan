#ifndef UTIL_H
#define UTIL_H

#include <gproshan/include.h>

#include <cassert>
#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


class partitions
{
	struct part;

	std::vector<index_t> splits;
	index_t * sorted = nullptr;

	public:
		part operator () (const index_t & i) const;
};

struct partitions::part
{
	struct iterator;

	index_t _begin = 0;
	index_t _end = 0;
	const index_t * sorted = nullptr;

	iterator begin() const;
	iterator end() const;
};

struct partitions::part::iterator
{
	index_t i = 0;
	const index_t * sorted = nullptr;

	iterator & operator ++ ();
	bool operator != (const iterator & it) const;
	const index_t & operator * ();
};


void copy_real_t_array(float * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(float * destination, const double * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const double * source, const size_t & n_elem);


} // namespace gproshan

#endif // UTIL_H

