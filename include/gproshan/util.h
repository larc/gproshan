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
		partitions(index_t * s = nullptr);
		void add(const index_t size);
		size_t size(const index_t i) const;
		part operator () (const index_t i) const;
		operator index_t * () const;
		operator index_t *& ();
};

struct partitions::part
{
	struct iterator
	{
		index_t i = 0;
		const index_t * sorted = nullptr;

		__host_device__
		iterator & operator ++ ()
		{
			++i;
			return *this;
		}

		__host_device__
		bool operator != (const iterator & it) const
		{
			return i != it.i;
		}

		__host_device__
		index_t operator * () const
		{
			return sorted ? sorted[i] : i;
		}
	};


	index_t _begin = 0;
	index_t _end = 0;
	const index_t * sorted = nullptr;


	__host_device__
	iterator begin() const
	{
		return {_begin, sorted};
	}

	__host_device__
	iterator end() const
	{
		return {_end, sorted};
	}
};


void copy_real_t_array(float * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(float * destination, const double * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const double * source, const size_t & n_elem);


} // namespace gproshan

#endif // UTIL_H

