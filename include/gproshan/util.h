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

	public:
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


template <class T>
T normalize(T * data, const size_t n_elem)
{
	T max = 0;

	#pragma omp parallel for reduction(max: max)
	for(index_t i = 0; i < n_elem; ++i)
		max = std::max(max, data[i]);

	if(max <= 1) return 1;

	#pragma omp parallel for
	for(index_t i = 0; i < n_elem; ++i)
		data[i] /= max;

	return max;
}


} // namespace gproshan

#endif // UTIL_H

