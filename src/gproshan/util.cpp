#include <gproshan/util.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan {


partitions::partitions(index_t * s): sorted(s)
{
	splits.push_back(0);
}

void partitions::add(const index_t size)
{
	return splits.push_back(size + splits.back());
}

size_t partitions::size(const index_t i) const
{
	return splits[i + 1] - splits[i];
}

partitions::part partitions::operator () (const index_t i) const
{
	assert(i > 0 && i < std::size(splits));
	return {splits[i], splits[i + 1], sorted};
}

partitions::operator index_t * () const
{
	return sorted;
}

partitions::operator index_t *& ()
{
	return sorted;
}


} // namespace gproshan

