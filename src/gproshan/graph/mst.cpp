#include <gproshan/graph/mst.h>


// geometry processing and shape analysis framework
namespace gproshan {


union_find::union_find(unsigned n) : n(n)
{
	sets = new unsigned[n];

	for(unsigned i = 0; i < n; ++i)
		sets[i] = i;
}

union_find::~union_find()
{
	delete sets;
}

unsigned union_find::find(const unsigned x)
{
	return sets[x] == x ? x : sets[x] = find(sets[x]);
}

bool union_find::merge(unsigned x, unsigned y)
{
	x = find(x);
	y = find(y);

	if(x == y) return false;

	sets[x] = y;

	return true;
}


} // namespace gproshan

