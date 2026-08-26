#include <gproshan/mesh/toplesets.h>


// geometry processing and shape analysis framework
namespace gproshan {


toplesets::toplesets(const che * mesh, const std::vector<index_t> & sources, const index_t max_level)
{
	reset(mesh, sources, max_level);
}

index_t toplesets::operator [] (const index_t i) const
{
	return level[i];
}

toplesets::operator const std::vector<index_t> & () const
{
	return tsorted;
}

void toplesets::reset(const che * mesh, const std::vector<index_t> & sources, const index_t max_level)
{
	level.assign(mesh->n_vertices, NIL);

	tsorted.clear();
	tsorted.reserve(mesh->n_vertices);

	index_t l = 0;
	for(index_t s: sources)
	{
		level[s] = l;
		tsorted.push_back(s);
	}

	splits.clear();
	splits.push_back(0);
	for(index_t i = 0; i < std::size(tsorted); ++i)
	{
		const index_t v = tsorted[i];

		if(level[v] > l)
		{
			if(++l > max_level) break;
			splits.push_back(i);
		}

		for(index_t u: mesh->link(v))
			if(level[u] == NIL)
			{
				level[u] = level[v] + 1;
				tsorted.push_back(u);
			}
	}
	splits.push_back(std::size(tsorted));

	che::rw(n_levels) = std::size(splits) - 1;
	che::rw(n_vertices) = mesh->n_vertices;
	sorted = tsorted.data();
}


} // namespace gproshan

