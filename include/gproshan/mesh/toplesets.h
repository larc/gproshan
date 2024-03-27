#ifndef TOPLESETS_H
#define TOPLESETS_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class toplesets: public partitions
{
	std::vector<index_t> level;

	public:
		const size_t n_levels = 0;

	public:
		toplesets(const che * mesh, const std::vector<index_t> & sources, const index_t max_level = NIL)
		{
			reset(mesh, sources, max_level);
		}

		~toplesets()
		{
			delete [] sorted;
		}

		index_t operator [] (const index_t i) const
		{
			return level[i];
		}

		void reset(const che * mesh, const std::vector<index_t> & sources, const index_t max_level = NIL)
		{
			if(std::size(level) < mesh->n_vertices)
			{
				delete [] sorted;
				sorted = new index_t[mesh->n_vertices];
			}
			level.assign(mesh->n_vertices, NIL);

			index_t l = 0;
			index_t n = 0;
			for(index_t s: sources)
			{
				sorted[n++] = s;
				level[s] = l;
			}

			splits.clear();
			splits.push_back(0);
			for(index_t i = 0; i < n; ++i)
			{
				const index_t v = sorted[i];

				if(level[v] > l)
				{
					if(++l > max_level) break;
					splits.push_back(i);
				}

				for(index_t u: mesh->link(v))
					if(level[u] == NIL)
					{
						level[u] = level[v] + 1;
						sorted[n++] = u;
					}
			}
			splits.push_back(n);

			che::rw(n_levels) = std::size(splits) - 1;
		}
};


} // namespace gproshan

#endif // TOPLESETS_H

