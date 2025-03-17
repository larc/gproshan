#ifndef TOPLESETS_H
#define TOPLESETS_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class toplesets: public partitions
{
	std::vector<index_t> level;
	std::vector<index_t> tsorted;

	public:
		const size_t n_levels = 0;

	public:
		toplesets(const che * mesh, const std::vector<index_t> & sources, const index_t max_level = NIL);

		index_t operator [] (const index_t i) const;
		operator const std::vector<index_t> & () const;

		void reset(const che * mesh, const std::vector<index_t> & sources, const index_t max_level = NIL);
};


} // namespace gproshan

#endif // TOPLESETS_H

