#ifndef SYNTHESIS_H
#define SYNTHESIS_H

#include "dictionary.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class synthesis : public dictionary
{
	public:
		synthesis(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, const bool & _plot = true);
		virtual ~synthesis() = default;

		distance_t execute();
};


} // namespace gproshan::mdict

#endif // SYNTHESIS_H

