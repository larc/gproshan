#ifndef SUPER_RESOLUTION_H
#define SUPER_RESOLUTION_H

#include "dictionary.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class super_resolution : public dictionary
{
	public:
		super_resolution(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool &_learn, const bool & _plot = true);
		virtual ~super_resolution() = default;

		distance_t execute();
};


} // namespace gproshan::mdict

#endif // SUPER_RESOLUTION_H

