#ifndef SUPER_RESOLUTION_H
#define SUPER_RESOLUTION_H

#include "dictionary.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class super_resolution : public dictionary
{
	public:
		super_resolution(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot = true);
		virtual ~super_resolution() = default;

		void execute();
};

} // mdict

#endif // SUPER_RESOLUTION_H

