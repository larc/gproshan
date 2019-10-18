#ifndef INPAINTING_H
#define INPAINTING_H

#include "dictionary.h"
#include "../che_poisson.h"
#include "../che_fill_hole.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class inpainting : public dictionary
{
	public:
		inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, const bool & _plot = true);
		virtual ~inpainting() = default;

		distance_t execute();
};


} // namespace gproshan::mdict

#endif // INPAINTING_H
