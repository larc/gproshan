#ifndef INPAINTING_H
#define INPAINTING_H

#include "mdict/dictionary.h"
#include "../mesh/che_poisson.h"
#include "../mesh/che_fill_hole.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class inpainting : public dictionary
{
	public:
		inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const real_t & _f, const bool & _plot = true);
		virtual ~inpainting() = default;

		void execute();
};


} // namespace gproshan::mdict

#endif // INPAINTING_H
