#ifndef DENOISING_H
#define DENOISING_H

#include "dictionary.h"

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class denoising : public dictionary
{
	public:
		denoising(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot = true);
		virtual ~denoising() = default;

		void execute();
};

} // mdict

#endif // DENOISING_H

