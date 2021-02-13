#include "features/descriptor.h"

#include "laplacian/laplacian.h"

// geometry processing and shape analysis framework
namespace gproshan {


descriptor::descriptor(const signature & sig, const che * mesh, const size_t & n_eigs)
{
	if(eigs_laplacian(mesh, eigval, eigvec, L, A, n_eigs) == 0)
		return;
}


} // namespace gproshan

