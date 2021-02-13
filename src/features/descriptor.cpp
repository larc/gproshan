#include "features/descriptor.h"

#include "laplacian/laplacian.h"

// geometry processing and shape analysis framework
namespace gproshan {


descriptor::descriptor(const signature & sig, const che * mesh, const size_t & n_eigs)
{
	if(!compute_eigs(mesh, n_eigs))
		return;
}

bool descriptor::compute_eigs(const che * mesh, const size_t & k)
{
	if(eigs_laplacian(mesh, eigval, eigvec, L, A, k) == 0)
		return false;
	

}


} // namespace gproshan

