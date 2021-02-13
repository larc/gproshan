#include "features/descriptor.h"

#include "laplacian/laplacian.h"

// geometry processing and shape analysis framework
namespace gproshan {


descriptor::descriptor(const signature & sig, const che * mesh, const size_t & n_eigs)
{
	if(!eigs_laplacian(mesh, eigval, eigvec, L, A, n_eigs))
		return;
	
	switch(sig)
	{
		case GPS: compute_gps(1); break;
		case HKS: break;
		case WKS: break;
	}
}

descriptor::operator bool () const
{
	return features.size() > 0;
}

real_t descriptor::operator () (const index_t & v) const
{
	return norm(features.row(v));
}

void descriptor::compute_gps(const size_t & T)
{
	features = eigvec.tail_cols(eigvec.n_cols - 1);
	for(index_t i = 1; i < eigval.size(); i++)
		features.col(i - 1) /= sqrt(eigval(i));
}

void descriptor::compute_hks(const size_t & T)
{
	features.zeros();
}

void descriptor::compute_wks(const size_t & T)
{
	features.zeros();
}


} // namespace gproshan

