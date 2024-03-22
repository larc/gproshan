#include <gproshan/features/descriptor.h>

#include <gproshan/laplacian/laplacian.h>


// geometry processing and shape analysis framework
namespace gproshan {


descriptor::descriptor(const signature & sig, const che * mesh, const size_t n_eigs)
{
	if(!eigs_laplacian(mesh, eigval, eigvec, L, A, n_eigs))
		return;

	switch(sig)
	{
		case GPS: compute_gps(); break;
		case HKS: compute_hks(); break;
		case WKS: compute_wks(); break;
	}
}

size_t descriptor::n_eigs()
{
	return eigval.n_elem;
}

descriptor::operator bool () const
{
	return features.size() > 0;
}

real_t descriptor::operator () (const index_t v) const
{
	return norm(features.row(v));
}

void descriptor::compute_gps()
{
	features = eigvec.tail_cols(eigvec.n_cols - 1);
	for(index_t i = 1; i < eigval.size(); ++i)
		features.col(i - 1) /= sqrt(eigval(i));
}

void descriptor::compute_hks()
{
	eigvec = eigvec % eigvec; 						// element wise product

	features.zeros(eigvec.n_rows, eigvec.n_cols);

	#pragma omp parallel for
	for(index_t t = 0; t < features.n_cols; ++t)
		features.col(t) = eigvec * exp(-eigval * t);
}

///< http://imagine.enpc.fr/~aubrym/projects/wks/index.html
void descriptor::compute_wks()
{
	eigvec = eigvec % eigvec; 						// element wise product
	eigval = log(eigval);

	arma::vec e = arma::linspace<arma::vec>(eigval(1), eigval(eigval.n_elem - 1), eigval.n_elem);
	real_t sigma = (e(1) - e(0)) * 6;				// 6 is wks variance see reference
	real_t sigma_2 = 2 * sigma * sigma;

	features.zeros(eigvec.n_rows, e.n_elem);

	#pragma omp parallel for
	for(index_t t = 0; t < features.n_cols; ++t)
		features.col(t) = eigvec * exp(-pow(e(t) - eigval, 2) / sigma_2) /
							sum(exp(-pow(e(t) - eigval, 2) / sigma_2));
}


} // namespace gproshan

