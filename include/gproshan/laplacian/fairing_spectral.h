#ifndef FAIRING_SPECTRAL_H
#define FAIRING_SPECTRAL_H

#include <gproshan/laplacian/fairing.h>


// geometry processing and shape analysis framework
namespace gproshan {


class fairing_spectral : public fairing
{
	public:
		size_t n_eigs;

	public:
		fairing_spectral(const size_t n_eigs_ = 100);
		virtual ~fairing_spectral() = default;

	private:
		void compute(che * mesh);
};


} // namespace gproshan

#endif // FAIRING_SPECTRAL_H

