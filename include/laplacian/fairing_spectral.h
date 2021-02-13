#ifndef FAIRING_SPECTRAL_H
#define FAIRING_SPECTRAL_H

#include "laplacian/fairing.h"


// geometry processing and shape analysis framework
namespace gproshan {


class fairing_spectral : public fairing
{
	private:
		size_t k;

	public:
		fairing_spectral(const size_t & k_ = 10);
		virtual ~fairing_spectral();

	private:
		void compute(che * mesh);
};


} // namespace gproshan

#endif // FAIRING_SPECTRAL_H
