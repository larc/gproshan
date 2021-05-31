#ifndef FAIRING_TAUBIN_H
#define FAIRING_TAUBIN_H

#include "laplacian/fairing.h"


// geometry processing and shape analysis framework
namespace gproshan {


class fairing_taubin : public fairing
{
	private:
		real_t step;

	public:
		fairing_taubin(const real_t & step_ = 0.01);
		virtual ~fairing_taubin() = default;

	private:
		void compute(che * mesh);
};


} // namespace gproshan

#endif // FAIRING_TAUBIN_H

