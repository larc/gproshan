#ifndef FAIRING_TAUBIN_H
#define FAIRING_TAUBIN_H

#include <gproshan/laplacian/fairing.h>


// geometry processing and shape analysis framework
namespace gproshan {


class fairing_taubin : public fairing
{
	public:
		float step;

	public:
		fairing_taubin(const float step_ = 0.001);
		virtual ~fairing_taubin() = default;

	private:
		void compute(che * mesh);
};


} // namespace gproshan

#endif // FAIRING_TAUBIN_H

