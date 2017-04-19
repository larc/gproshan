#ifndef FAIRING_SPECTRAL_H
#define FAIRING_SPECTRAL_H

#include "fairing.h"

class fairing_spectral : public fairing
{
	private:
		size_t k;

	public:
		fairing_spectral(size_t k_ = 10);
		virtual ~fairing_spectral();

	private:
		void compute(che * shape);
};

#endif // FAIRING_SPECTRAL_H
