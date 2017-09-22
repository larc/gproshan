#ifndef FAIRING_TAUBIN_H
#define FAIRING_TAUBIN_H

#include "fairing.h"

class fairing_taubin : public fairing
{
	private:
		matrix_t step;

	public:
		fairing_taubin(matrix_t step_ = 0.002);
		virtual ~fairing_taubin();

	private:
		void compute(che * shape);
};

#endif // FAIRING_TAUBIN_H

