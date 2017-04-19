#ifndef FAIRING_H
#define FAIRING_H

#include "che.h"

class fairing
{
	protected:
		vertex * positions;

	public:
		fairing();
		virtual ~fairing();
		void run(che * shape);
		vertex * get_postions();

	protected:
		virtual void compute(che * shape) = 0;
};

#endif // FAIRING_H

