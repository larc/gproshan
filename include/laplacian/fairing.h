#ifndef FAIRING_H
#define FAIRING_H

#include "mesh/che.h"


// geometry processing and shape analysis framework
namespace gproshan {


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


} // namespace gproshan

#endif // FAIRING_H

