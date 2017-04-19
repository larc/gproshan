#ifndef REPEATABILITY_H
#define REPEATABILITY_H

#include "off.h"

typedef double rep_t;

class repeatability
{
	private:
		rep_t rep;

	public:
		repeatability(off & shape_x, off & shape_y);

	private:
		void get_correspondence(off & shape_x, off & shape_y);

};

#endif // REPEATABILITY_H
