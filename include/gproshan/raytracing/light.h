#ifndef LIGHT_H
#define LIGHT_H

#include <gproshan/geometry/vec.h>


// geometry processing and shape analysis framework
namespace gproshan {


struct light
{
	vec3 pos = 0;
	vec3 color = 1;
	float power = 10;
};


} // namespace gproshan

#endif // LIGHT_H

