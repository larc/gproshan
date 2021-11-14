#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include "mesh/vertex.h"

#include <optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
	struct
	{
		uint32_t * colorBuffer = nullptr;
		uint32_t width, height;
	} frame;

	struct
	{
		vertex_cu position;
		vertex_cu direction;
		vertex_cu horizontal;
		vertex_cu vertical;
	} camera;

	OptixTraversableHandle traversable;
};


} // namespace gproshan

#endif // RT_OPTIX_PARAMS_H

#endif // GPROSHAN_OPTIX

