#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include <optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
	struct
	{
		void * color_buffer = nullptr;
		uint32_t width, height;
	} frame;

	void * camera = nullptr;

	OptixTraversableHandle traversable;
};


} // namespace gproshan

#endif // RT_OPTIX_PARAMS_H

#endif // GPROSHAN_OPTIX

