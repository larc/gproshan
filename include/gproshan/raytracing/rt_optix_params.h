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
	
	float light[3];
	float cam_pos[3];
	float inv_proj_view[16];

	OptixTraversableHandle traversable;
};


} // namespace gproshan

#endif // RT_OPTIX_PARAMS_H

#endif // GPROSHAN_OPTIX

