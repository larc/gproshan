#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include <gproshan/include.h>


#ifdef GPROSHAN_OPTIX

#include <optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
	void * color_buffer = nullptr;
	uint32_t window_width, window_height;
	uint32_t viewport_width, viewport_height;
	uint32_t viewport_x, viewport_y;
	uint32_t n_lights = 0;
	float lights[3][16];
	float cam_pos[3];
	float inv_proj_view[16];
	bool flat;

	OptixTraversableHandle traversable;
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_PARAMS_H

