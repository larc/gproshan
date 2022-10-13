#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>


#ifdef GPROSHAN_OPTIX

#include <optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
	vec4 * color_buffer = nullptr;
	int buffer_size = 0;
	int window_width = 0;
	int window_height = 0;
	int viewport_x = 0;
	int viewport_y = 0;
	int n_samples = 0;
	int n_lights = 0;
	vec3 lights[16];
	vec3 cam_pos;
	mat4 inv_proj_view;
	bool flat;

	OptixTraversableHandle traversable;
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_PARAMS_H

