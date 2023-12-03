#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>
#include <gproshan/scenes/scene.h>
#include <gproshan/raytracing/render_params.h>


#ifdef GPROSHAN_OPTIX

#include <optix.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
	vec4 * color_buffer = nullptr;
	unsigned int buffer_size	= 0;
	uvec2 window_size;
	uvec2 viewport_pos;
	unsigned int depth			= 1;
	unsigned int n_frames		= 0;
	unsigned int n_lights		= 0;
	light lights[NL];
	light ambient;
	vec3 cam_pos;
	mat4 inv_proj_view;
	bool flat;

	OptixTraversableHandle traversable;

	scene_data sc;

	void * other = nullptr;
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_PARAMS_H

