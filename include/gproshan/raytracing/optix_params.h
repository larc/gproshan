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


struct launch_params: public base_params
{
	OptixTraversableHandle traversable;

	scene_data sc;

	bool flat;
	void * other = nullptr;
	vec4 * color_buffer = nullptr;
	unsigned int buffer_size = 0;
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_OPTIX_PARAMS_H

