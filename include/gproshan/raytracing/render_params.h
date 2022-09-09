#ifndef RENDER_PARAMS_H
#define RENDER_PARAMS_H


#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct render_params
{
	int window_width = 0;
	int window_height = 0;
	int viewport_width = 0;
	int viewport_height = 0;
	int viewport_x = 0;
	int viewport_y = 0;
	bool restart = false;
	bool viewport_is_window = true;
	vertex cam_pos;
	vertex lights[16];
	unsigned int n_lights = 0;
	mat4 inv_proj_view;

	bool add_light(const vertex & light)
	{
		if(n_lights == sizeof(lights) / sizeof(vertex))
			return false;

		lights[n_lights++] = light;
		return true;
	}
};


} // namespace gproshan


#endif // RENDER_PARAMS_H

