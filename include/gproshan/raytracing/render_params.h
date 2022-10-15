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
	int n_lights = 0;
	vertex lights[16];
	vertex cam_pos;
	mat4 inv_proj_view;
	bool restart = false;
	bool viewport_is_window = true;

	bool add_light(const vertex & light)
	{
		if(n_lights == sizeof(lights) / sizeof(vertex))
			return false;

		lights[n_lights++] = light;
		return true;
	}

	void log()
	{
		gproshan_log_var(window_width);
		gproshan_log_var(window_height);
		gproshan_log_var(viewport_width);
		gproshan_log_var(viewport_height);
		gproshan_log_var(viewport_x);
		gproshan_log_var(viewport_y);
		gproshan_log_var(n_lights);
		gproshan_log_var(cam_pos);
		gproshan_log_var(restart);
		gproshan_log_var(viewport_is_window);
	}
};


} // namespace gproshan


#endif // RENDER_PARAMS_H

