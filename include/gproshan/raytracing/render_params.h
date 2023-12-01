#ifndef RT_RENDER_PARAMS_H
#define RT_RENDER_PARAMS_H

#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>
#include <gproshan/raytracing/light.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


const size_t NL = 16;	// number of lights

struct render_params
{
	int depth = 1;
	int window_width = 0;
	int window_height = 0;
	int viewport_width = 0;
	int viewport_height = 0;
	int viewport_x = 0;
	int viewport_y = 0;
	int n_lights = 0;
	light lights[NL];
	light ambient = {0, 1, 0.1};
	vertex cam_pos;
	mat4 inv_proj_view;
	bool restart = false;
	bool viewport_is_window = true;

	bool add_light(const light & l)
	{
		if(n_lights == NL)
			return false;

		lights[n_lights++] = l;
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


#endif // RT_RENDER_PARAMS_H

