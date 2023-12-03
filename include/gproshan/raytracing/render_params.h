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
	uvec2 window_size;
	uvec2 viewport_size;
	uvec2 viewport_pos;
	unsigned int depth = 1;
	unsigned int n_frames = 0;
	unsigned int n_lights = 0;
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
		gproshan_log_var(window_size);
		gproshan_log_var(viewport_size);
		gproshan_log_var(viewport_pos);
		gproshan_log_var(n_lights);
		gproshan_log_var(cam_pos);
		gproshan_log_var(restart);
		gproshan_log_var(viewport_is_window);
	}
};


} // namespace gproshan


#endif // RT_RENDER_PARAMS_H

