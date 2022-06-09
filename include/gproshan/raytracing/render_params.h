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
	mat4 proj_view_mat;
	vertex cam_pos;
	std::vector<vertex> lights;
};


} // namespace gproshan


#endif // RENDER_PARAMS_H

