#ifndef RENDER_PARAMS_H
#define RENDER_PARAMS_H


#include <gproshan/include.h>

#include <glm/glm.hpp>


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
	glm::mat4 proj_view_mat;
	glm::vec3 cam_pos;
	std::vector<glm::vec3> lights;
	bool viewport_is_window = true;
};


} // namespace gproshan


#endif // RENDER_PARAMS_H

