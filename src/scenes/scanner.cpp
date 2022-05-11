#include "scenes/scanner.h"


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


che * scanner_ptx(const raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam)
{
	std::vector<vertex> ptx;

	return new che(ptx.data(), ptx.size(), nullptr, 0);
}

che * scanner_ptx(const che * mesh, const raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam)
{
	std::vector<vertex> ptx;
	//index_t v = rt->cast_ray(cam_pos, glm::normalize(p - cam_pos));
	//cast_ray_intersect_depth

	// use the raytracing to find the hit points
	// save them in ptx.data

	glm::vec3 cam_pos = glm::vec3(cam.x, cam.y, cam.z);
	glm::vec3 p;

	std::vector<vertex> vertices;
	const real_t  r = 1;
	index_t v;
	float distance;

	const real_t delta_phi = M_PI / n_rows;
	const real_t delta_theta = (2 * M_PI) / n_cols;

	for(real_t phi = 0; phi < 2 * M_PI - 0.5 * delta_phi; phi += delta_phi)
	for(real_t theta = delta_theta; theta < M_PI - 0.5 * delta_theta; theta += delta_theta)
	{
		p = glm::vec3( r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta) );
		auto [v, distance] = rt->cast_ray_intersect_depth(cam_pos, glm::normalize(p - cam_pos));
	}	

	return new che(ptx.data(), ptx.size(), nullptr, 0);
}


} // namespace gproshan

