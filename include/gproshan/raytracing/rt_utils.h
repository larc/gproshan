#ifndef RT_UTILS_H
#define RT_UTILS_H

#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


template <class T>
__host__ __device__
vec<T, 3> ray_view_dir(const ivec2 & pos, const ivec2 & windows_size, const mat<T, 4> & inv_proj_view, const vec<T, 3> & cam_pos)
{
	vec2 screen = {	(float(pos.x()) + 0.5f) / windows_size.x(),
					(float(pos.y()) + 0.5f) / windows_size.y()
					};
	vec<T, 4> view = {screen.x() * 2 - 1, screen.y() * 2 - 1, 1, 1};
	vec<T, 4> q = inv_proj_view * view;
	vec<T, 3> p = q / q[3];

	return normalize(p - cam_pos);
}


} // namespace gproshan

#endif // RT_UTILS_H

