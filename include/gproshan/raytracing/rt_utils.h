#ifndef RT_UTILS_H
#define RT_UTILS_H

#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


template <class T, uint32_t N = 16>
struct random
{
	uint32_t previous;
	
	__host__ __device__
	random(uint32_t v0, uint32_t v1)
	{
		uint32_t s = 0;
		for(uint32_t i = 0; i < N; ++i)
		{
			s += 0x9e3779b9;
			v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s) ^ ((v1 >> 5) + 0xc8013ea4);
			v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s) ^ ((v0 >> 5) + 0x7e95761e);
		}
		previous = v0;
	}

	__host__ __device__
	T operator () ()
	{
		previous = previous * 1664525 + 1013904223;
		return T(previous & 0x00FFFFFF) / T(0x01000000);
	}
};


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

