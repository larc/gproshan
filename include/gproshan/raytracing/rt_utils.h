#ifndef RT_UTILS_H
#define RT_UTILS_H

#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>

#include <gproshan/mesh/che.cuh>


// geometry processing and shape analysis framework
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
struct t_eval_hit
{
	unsigned int primID = NIL;
	T dist = 0;
	T u = 0, v = 0;
	vec<T, 3> position;
	vec<T, 3> normal;
	vec<T, 3> color;

	__host__ __device__
	t_eval_hit() {}

	__host__ __device__
	t_eval_hit(const CHE & mesh, const unsigned int & aprimID, const T & au, const T & av)
	{
		primID = aprimID;
		u = au;
		v = av;

		const int he = primID * che::mtrig;

		const int a = mesh.VT[he];
		const int b = mesh.VT[he + 1];
		const int c = mesh.VT[he + 2];

		const vertex ca = {float(mesh.VC[a].r), float(mesh.VC[a].g), float(mesh.VC[a].b)};
		const vertex cb = {float(mesh.VC[b].r), float(mesh.VC[b].g), float(mesh.VC[b].b)};
		const vertex cc = {float(mesh.VC[c].r), float(mesh.VC[c].g), float(mesh.VC[c].b)};

		color = ((1.f - u - v) * ca + u * cb + v * cc) / 255;
		normal = (1.f - u - v) * mesh.VN[a] + u * mesh.VN[b] + v * mesh.VN[c];
	}
};

template <class T, class Occluded>
__host__ __device__
vec<T, 3> eval_li(const t_eval_hit<T> & hit, const vec<T, 3> * lights, const int & n_lights, Occluded occluded)
{
	vec<T, 3> li, wi;
	float light_dist, dot_wi_normal;

	for(int i = 0; i < n_lights; ++i)
	{
		wi = lights[i] - hit.position;
		light_dist = length(wi);

		wi /= light_dist;
		dot_wi_normal = dot(wi, hit.normal);

		li += (dot_wi_normal < 0 ? -dot_wi_normal : dot_wi_normal) * (occluded(hit.position, wi, light_dist) ? 0.4f : 1.0f) * hit.color;
	}

	return li / n_lights;
}


using eval_hit = t_eval_hit<float>;


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

