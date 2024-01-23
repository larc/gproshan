#ifndef RT_UTILS_H
#define RT_UTILS_H

#include <gproshan/include.h>
#include <gproshan/geometry/mat.h>

#include <gproshan/mesh/che.cuh>
#include <gproshan/scenes/scene.h>
#include <gproshan/raytracing/light.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


template <class T, unsigned int N = 16>
struct random
{
	unsigned int previous;

	__host_device__
	random(unsigned int p): previous(p) {}

	__host_device__
	random(unsigned int v0, unsigned int v1)
	{
		unsigned int s = 0;
		for(unsigned int i = 0; i < N; ++i)
		{
			s += 0x9e3779b9;
			v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s) ^ ((v1 >> 5) + 0xc8013ea4);
			v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s) ^ ((v0 >> 5) + 0x7e95761e);
		}
		previous = v0;
	}

	__host_device__
	T operator () ()
	{
		previous = previous * 1664525 + 1013904223;
		return T(previous & 0x00FFFFFF) / T(0x01000000);
	}

	__host_device__
	operator unsigned int & ()
	{
		return previous;
	}
};

template<class T>
__host_device__
vec<T, 3> texture(const scene::texture & tex, const vec<T, 2> & coord)
{
	const int i = (tex.width + int(coord.x() * (tex.width - 1))) % tex.width;
	const int j = (tex.height + int(coord.y() * (tex.height -1))) % tex.height;
	const int k = j * tex.width + i;

	che::rgb_t color;
	if(tex.spectrum == 3)
	{
		const che::rgb_t * img = (const che::rgb_t *) tex.data;
		color = img[k];
	}
	if(tex.spectrum == 1)
	{
		color.r = color.g = color.b = tex.data[k];
	}

	return {T(color.r) / 255, T(color.g) / 255, T(color.b) / 255};
}

template <class T>
struct t_eval_hit
{
	index_t primID = NIL;
	int illum = 1;
	T u = 0, v = 0;
	T dist = 0;
	T heatmap = 0;
	vec<T, 3> position;
	vec<T, 3> normal;
	vec<T, 3> Ka = 1;
	vec<T, 3> Kd = 0.5;
	vec<T, 3> Ks = 0.2;
	T Ns = 10;
	T Ni = 0;
	T d = 1;

	__host_device__
	t_eval_hit() {}

	__host_device__
	t_eval_hit(const CHE & mesh, const index_t aprimID, const T & au, const T & av, const scene_data & sc)
	{
		primID = aprimID;
		u = au;
		v = av;

		if(!mesh.n_trigs) // pointcloud
		{
			Kd = {T(mesh.VC[aprimID].r), T(mesh.VC[aprimID].g), T(mesh.VC[aprimID].b)};
			Kd /= 255;
			normal = mesh.VN[aprimID];
			heatmap = mesh.VHC[aprimID];
			return;
		}

		const index_t he = primID * che::mtrig;

		const index_t a = mesh.VT[he];
		const index_t b = mesh.VT[he + 1];
		const index_t c = mesh.VT[he + 2];

		const vec<T, 3> ca = {T(mesh.VC[a].r), T(mesh.VC[a].g), T(mesh.VC[a].b)};
		const vec<T, 3> cb = {T(mesh.VC[b].r), T(mesh.VC[b].g), T(mesh.VC[b].b)};
		const vec<T, 3> cc = {T(mesh.VC[c].r), T(mesh.VC[c].g), T(mesh.VC[c].b)};

		Kd = ((1.f - u - v) * ca + u * cb + v * cc) / 255;
		normal = (1.f - u - v) * mesh.VN[a] + u * mesh.VN[b] + v * mesh.VN[c];
		heatmap = (1.f - u - v) * mesh.VHC[a] + u * mesh.VHC[b] + v * mesh.VHC[c];

		if(!sc.trig_mat) return;
		if(sc.trig_mat[primID] == NIL) return;

		const scene::material & mat = sc.materials[sc.trig_mat[primID]];
		vec<T, 2> texcoord;
		if(sc.texcoords)
			texcoord = (1.f - u - v) * sc.texcoords[a] + u * sc.texcoords[b] + v * sc.texcoords[c];

		Ka = mat.Ka;
		if(mat.map_Ka != -1)
			Ka = texture(sc.textures[mat.map_Ka], texcoord);

		Kd = mat.Kd;
		if(mat.map_Kd != -1)
			Kd = texture(sc.textures[mat.map_Kd], texcoord);

		Ks = mat.Ks;
		if(mat.map_Ks != -1)
			Ks = texture(sc.textures[mat.map_Ks], texcoord);

		Ns = mat.Ns;
		Ni = mat.Ni;

		d = mat.d;
//		if(mat.map_d != -1)
//			d = texture(sc.textures[mat.map_d], texcoord);

		illum = mat.illum;
	}

	//	PTX symbols of certain types (e.g. pointers to functions) cannot be used to initialize array
	__host_device__
	bool scatter_mat(vec<T, 3> & dir, random<T> & rnd)		// dir in: v (view), out: scattered
	{
		switch(illum)
		{
			case 3:
			case 5:
			case 6:
			case 7:
				return scatter_reflect(dir, rnd);
		}

		return scatter_diffuse(dir, rnd);
	}

	__host_device__
	bool scatter_reflect(vec<T, 3> & dir, random<T> & )
	{
		dir = normalize(dir - 2 * dot(dir, normal) * normal);
		return dot(dir, normal) > 0;
	}

	__host_device__
	bool scatter_refract(vec<T, 3> & dir, random<T> & )
	{
		const float dvn = dot(dir, normal);
		const float d = 1 - Ni * Ni * (1 - dvn * dvn);

		if(d <= 0) return false;

		dir = Ni * (dir - dvn * normal) - normal * sqrtf(d);
		return true;
	}

	__host_device__
	bool scatter_diffuse(vec<T, 3> & dir, random<T> & rnd)
	{
		// random unit sphere
		const T & theta = rnd() * 2 * M_PI;
		const T & phi = acosf(2 * rnd() - 1);
		const T & r = cbrtf(rnd());

		const vec<T, 3> p = { r * sinf(phi) * cosf(theta)
							, r * sinf(phi) * sinf(theta)
							, r * cosf(phi)
							};

		dir = normalize(normal + p);

		return true;
	}
};

template <class T, class Occluded>
__host_device__
vec<T, 3> eval_li(const t_eval_hit<T> & hit, const light & ambient, const light * lights, const int n_lights, const vec<T, 3> & eye, Occluded occluded)
{
	const vec<T, 3> & v = normalize(eye - hit.position);
	const vec<T, 3> & n = hit.normal;

	T lambertian;
	T specular;
	vec<T, 3> li, l, h;

	for(int i = 0; i < n_lights; ++i)
	{
		const light & L = lights[i];

		l = L.pos - hit.position;
		const T & r = length(l);

		l /= r;
		h = normalize(l + v);

	#ifdef __CUDACC__
		lambertian = max(dot(l, n), 0.f);
		specular = powf(max(dot(h, n), 0.f), hit.Ns);
	#else
		lambertian = std::max(dot(l, n), 0.f);
		specular = powf(std::max(dot(h, n), 0.f), hit.Ns);
	#endif // __CUDACC__

		const vec<T, 3> & color = hit.Ka * ambient.color * ambient.power +
									(lambertian * hit.Kd + specular * hit.Ks) * L.color * L.power / (r * r);

		li += (dot(v, n) < 0 || occluded(hit.position, l, r) ? 0.4f : 1.0f) * color;
	}

	return li / n_lights;
}


using eval_hit = t_eval_hit<float>;


template <class T>
__host_device__
vec<T, 3> ray_view_dir(	const uvec2 & pos,
						const uvec2 & windows_size,
						const mat<T, 4> & inv_proj_view,
						const vec<T, 3> & cam_pos,
						random<T> & rnd
						)
{
	vec2 screen = {	(float(pos.x()) + rnd()) / windows_size.x(),
					(float(pos.y()) + rnd()) / windows_size.y()
					};
	vec<T, 4> view = {screen.x() * 2 - 1, screen.y() * 2 - 1, 1, 1};
	vec<T, 4> q = inv_proj_view * view;
	vec<T, 3> p = q / q[3];

	return normalize(p - cam_pos);
}


template <class H>
__host_device__
index_t closest_hit_vertex(const CHE & mesh, const H & hit)
{
	if(!mesh.n_trigs) return hit.primID;

	index_t he = 0;
	real_t w = 1 - hit.u - hit.v;

	if(w < hit.u)
	{
		he = 1;
		w = hit.u;
	}

	if(w < hit.v)
	{
		he = 2;
		w = hit.v;
	}

	return mesh.VT[hit.primID * 3 + he];
}


} // namespace gproshan

#endif // RT_UTILS_H

