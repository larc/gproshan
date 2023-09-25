#ifndef SPLAT_UTILS_H
#define SPLAT_UTILS_H

#include <gproshan/mesh/che.h>
#include <gproshan/geometry/mat.h>
#include <gproshan/raytracing/utils.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {

template <class T>
__host_device__
unsigned int expand_bits(T f);

template <class T>
__host_device__
unsigned int morton_2d(T x, T y);


template <class T>
__host_device__
T gaussian(const T & x, const T & std)
{
	return exp(- x * x / std);
}


template <class T>
struct splat_t
{
	index_t begin = 0;
	index_t end = 0;
	T radius = 1;
	vec<T, 3> center;
	mat<T, 3> tbn;

	__host_device__
	unsigned int morton2d(vec<T, 3> p) const
	{
		p = (tbn * (p - center)) / radius;
		p = (p + 1) / 2;
		unsigned int xx = expand_bits(p.x());
		unsigned int yy = expand_bits(p.y());
		return (xx << 1) | yy;
	}
};


struct splats_data
{
	unsigned int * morton_codes = nullptr;
	index_t * primID_splat = nullptr;

	splat_t<real_t> * splats = nullptr;
	size_t n_splats = 0;

	splats_data() = default;

	splats_data(const size_t & np, const size_t & ns): n_splats(ns)
	{
		morton_codes = new unsigned int[np];
		splats = new splat_t<real_t>[n_splats];
	}

	splats_data(splats_data && sd)
	{
		morton_codes = sd.morton_codes;
		primID_splat = sd.primID_splat;
		splats = sd.splats;
		n_splats = sd.n_splats;

		sd.morton_codes = nullptr;
		sd.primID_splat = nullptr;
		sd.splats = nullptr;
		sd.n_splats = 0;
	}

	~splats_data()
	{
		delete morton_codes;
		delete primID_splat;
		delete splats;
	}
};


template <class T>
__host_device__
index_t binary_search(const T * data, index_t i, index_t j, const T & value)
{
	while(i < j)
	{
		const index_t & m = (i + j) >> 1;
		if(data[m] == value)
			return m;

		data[m] < value ? i = m + 1 : j = m - 1;
	}

	return i;
}


template <class T>
__host_device__
void splat_hit(t_eval_hit<T> & hit, const CHE & pc, const splats_data & sd, const index_t & aprimID, const vec<T, 3> & x, const vec<T, 3> & d, const int & k)
{
	hit.primID = aprimID;
	const index_t sid = sd.primID_splat[hit.primID];
	const splat_t<T> & s = sd.splats[sid];

	index_t begin = s.begin;
	index_t end = s.end;

	const index_t & h = binary_search(sd.morton_codes, begin, end - 1, s.morton2d(x));
	const real_t & sigma2 = 0.1;

	vec<T, 3> & color = hit.Kd = {0, 0, 0};
	vec<T, 3> & normal = hit.normal = {0, 0, 0};
	vec<T, 3> & position = hit.position = {0, 0, 0};

	begin = h - k;
	end = h + k;
	begin = begin < s.begin || begin > end ? s.begin : begin;
	end = end > s.end ? s.end : end;

	T w, sum_w = 1e-5;
	for(index_t v = s.begin; v < s.end; ++v)
	{
		vec<T, 3> p = pc.GT[v];
		p.z() = p.x() * p.x() + p.y() * p.y();
		const vec<T, 3> & q = dot(d, p - x) * d + x;

		w = 1 - length(p - q);
		if(w < 0) w = 0;
		sum_w += w;

		const che::rgb_t & c = pc.VC[v];
		vec<T, 3> vc = {T(c.r), T(c.g), T(c.b)};
		vc /= 255;
		
		normal += w * pc.VN[v];
		color += w * vc;
		position += w * q;
	}

	normal = normalize(normal);
	color /= sum_w;
	position /= sum_w;

	return;

/*
	color = length(x - s.center) / s.radius;

	static int a = 0;
	if(!sid && a < 100)
	{
		gproshan_log_var(length(x - s.center) / s.radius);
		gproshan_log_var(s.morton2d(x));
		vertex p = (s.tbn * (x - s.center)) / s.radius;
		gproshan_log_var(((p + 1)/2) * 1024);
		gproshan_log_var((p.x() + 1) / 2);
		unsigned int aa, bb;
		gproshan_log_var(aa = expand_bits((p.x() + 1) / 2));
		gproshan_log_var(bb = expand_bits((p.y() + 1) / 2));
		gproshan_log_var((aa << 1) | bb);
		gproshan_log_var(morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2));
		++a;
	}
*/
}


// FROM: https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
template <class T>
__host_device__
unsigned int expand_bits(T f)
{
	f *= 1024;
	f = f < 0 ? 0 : f;
	f = f > 1023 ? 1023 : f;

	unsigned int v = (unsigned int) f;
    v = (v | (v << 8)) & 0x00FF00FFu;
    v = (v | (v << 4)) & 0x0F0F0F0Fu;
    v = (v | (v << 2)) & 0x33333333u;
    v = (v | (v << 1)) & 0x55555555u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 2D point located within the unit square [0,1].
template <class T>
__host_device__
unsigned int morton_2d(T x, T y)
{
	unsigned int xx = expand_bits(x);
	unsigned int yy = expand_bits(y);
	return (xx << 1) | yy;
}


} // namespace gproshan

#endif // SPLAT_UTILS_H

