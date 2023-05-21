#ifndef SPLAT_UTILS_H
#define SPLAT_UTILS_H

#include <gproshan/mesh/che.h>
#include <gproshan/geometry/mat.h>
#include <gproshan/raytracing/utils.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


template <class T>
__host_device__
unsigned int morton_2d(T x, T y);


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
		return morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2);
	}
};


struct splats_data
{
	CHE * pc = nullptr;
	unsigned int * morton_codes = nullptr;
	index_t * primID_splat = nullptr;

	splat_t<real_t> * splats = nullptr;
	size_t n_splats = 0;

	splats_data(const size_t & np, const size_t & ns): n_splats(ns)
	{
		morton_codes = new unsigned int[np];
		splats = new splat_t<real_t>[n_splats];
	}

	~splats_data()
	{
		delete pc;
		delete morton_codes;
		delete primID_splat;
		delete splats;
	}
};


template <class T>
__host_device__
int binary_search(const T * data, int i, int j, const T & value)
{
	int m = 0;
	while(i <= j)
	{
		m = (i + j) >> 1;
		if(data[m] == value)
			break;

		data[m] < value ? i = m + 1 : j = m - 1;
	}

	return m;
}


template <class T>
__host_device__
void splat_hit(t_eval_hit<T> & hit, const splats_data * sd, const index_t & aprimID, const vec<T, 3> & x, const int & k)
{
	hit.primID = aprimID;
	const index_t sid = sd->primID_splat[hit.primID];
	const splat_t<T> & s = sd->splats[sid];

	index_t begin = s.begin;
	index_t end = s.end;

	const int h = binary_search(sd->morton_codes, begin, end - 1, s.morton2d(x));
	T sigma = length(x - sd->pc->GT[h]);
	sigma *= sigma;

	vec<T, 3> & color = hit.Kd = {0, 0, 0};
	vec<T, 3> & normal = hit.normal = {0, 0, 0};

	begin = h - k >= begin ? h - k : begin;
	end = h + k <= end ? h + k : end;

	T w, sum_w = 0;
	for(index_t v = h; v < h + 1; ++v)
	{
		w = length(x - sd->pc->GT[v]);
		w = exp(-0.5 * w * w / sigma);
		sum_w += w;

		const che::rgb_t & c = sd->pc->VC[v];
		vec<T, 3> vc = {T(c.r), T(c.g), T(c.b)};
		vc /= 255;
		normal += w * sd->pc->VN[v];
		color += w * vc;
	}

	normal /= sum_w;
	color /= sum_w;

	hit.position = x;
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

