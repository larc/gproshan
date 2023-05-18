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


struct splats_data
{
	CHE * pc = nullptr;
	unsigned int * morton_codes = nullptr;
	unsigned int * primID_splat = nullptr;

	unsigned int n_splats = 0;
	unsigned int * idx_splats = nullptr;
	vertex * centers = nullptr;
	mat3 * tbns = nullptr;
};

template <class T>
__host_device__
index_t binary_search(const T * data, index_t i, index_t j, const T & value)
{
	index_t m = 0;
	while(i < j)
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
void splat_hit(t_eval_hit<T> & hit, const splats_data * splat, const index_t & aprimID, const vec<T, 3> & x)
{
	hit.primID = aprimID;
	const index_t s = splat->primID_splat[hit.primID];
		
	const index_t begin = splat->idx_splats[s];
	const index_t end = splat->idx_splats[s + 1];

	const vec<T, 3> & center = splat->centers[s];
	const mat<T, 3> & tbn = splat->tbns[s];

	const vec<T, 3> & p = tbn * (x - center);  
	unsigned int code = morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2);

	const index_t h = binary_search(splat->morton_codes, begin, end - 1, code);
	
	hit.normal = splat->pc->VN[h];

	const che::rgb_t & color = splat->pc->VC[h];
	hit.Kd = {T(color.r), T(color.g), T(color.b)};
	hit.Kd /= 255;

	hit.position = x;
}


// FROM: https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
template <class T>
__host_device__
unsigned int expand_bits(const T & fv)
{
	unsigned int v = (unsigned int) fv;
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 2D point located within the unit square [0,1].
template <class T>
__host_device__
unsigned int morton_2d(T x, T y)
{
	unsigned int xx = expand_bits(x * 1023 + 0.5);
	unsigned int yy = expand_bits(y * 1023 + 0.5);
	return (xx >> 1) + yy;
}


} // namespace gproshan

#endif // SPLAT_UTILS_H

