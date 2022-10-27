#ifndef SPLAT_UTILS_H
#define SPLAT_UTILS_H

#include <gproshan/mesh/che.h>
#include <gproshan/geometry/mat.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct splats_data
{
	CHE * pc = nullptr;
	unsigned int * morton_codes = nullptr;
	unsigned int * idx_splats = nullptr;
	unsigned int n_splats = 0;
	mat4 model_mat;
};


// FROM: https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
template <class T>
__host__ __device__
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
__host__ __device__
unsigned int morton_2d(T x, T y)
{
	unsigned int xx = expand_bits(x * 1023 + 0.5);
	unsigned int yy = expand_bits(y * 1023 + 0.5);
	return (xx >> 1) + yy;
}


} // namespace gproshan

#endif // SPLAT_UTILS_H

