#ifdef GPROSHAN_OPTIX


#include "mesh/che.h"
#include "mesh/vertex.cuh"
#include "raytracing/rt_optix_params.h"


#include <optix_device.h>
#include <cuda_runtime.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


__host__ __device__
vertex_cu operator * (const real_t & a, const vertex_cu & v)
{
	return vertex_cu(a * v.x, a * v.y, a * v.z);
}

__host__ __device__
vertex_cu operator + (const real_t & a, const vertex_cu & v)
{
	return vertex_cu(a + v.x, a + v.y, a + v.z);
}

__host__ __device__
vertex_cu normalize (const vertex_cu & v)
{
	return v / *v;
}


extern "C" __constant__ launch_params optixLaunchParams;

static __forceinline__ __device__
void * unpackPointer(uint32_t i0, uint32_t i1)
{
	const uint64_t uptr = static_cast<uint64_t>(i0) << 32 | i1;
	void * ptr = reinterpret_cast<void*>(uptr);
	return ptr;
}

static __forceinline__ __device__
void packPointer(void * ptr, uint32_t & i0, uint32_t & i1)
{
	const uint64_t uptr = reinterpret_cast<uint64_t>(ptr);
	i0 = uptr >> 32;
	i1 = uptr & 0x00000000ffffffff;
}

template<typename T>
static __forceinline__ __device__ T * getPRD()
{
	const uint32_t u0 = optixGetPayload_0();
	const uint32_t u1 = optixGetPayload_1();
	return reinterpret_cast<T*>(unpackPointer(u0, u1));
}


extern "C" __global__ void __closesthit__shadow() {}

extern "C" __global__ void __closesthit__radiance()
{
	const CHE & mesh = **(const CHE **) optixGetSbtDataPointer();

	const int primID = optixGetPrimitiveIndex();
	const int he = primID * che::mtrig;
	const float u = optixGetTriangleBarycentrics().x;
	const float v = optixGetTriangleBarycentrics().y;

	const int a = mesh.VT[he];
	const int b = mesh.VT[he + 1];
	const int c = mesh.VT[he + 2];

	const vertex_cu & A = mesh.GT[a];
	const vertex_cu & B = mesh.GT[b];
	const vertex_cu & C = mesh.GT[c];

	const vertex_cu Ng = normalize((B - A) * (C - A));
	const vertex_cu normal = mesh.VN ? (1.f - u - v) * mesh.VN[a] + u * mesh.VN[b] + v * mesh.VN[c] : Ng;

	const vertex_cu ca(mesh.VC[a].r, mesh.VC[a].g, mesh.VC[a].b);
	const vertex_cu cb(mesh.VC[b].r, mesh.VC[b].g, mesh.VC[b].b);
	const vertex_cu cc(mesh.VC[c].r, mesh.VC[c].g, mesh.VC[c].b);

	const vertex_cu color = ((1.f - u - v) * ca + u * cb + v * cc) / 255;

	const vertex_cu & light = *(vertex_cu *) optixLaunchParams.light;
	const vertex_cu position = (1.f - u - v) * A + u * B + v * C;
	const vertex_cu wi = normalize(light - position);
	const float dot_wi_normal = (wi, normal);

	vertex_cu & L = *getPRD<vertex_cu>();
	L = (dot_wi_normal < 0 ? -dot_wi_normal : dot_wi_normal) * color;

	bool occluded = true;
	uint32_t u0, u1;
	packPointer(&occluded, u0, u1);
	optixTrace( optixLaunchParams.traversable,
				position + 1e-5f * Ng,
				wi,
				1e-5f,		// tmin
				1e20f,		// tmax
				0.0f,		// rayTime
				OptixVisibilityMask(255),
				// For shadow rays: skip any/closest hit shaders and terminate on first
				// intersection with anything. The miss shader is used to mark if the
				// light was visible.
				OPTIX_RAY_FLAG_DISABLE_ANYHIT
				| OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT
				| OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT,
				1,	// SBT offset
				2,	// SBT stride
				1,	// missSBTIndex
				u0, u1);

	if(occluded) L = 0.6 * L;
}


extern "C" __global__ void __anyhit__radiance() {}

extern "C" __global__ void __anyhit__shadow() {}


extern "C" __global__ void __miss__radiance()
{
	bool & prd = *getPRD<bool>();
	prd = false;
}

extern "C" __global__ void __miss__shadow()
{
	bool & prd = *getPRD<bool>();
	prd = false;
}


extern "C" __global__ void __raygen__render_frame()
{
	const int ix = optixGetLaunchIndex().x;
	const int iy = optixGetLaunchIndex().y;

	const float sx = (float(ix) + .5f) / optixLaunchParams.frame.width;
	const float sy = (float(iy) + .5f) / optixLaunchParams.frame.height;

	vertex_cu & cam_pos = *(vertex_cu *) optixLaunchParams.cam_pos;

	vertex_cu ipv[3];
	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
		ipv[i][j] = optixLaunchParams.inv_proj_view[i + j * 4];

	vertex_cu d = { optixLaunchParams.inv_proj_view[0 * 4 + 3],
					optixLaunchParams.inv_proj_view[1 * 4 + 3],
					optixLaunchParams.inv_proj_view[2 * 4 + 3]
					};
	vertex_cu e = { optixLaunchParams.inv_proj_view[3 * 4 + 0],
					optixLaunchParams.inv_proj_view[3 * 4 + 1],
					optixLaunchParams.inv_proj_view[3 * 4 + 2]
					};

	float & de = optixLaunchParams.inv_proj_view[15];

	vertex_cu view = {sx * 2 - 1, sy * 2 - 1, 1};
	vertex_cu q = vertex_cu{(ipv[0], view), (ipv[1], view), (ipv[2], view)} + e;
	vertex_cu p = (1.f / ((d, view) + de)) * q;
	vertex_cu ray_dir = p - cam_pos;
	ray_dir /= *ray_dir;

	vertex_cu pixelColorPRD;
	uint32_t u0, u1;
	packPointer(&pixelColorPRD, u0, u1);
	optixTrace(	optixLaunchParams.traversable,
				cam_pos,
				ray_dir,
				0.f,	// tmin
				1e20f,	// tmax
				0.0f,	// rayTime
				OptixVisibilityMask(255),
				OPTIX_RAY_FLAG_DISABLE_ANYHIT, //OPTIX_RAY_FLAG_NONE,
				0,	// SBT offset
				2,	// SBT stride
				0,	// missSBTIndex
				u0, u1);

	const uint32_t fbIndex = ix + iy * optixLaunchParams.frame.width;

	float4 * frame = (float4 *) optixLaunchParams.frame.color_buffer;
	frame[fbIndex].x = pixelColorPRD.x;
	frame[fbIndex].y = pixelColorPRD.y;
	frame[fbIndex].z = pixelColorPRD.z;
	frame[fbIndex].w = 1;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

