#include <gproshan/mesh/che.h>
#include <gproshan/mesh/vertex.cuh>
#include <gproshan/raytracing/rt_optix_params.h>


#include <optix_device.h>
#include <cuda_runtime.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


static __forceinline__ __device__
vertex_cu operator * (const real_t & a, const vertex_cu & v)
{
	return vertex_cu(a * v.x, a * v.y, a * v.z);
}

static __forceinline__ __device__
vertex_cu normalize (const vertex_cu & v)
{
	return v / *v;
}


extern "C" __constant__ launch_params optix_params;

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

	const unsigned int primID = optixGetPrimitiveIndex();

	const int he = primID * che::mtrig;
	const float u = optixGetTriangleBarycentrics().x;
	const float v = optixGetTriangleBarycentrics().y;

	const int a = mesh.VT[he];
	const int b = mesh.VT[he + 1];
	const int c = mesh.VT[he + 2];

	OptixTraversableHandle gas = optixGetGASTraversableHandle();
	const unsigned int sbtID = optixGetSbtGASIndex();
	const float time = optixGetRayTime();

	vertex_cu data[3];
	optixGetTriangleVertexData(gas, primID, sbtID, time, (float3 *) data);

	const vertex_cu & A = data[0];
	const vertex_cu & B = data[1];
	const vertex_cu & C = data[2];

	const vertex_cu Ng = normalize((B - A) * (C - A));
	const vertex_cu normal = optix_params.flat ? Ng : (1.f - u - v) * mesh.VN[a] + u * mesh.VN[b] + v * mesh.VN[c];

	const vertex_cu ca(mesh.VC[a].r, mesh.VC[a].g, mesh.VC[a].b);
	const vertex_cu cb(mesh.VC[b].r, mesh.VC[b].g, mesh.VC[b].b);
	const vertex_cu cc(mesh.VC[c].r, mesh.VC[c].g, mesh.VC[c].b);

	const vertex_cu & light = *(vertex_cu *) optix_params.light;
	const vertex_cu color = ((1.f - u - v) * ca + u * cb + v * cc) / 255;
	const vertex_cu position = (1.f - u - v) * A + u * B + v * C;

	const vertex_cu wi = normalize(light - position);
	const float dot_wi_normal = (wi, normal);

	vertex_cu & L = *getPRD<vertex_cu>();
	L = (dot_wi_normal < 0 ? -dot_wi_normal : dot_wi_normal) * color;

	unsigned int occluded = 1;
	optixTrace( optix_params.traversable,
				position,
				wi,
				1e-5f,		// tmin
				1e20f,		// tmax
				0.0f,		// rayTime
				OptixVisibilityMask(255),
				OPTIX_RAY_FLAG_DISABLE_ANYHIT
				| OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT
				| OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
				1,	// SBT offset
				2,	// SBT stride
				1,	// missSBTIndex
				occluded);

	if(occluded) L = 0.4 * L;
}


extern "C" __global__ void __anyhit__radiance() {}

extern "C" __global__ void __anyhit__shadow() {}


extern "C" __global__ void __miss__radiance()
{
	vertex_cu & prd = *getPRD<vertex_cu>();
	prd = {0, 0, 0};
}

extern "C" __global__ void __miss__shadow()
{
	optixSetPayload_0(0);
}


extern "C" __global__ void __raygen__render_frame()
{
	const int ix = optixGetLaunchIndex().x;
	const int iy = optixGetLaunchIndex().y;

	const float sx = (float(ix + optix_params.viewport_x) + .5f) / optix_params.window_width;
	const float sy = (float(iy + optix_params.viewport_y) + .5f) / optix_params.window_height;

	vertex_cu & cam_pos = *(vertex_cu *) optix_params.cam_pos;

	vertex_cu ipv[3];
	for(int i = 0; i < 3; ++i)
	for(int j = 0; j < 3; ++j)
		ipv[i][j] = optix_params.inv_proj_view[i + j * 4];

	vertex_cu d = { optix_params.inv_proj_view[0 * 4 + 3],
					optix_params.inv_proj_view[1 * 4 + 3],
					optix_params.inv_proj_view[2 * 4 + 3]
					};
	vertex_cu e = { optix_params.inv_proj_view[3 * 4 + 0],
					optix_params.inv_proj_view[3 * 4 + 1],
					optix_params.inv_proj_view[3 * 4 + 2]
					};

	float & de = optix_params.inv_proj_view[15];

	vertex_cu view = {sx * 2 - 1, sy * 2 - 1, 1};
	vertex_cu q = vertex_cu{(ipv[0], view), (ipv[1], view), (ipv[2], view)} + e;
	vertex_cu p = (1.f / ((d, view) + de)) * q;
	vertex_cu ray_dir = p - cam_pos;
	ray_dir /= *ray_dir;

	vertex_cu pixelColorPRD;
	uint32_t u0, u1;
	packPointer(&pixelColorPRD, u0, u1);
	optixTrace(	optix_params.traversable,
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

	const uint32_t fbIndex = ix + iy * optix_params.viewport_width;

	float4 * frame = (float4 *) optix_params.color_buffer;
	frame[fbIndex].x = pixelColorPRD.x;
	frame[fbIndex].y = pixelColorPRD.y;
	frame[fbIndex].z = pixelColorPRD.z;
	frame[fbIndex].w = 1;
}


} // namespace gproshan

