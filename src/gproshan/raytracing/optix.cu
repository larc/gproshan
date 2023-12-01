#include <gproshan/mesh/che.cuh>
#include <gproshan/raytracing/utils.h>
#include <gproshan/raytracing/optix_params.h>

#include <optix_device.h>
#include <cuda_runtime.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


extern "C" __constant__ launch_params optix_params;

static __forceinline__ __device__
void * unpack_pointer(uint32_t i0, uint32_t i1)
{
	return (void *) (uint64_t(i0) << 32 | i1);
}

static __forceinline__ __device__
void pack_pointer(void * ptr, uint32_t & i0, uint32_t & i1)
{
	const uint64_t uptr = uint64_t(ptr);
	i0 = uptr >> 32;
	i1 = uptr & 0x00000000ffffffff;
}

template<typename T>
static __forceinline__ __device__
T * ray_data()
{
	return (T *) unpack_pointer(optixGetPayload_0(), optixGetPayload_1());
}


extern "C" __global__ void __closesthit__shadow() {}

extern "C" __global__ void __closesthit__radiance()
{
	const CHE & mesh = **(const CHE **) optixGetSbtDataPointer();

	const int primID = optixGetPrimitiveIndex();
	float2 bar = optixGetTriangleBarycentrics();

	OptixTraversableHandle gas = optixGetGASTraversableHandle();
	const index_t sbtID = optixGetSbtGASIndex();
	const float time = optixGetRayTime();

	vertex data[3];
	optixGetTriangleVertexData(gas, primID, sbtID, time, (float3 *) data);

	const vertex & A = data[0];
	const vertex & B = data[1];
	const vertex & C = data[2];

	eval_hit hit(mesh, primID, bar.x, bar.y, optix_params.sc);
	hit.normal = optix_params.flat ? normalize(cross(B - A, C - A)) : hit.normal;
	hit.position = (1.f - hit.u - hit.v) * A + hit.u * B + hit.v * C;

	vec3 * trace = ray_data<vec3>();
	vec3 & color		= trace[0];
	vec3 & position		= trace[1];
	vec3 & scattered	= trace[2];	// in ray_dir / out scattered
	vec3 & attenuation	= trace[3];

	color = eval_li(hit, optix_params.ambient, optix_params.lights, optix_params.n_lights, optix_params.cam_pos,
					[&](const vec3 & position, const vec3 & wi, const float & light_dist) -> bool
					{
						uint32_t occluded = 1;
						optixTrace( optix_params.traversable,
									* (float3 *) &position,
									* (float3 *) &wi,
									1e-3f,					// tmin
									light_dist - 1e-3f,		// tmax
									0.0f,					// rayTime
									OptixVisibilityMask(255),
										OPTIX_RAY_FLAG_DISABLE_ANYHIT
										| OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT
										| OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT,
										1,	// SBT offset
										2,	// SBT stride
										1,	// missSBTIndex
										occluded);

							return occluded != 0;
						});

	color *= attenuation;
	position = hit.position;

	if(!hit.scatter_mat(scattered, scattered))
		attenuation = 0;
}


extern "C" __global__ void __anyhit__radiance() {}

extern "C" __global__ void __anyhit__shadow() {}


extern "C" __global__ void __miss__radiance()
{
	vec4 & pixel_color = *ray_data<vec4>();
	pixel_color = {0, 0, 0, 0};
}

extern "C" __global__ void __miss__shadow()
{
	optixSetPayload_0(0);
}


extern "C" __global__ void __raygen__render_frame()
{
	const int ix = optixGetLaunchIndex().x;
	const int iy = optixGetLaunchIndex().y;

	vec3 trace[4];
	vec3 & color		= trace[0];
	vec3 & position		= trace[1];
	vec3 & ray_dir		= trace[2];
	vec3 & attenuation	= trace[3] = 1;

	position = optix_params.cam_pos;
	ray_dir = ray_view_dir(	{ix + optix_params.viewport_x, iy + optix_params.viewport_y},
							{optix_params.window_width, optix_params.window_height},
							optix_params.inv_proj_view,
							optix_params.cam_pos
							);

	uint32_t u0, u1;
	pack_pointer(trace, u0, u1);

	int depth = 3;
	vec3 color_acc = 0;

	while(--depth)
	{
		optixTrace(	optix_params.traversable,
					* (float3 *) &position,
					* (float3 *) &ray_dir,
					1e-5f,	// tmin
					1e20f,	// tmax
					0.0f,	// rayTime
					OptixVisibilityMask(255),
					OPTIX_RAY_FLAG_DISABLE_ANYHIT, //OPTIX_RAY_FLAG_NONE,
					0,	// SBT offset
					2,	// SBT stride
					0,	// missSBTIndex
					u0, u1);

		color_acc += color;
	}

	vec4 & pixel_color = optix_params.color_buffer[ix + iy * optixGetLaunchDimensions().x];
	pixel_color = (pixel_color * optix_params.n_samples + (color_acc, 1)) / (optix_params.n_samples + 1);
}


} // namespace gproshan

