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
	const float2 bar = optixGetTriangleBarycentrics();

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
	vec3 & attenuation	= trace[1];
	vec3 & position		= trace[2];
	vec3 & ray_dir		= trace[3];

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

	random<float> rnd = optixGetPayload_2();
	color *= attenuation;
	position = hit.position;

	if(!hit.scatter_mat(ray_dir, rnd))
		attenuation = 0;

	attenuation /= 2;
	optixSetPayload_2(rnd);
}


extern "C" __global__ void __anyhit__radiance() {}

extern "C" __global__ void __anyhit__shadow() {}


extern "C" __global__ void __miss__radiance()
{
	optixSetPayload_0(0);
}

extern "C" __global__ void __miss__shadow()
{
	optixSetPayload_0(0);
}


extern "C" __global__ void __raygen__render_frame()
{
	const uvec2 & id = {optixGetLaunchIndex().x,
						optixGetLaunchIndex().y
						};

	const uvec2 & pos = id + optix_params.viewport_pos;

	random<float> rnd(pos.x() + optix_params.window_size.x() * pos.y(), optix_params.n_frames);

	unsigned int depth = optix_params.depth;
	unsigned int samples = optix_params.n_samples;

	vec3 color_acc = 0;

	vec3 trace[4];
	vec3 & color		= trace[0];
	vec3 & attenuation	= trace[1];
	vec3 & position		= trace[2];
	vec3 & ray_dir		= trace[3];

	uint32_t u0, u1;

	do
	{
		color		= 0;
		attenuation = 1;
		position	= optix_params.cam_pos;
		ray_dir		= ray_view_dir(pos, optix_params.window_size, optix_params.inv_proj_view, optix_params.cam_pos, rnd);

		depth = optix_params.depth;
		do
		{
			pack_pointer(trace, u0, u1);
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
						u0, u1, (unsigned int &) rnd);

			if(!u0) break;	// miss
			color_acc += color;
		}
		while(--depth);
	}
	while(--samples);

	color_acc /= optix_params.n_samples;

	vec4 & pixel_color = optix_params.color_buffer[id.x() + id.y() * optixGetLaunchDimensions().x];
	pixel_color = (pixel_color * optix_params.n_frames + (color_acc, 1)) / (optix_params.n_frames + 1);
}


} // namespace gproshan

