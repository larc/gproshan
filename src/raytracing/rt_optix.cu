#ifdef GPROSHAN_OPTIX


#include "mesh/che.h"
#include "mesh/vertex.cuh"
#include "raytracing/rt_optix_params.h"


#include <optix_device.h>
#include <cuda_runtime.h>

#include <glm/glm.hpp>


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

//------------------------------------------------------------------------------
// closest hit and anyhit programs for radiance-type rays.
//
// Note eventually we will have to create one pair of those for each
// ray type and each geometry type we want to render; but this
// simple example doesn't use any actual geometries yet, so we only
// create a single, dummy, set of them (we do have to have at least
// one group of them to set up the SBT)
//------------------------------------------------------------------------------

extern "C" __global__ void __closesthit__shadow() {}

extern "C" __global__ void __closesthit__radiance()
{
	const CHE & sbtData = **(const CHE **) optixGetSbtDataPointer();

	// ------------------------------------------------------------------
	// gather some basic hit information
	// ------------------------------------------------------------------
	const int primID = optixGetPrimitiveIndex();
	const index_t he = primID * che::mtrig;
	const float u = optixGetTriangleBarycentrics().x;
	const float v = optixGetTriangleBarycentrics().y;

	// ------------------------------------------------------------------
	// compute normal, using either shading normal (if avail), or
	// geometry normal (fallback)
	// ------------------------------------------------------------------
	const vertex_cu & A = sbtData.GT[sbtData.VT[he]];
	const vertex_cu & B = sbtData.GT[sbtData.VT[he + 1]];
	const vertex_cu & C = sbtData.GT[sbtData.VT[he + 2]];

	vertex_cu Ng = (B - A) * (C - A);
	vertex_cu Ns = (sbtData.VN)
		? ((1.f - u - v) * sbtData.VN[sbtData.VT[he]]
			 +			 u * sbtData.VN[sbtData.VT[he + 1]]
			 +			 v * sbtData.VN[sbtData.VT[he + 2]])
		: Ng;

	// ------------------------------------------------------------------
	// face-forward and normalize normals
	// ------------------------------------------------------------------
	const vertex_cu rayDir = (vertex_cu) optixGetWorldRayDirection();

	if((rayDir , Ng) > 0.f) Ng = -Ng;
	Ng /= *Ng;

	if((Ng , Ns) < 0.f)
		Ns -= 2.f * (Ng , Ns) * Ng;
	Ns /= *Ns;

	// ------------------------------------------------------------------
	// compute diffuse material color, including diffuse texture, if
	// available
	// ------------------------------------------------------------------
	vertex_cu diffuseColor(230.0/255, 240.0/255, 250.0/255);

	// ------------------------------------------------------------------
	// compute shadow
	// ------------------------------------------------------------------
	const vertex_cu surfPos = (1.f - u - v) * A + u * B + v * C;
	const vertex_cu lightPos = { optixLaunchParams.light[0], optixLaunchParams.light[1], optixLaunchParams.light[2] };
	const vertex_cu lightDir = lightPos - surfPos;

	// trace shadow ray:
	vertex_cu lightVisibility = 0.f;
	// the values we store the PRD pointer in:
	uint32_t u0, u1;
	packPointer(&lightVisibility, u0, u1);
	optixTrace(optixLaunchParams.traversable,
						 surfPos + 1e-3f * Ng,
						 lightDir,
						 1e-3f,			// tmin
						 1.f-1e-3f,	// tmax
						 0.0f,			 // rayTime
						 OptixVisibilityMask(255),
						 // For shadow rays: skip any/closest hit shaders and terminate on first
						 // intersection with anything. The miss shader is used to mark if the
						 // light was visible.
						 OPTIX_RAY_FLAG_DISABLE_ANYHIT
						 | OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT
						 | OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT,
						 1,						// SBT offset
						 2,							 // SBT stride
						 1,						// missSBTIndex
						 u0, u1);

	// ------------------------------------------------------------------
	// final shading: a bit of ambient, a bit of directional ambient,
	// and directional component based on shadowing
	// ------------------------------------------------------------------
	const float cosDN = 0.1f + .8f * fabsf((rayDir , Ns));

	vertex_cu & prd = *(vertex_cu * ) getPRD<vertex_cu>();
	prd = (.1f + (.2f + .8f * lightVisibility) * cosDN) * diffuseColor;
}

extern "C" __global__ void __anyhit__radiance() {}

extern "C" __global__ void __anyhit__shadow() {}

//------------------------------------------------------------------------------
// miss program that gets called for any ray that did not have a
// valid intersection
//
// as with the anyhit/closest hit programs, in this example we only
// need to have _some_ dummy function to set up a valid SBT
// ------------------------------------------------------------------------------

extern "C" __global__ void __miss__radiance()
{
	vertex_cu &prd = *(vertex_cu *) getPRD<vertex_cu>();
	// set to constant white as background color
//	prd = {0, 0, 0};
}

extern "C" __global__ void __miss__shadow()
{
	// we didn't hit anything, so the light is visible
	vertex_cu &prd = *(vertex_cu *) getPRD<vertex_cu>();
//	prd = {1, 1, 1};
}

//------------------------------------------------------------------------------
// ray gen program - the actual rendering happens in here
//------------------------------------------------------------------------------
extern "C" __global__ void __raygen__render_frame()
{
	const int ix = optixGetLaunchIndex().x;
	const int iy = optixGetLaunchIndex().y;

	vertex_cu pixelColorPRD;

	uint32_t u0, u1;
	packPointer(&pixelColorPRD, u0, u1);

	glm::vec2 screen = glm::vec2(	(float(ix) + .5f) / optixLaunchParams.frame.width,
									(float(iy) + .5f) / optixLaunchParams.frame.height
									);

	glm::vec3 cam_pos = glm::vec3(	optixLaunchParams.cam_pos[0],
									optixLaunchParams.cam_pos[1],
									optixLaunchParams.cam_pos[2]
									);
	glm::mat4 inv_proj_view;
	for(int i = 0; i < 4; ++i)
	for(int j = 0; j < 4; ++j)
		//if(ix + iy == 0)printf("m %f\n", optixLaunchParams.inv_proj_view[i * 4 + j]);
		inv_proj_view[i][j] = optixLaunchParams.inv_proj_view[i * 4 + j];

	glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
	glm::vec4 q = inv_proj_view * view;
	glm::vec3 p = glm::vec3(q * (1.f / q.w));
	glm::vec3 r = p - cam_pos;

	vertex_cu position = {cam_pos.x, cam_pos.y, cam_pos.z}; 
	vertex_cu ray_dir = {r.x, r.y, r.z};
	ray_dir /= *ray_dir;

	if(ix + iy == 0)
		printf("%f %f %f\n", position.x, position.y, position.z);

	optixTrace(	optixLaunchParams.traversable,
				position,
				ray_dir,
				0.f,	// tmin
				1e20f,	// tmax
				0.0f,	 // rayTime
				OptixVisibilityMask(255),
				OPTIX_RAY_FLAG_DISABLE_ANYHIT, //OPTIX_RAY_FLAG_NONE,
				0,						// SBT offset
				2,							 // SBT stride
				0,						// missSBTIndex
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

