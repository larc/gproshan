#ifdef GPROSHAN_OPTIX_FAIL


#include "raytracing/rt_optix_params.h"


#include <optix_device.h>


// geometry processing and shape analysis framework                            
namespace gproshan::rt {


extern "C" __constant__ launch_params params;

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

extern "C" __global__ void __closesthit__shadow()
{
	/* not going to be used ... */
}

extern "C" __global__ void __closesthit__radiance()
{
	const TriangleMeshSBTData & sbtData = *(const TriangleMeshSBTData *) optixGetSbtDataPointer();

	// ------------------------------------------------------------------
	// gather some basic hit information
	// ------------------------------------------------------------------
	const int primID = optixGetPrimitiveIndex();
	const vec3i index = sbtData.index[primID];
	const float u = optixGetTriangleBarycentrics().x;
	const float v = optixGetTriangleBarycentrics().y;

	// ------------------------------------------------------------------
	// compute normal, using either shading normal (if avail), or
	// geometry normal (fallback)
	// ------------------------------------------------------------------
	const vertex_cu & A = sbtData.vertex[index.x];
	const vertex_cu & B = sbtData.vertex[index.y];
	const vertex_cu & C = sbtData.vertex[index.z];
	vertex_cu Ng = (B - A) * (C - A);
	vertex_cu Ns = (sbtData.normal)
		? ((1.f-u-v) * sbtData.normal[index.x]
			 +			 u * sbtData.normal[index.y]
			 +			 v * sbtData.normal[index.z])
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
	vertex_cu diffuseColor = sbtData.color;

	// ------------------------------------------------------------------
	// compute shadow
	// ------------------------------------------------------------------
	const vertex_cu surfPos
		= (1.f-u-v) * sbtData.vertex[index.x]
		+				 u * sbtData.vertex[index.y]
		+				 v * sbtData.vertex[index.z];
	const vertex_cu lightPos(-907.108f, 2205.875f, -400.0267f);
	const vertex_cu lightDir = lightPos - surfPos;

	// trace shadow ray:
	vertex_cu lightVisibility = 0.f;
	// the values we store the PRD pointer in:
	uint32_t u0, u1;
	packPointer(&lightVisibility, u0, u1);
	optixTrace(params.traversable,
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

	vertex_cu &prd = *(vertex_cu * ) getPRD<vertex_cu>();
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
	prd = vertex_cu(1.f);
}

extern "C" __global__ void __miss__shadow()
{
	// we didn't hit anything, so the light is visible
	vertex_cu &prd = *(vertex_cu *)getPRD<vertex_cu>();
	prd = vertex_cu(1.f);
}

//------------------------------------------------------------------------------
// ray gen program - the actual rendering happens in here
//------------------------------------------------------------------------------
extern "C" __global__ void __raygen__render_frame()
{
	// compute a test pattern based on pixel ID
	const int ix = optixGetLaunchIndex().x;
	const int iy = optixGetLaunchIndex().y;

	const auto &camera = params.camera;

	// our per-ray data for this example. what we initialize it to
	// won't matter, since this value will be overwritten by either
	// the miss or hit program, anyway
	vertex_cu pixelColorPRD = vertex_cu(0.f);

	// the values we store the PRD pointer in:
	uint32_t u0, u1;
	packPointer(&pixelColorPRD, u0, u1);

	// normalized screen plane position, in [0,1]^2
	const vec2f screen(vec2f(ix+.5f,iy+.5f)
										 / vec2f(params.frame.size));

	// generate ray direction
	vertex_cu rayDir = normalize(camera.direction
													 + (screen.x - 0.5f) * camera.horizontal
													 + (screen.y - 0.5f) * camera.vertical);

	optixTrace(params.traversable,
						 camera.position,
						 rayDir,
						 0.f,		// tmin
						 1e20f,	// tmax
						 0.0f,	 // rayTime
						 OptixVisibilityMask(255),
						 OPTIX_RAY_FLAG_DISABLE_ANYHIT,//OPTIX_RAY_FLAG_NONE,
						 0,						// SBT offset
						 2,							 // SBT stride
						 0,						// missSBTIndex
						 u0, u1);

	const int r = int(255.99f*pixelColorPRD.x);
	const int g = int(255.99f*pixelColorPRD.y);
	const int b = int(255.99f*pixelColorPRD.z);

	// convert to 32-bit rgba value (we explicitly set alpha to 0xff
	// to make stb_image_write happy ...
	const uint32_t rgba = 0xff000000
		| (r<<0) | (g<<8) | (b<<16);

	// and write to frame buffer ...
	const uint32_t fbIndex = ix+iy*params.frame.size.x;
	params.frame.colorBuffer[fbIndex] = rgba;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

