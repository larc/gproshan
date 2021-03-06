#ifdef GPROSHAN_OPTIX

#include <optix_device.h>

extern "C" __constant__ void * optix_launch_params;

extern "C" __global__ void closest_hit()
{
	const int f = optixGetPrimitiveIndex();
}

extern "C" __global__ void any_hit()
{
}

extern "C" __global__ void miss_hit()
{
}

extern "C" __global__ void __raygen__rt()
{
}

#endif // GPROSHAN_OPTIX

