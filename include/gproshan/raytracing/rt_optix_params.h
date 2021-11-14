#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


#include "geodesics/vertex.cuh"


// geometry processing and shape analysis framework                            
namespace gproshan::rt {


struct launch_params
{
};


struct TriangleMeshSBTData
{
	vertex_cu color;
	vertex_cu * vertex;
	vertex_cu * normal;
};


} // namespace gproshan

#endif // RT_OPTIX_PARAMS_H

#endif // GPROSHAN_OPTIX

