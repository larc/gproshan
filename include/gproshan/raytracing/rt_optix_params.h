#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_PARAMS_H
#define RT_OPTIX_PARAMS_H


// geometry processing and shape analysis framework
namespace gproshan::rt {


struct launch_params
{
};


struct vertex_cu;

struct mesh_sbt_data
{
	vertex_cu * vertex;
	vertex_cu * normal;
	vertex_cu * color;
	vertex_cu * index;
};


} // namespace gproshan

#endif // RT_OPTIX_PARAMS_H

#endif // GPROSHAN_OPTIX

