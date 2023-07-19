#ifndef RT_SPLAT_OPTIX_H
#define RT_SPLAT_OPTIX_H

#include <gproshan/raytracing/splat.h>
#include <gproshan/raytracing/optix.h>


#ifdef GPROSHAN_OPTIX


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat_optix: public splat, public optix
{
	private:
		splats_data ** d_splats_pcs = nullptr;

	public:
		splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats);
		~splat_optix();
};


} // namespace gproshan

#endif // GPROSHAN_OPTIX

#endif // RT_SPLAT_OPTIX_H

