#ifdef GPROSHAN_OPTIX

#ifndef RT_OPTIX_H
#define RT_OPTIX_H

#include "che.h"
#include "raytracing.h"

#include <cuda_runtime.h>
#include <optix.h>
#include <optix_stubs.h>

#include <glm/glm.hpp>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class optix : public raytracing
{
	public:
		optix();
		~optix();

		void build_bvh();
		index_t add_sphere(const glm::vec4 & xyzr);
		index_t add_mesh(const che * mesh, const glm::mat4 & model_matrix = glm::mat4(1.f));
		index_t add_point_cloud(const che * mesh, const glm::mat4 & model_matrix = glm::mat4(1.f));

	private:

		const glm::vec4 intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light);
		const float intersect_depth(const glm::vec3 & org, const glm::vec3 & dir);
};


} // namespace gproshan

#endif // RT_OPTIX_H

#endif // GPROSHAN_OPTIX

