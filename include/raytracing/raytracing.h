#ifndef RAYTRACING_H
#define RAYTRACING_H

#include "mesh/che.h"

#include <vector>
#include <map>

#include <glm/glm.hpp>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class raytracing
{
	protected:
		std::map<index_t, const che *> geomID_mesh;

		size_t width;
		size_t height;
		size_t n_samples;
	
	public:
		glm::vec4 * img;

	public:
		raytracing();
		virtual ~raytracing();

		virtual bool rt_restart(const size_t & w, const size_t & h);
		virtual void pathtracing(	const glm::uvec2 & windows_size,
							const glm::mat4 & view_mat,
							const glm::mat4 & proj_mat,
							const std::vector<glm::vec3> & light,
							const bool & flat,
							const bool & restart = false
							);

		virtual float * raycaster(	const glm::uvec2 & windows_size,
									const glm::mat4 & view_mat,
									const glm::mat4 & proj_mat,
									const index_t & samples = 4
									);

	
	protected:
		virtual glm::vec4 intersect_li(	const glm::vec3 & org,
										const glm::vec3 & dir,
										const glm::vec3 & light,
										const bool & flat ) = 0;
		
		virtual float intersect_depth(	const glm::vec3 & org,
										const glm::vec3 & dir ) = 0;
};


} // namespace gproshan

#endif // RAYTRACING_H

