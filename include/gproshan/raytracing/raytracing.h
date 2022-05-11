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
		struct rt_mesh
		{
			che * mesh;
			bool pointcloud;

			che * operator -> () const
			{
				return mesh;
			}
		};

		std::map<index_t, rt_mesh> geomID_mesh;

		size_t n_samples = 0;

	public:
		raytracing() = default;
		virtual ~raytracing() = default;

		virtual void render(glm::vec4 * img,
							const glm::uvec2 & windows_size,
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

		virtual index_t cast_ray(	const glm::vec3 &,// org,
									const glm::vec3 &// dir
									) { return NIL; };
		
		virtual float intersect_depth(	const glm::vec3 &,// org,
										const glm::vec3 &// dir
										) { return 0; };

		virtual std::tuple<index_t, float> cast_ray_intersect_depth(	const glm::vec3 & origin,// org,
													const glm::vec3 & direction// dir
													) { return { cast_ray(origin, direction), intersect_depth( origin, direction) }; };
		
	protected:
		virtual glm::vec4 intersect_li(	const glm::vec3 &,// org,
										const glm::vec3 &,// dir,
										const glm::vec3 &,// light,
										const bool &// flat
										) { return glm::vec4(0); };

		
};


} // namespace gproshan

#endif // RAYTRACING_H

