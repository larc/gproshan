#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_H
#define RT_EMBREE_H

#include "mesh/che.h"
#include "raytracing/raytracing.h"

#include <embree3/rtcore.h>
#include <glm/glm.hpp>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class embree : public raytracing
{
	protected:
		struct ray_hit : public RTCRayHit
		{
			ray_hit(const glm::vec3 & p_org = glm::vec3(0.0f),
					const glm::vec3 & v_dir = glm::vec3(0.0f),
					float near = 1e-5f,
					float far = FLT_MAX);

			glm::vec3 org() const;
			glm::vec3 dir() const;
			glm::vec3 color(const rt_mesh & mesh) const;
			glm::vec3 normal(const rt_mesh & mesh, const bool & flat = false) const;
			glm::vec3 position() const;
		};


		RTCDevice device;
		RTCScene scene;	
		RTCIntersectContext intersect_context;
		
	public:
		static float pc_radius;

	public:
		embree();
		embree(const std::vector<che *> & meshes, const bool & pointcloud = false);
		virtual ~embree();
	
	protected:
		bool intersect(ray_hit & r);
		bool occluded(ray_hit & r);

		void build_bvh(const std::vector<che *> & meshes, const bool & pointcloud = false);
		index_t add_sphere(const glm::vec4 & xyzr);
		index_t add_mesh(const che * mesh);
		virtual index_t add_pointcloud(const che * mesh);

		float pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r);

		glm::vec4 li(const glm::vec3 & light, const glm::vec3 & position, const glm::vec3 & normal, const glm::vec3 & color, const float & near = 1e-5f);
		glm::vec4 li(ray_hit r, const glm::vec3 & light, const bool & flat);
		
		glm::vec4 intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light, const bool & flat);
		float intersect_depth(const glm::vec3 & org, const glm::vec3 & dir);
};


} // namespace gproshan

#endif // RT_EMBREE_H

#endif // GPROSHAN_EMBREE

