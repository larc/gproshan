#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_H
#define RT_EMBREE_H

#include "che.h"
#include "raytracing.h"

#include <embree3/rtcore.h>
#include <glm/glm.hpp>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class embree : public raytracing
{
	struct ray_hit: public RTCRayHit
	{
		ray_hit(const glm::vec3 & origin = glm::vec3(0.0f),
				const glm::vec3 & direction = glm::vec3(0.0f),
				float near = 0.0f,
				float far = FLT_MAX)
		{
			ray.org_x = origin.x;
			ray.org_y = origin.y;
			ray.org_z = origin.z;
			ray.tnear = near;

			ray.dir_x = direction.x;
			ray.dir_y = direction.y;
			ray.dir_z = direction.z;
			ray.time = 0.0f;

			ray.tfar = far;
			ray.mask = 0;
			ray.flags = 0;

			//hit.primID = RTC_INVALID_GEOMETRY_ID;
			hit.geomID = RTC_INVALID_GEOMETRY_ID;
			hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
		}

		const glm::vec3 org() const
		{
			return {ray.org_x, ray.org_y, ray.org_z};
		}
		
		const glm::vec3 dir() const
		{
			return {ray.dir_x, ray.dir_y, ray.dir_z};
		}
		
		const glm::vec3 geometry_normal() const
		{
			return glm::normalize(glm::vec3(hit.Ng_x, hit.Ng_y, hit.Ng_z));
		}

		const glm::vec3 shading_normal(const che * mesh) const
		{
			vertex n = mesh->shading_normal(hit.primID, 1.0 - hit.u - hit.v, hit.u, hit.v);

			return glm::normalize(glm::vec3(n.x, n.y, n.z));
		}

		const glm::vec3 position() const
		{
			return org() + ray.tfar * dir();
		}
	};


	RTCDevice device;
	RTCScene scene;	
	RTCIntersectContext intersect_context;

	public:
		embree(const std::vector<che *> & meshes);
		~embree();
	
	private:
		bool intersect(ray_hit & r);
		bool occluded(ray_hit & r);

		glm::vec4 intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light, const bool & flat);
		float intersect_depth(const glm::vec3 & org, const glm::vec3 & dir);
		
		void build_bvh(const std::vector<che *> & meshes);
		index_t add_sphere(const glm::vec4 & xyzr);
		index_t add_mesh(const che * mesh);
		index_t add_point_cloud(const che * mesh);

		glm::vec4 li(const ray_hit & r, const glm::vec3 & light, const bool & flat);
};


} // namespace gproshan

#endif // RT_EMBREE_H

#endif // GPROSHAN_EMBREE

