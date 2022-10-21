#ifndef RT_EMBREE_H
#define RT_EMBREE_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/raytracing.h>

#include <cfloat>

#include <embree3/rtcore.h>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class embree : public raytracing
{
	protected:
		struct ray_hit : public RTCRayHit
		{
			ray_hit(const vertex & p_org = {0, 0, 0},
					const vertex & v_dir = {0, 0, 0},
					float near = 1e-5f,
					float far = FLT_MAX);

			vertex org() const;
			vertex dir() const;
			vertex position() const;
			vertex normal() const;
		};


		RTCDevice device;
		RTCScene scene;
		RTCIntersectContext intersect_context;

	protected:
		float pc_radius = 1;

	public:
		embree();
		embree(	const std::vector<che *> & meshes,
				const std::vector<mat4> & model_mats,
				const bool & pointcloud = false,
				const float & pcr = 0.01
				);
		virtual ~embree();

		virtual index_t closest_vertex(const vertex & org, const vertex & dir);
		virtual eval_hit<float> intersect(const vertex & org, const vertex & dir);


	protected:
		void build_bvh(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats, const bool & pointcloud = false);
		index_t add_sphere(const vec4 & xyzr);
		index_t add_mesh(const che * mesh, const mat4 & model_mat);

		virtual index_t add_pointcloud(const che * mesh, const mat4 & model_mat);

		virtual vec3 closesthit_radiance(const vertex & org, const vertex & dir, const vertex * lights, const int & n_lights, const bool & flat);

		float intersect_depth(const vertex & org, const vertex & dir);

		bool intersect(ray_hit & r);
		bool occluded(ray_hit & r);
};


} // namespace gproshan

#endif // RT_EMBREE_H

