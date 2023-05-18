#ifndef RT_SPLAT_EMBREE_H
#define RT_SPLAT_EMBREE_H

#include <gproshan/raytracing/splat.h>
#include <gproshan/raytracing/embree.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat_embree: public splat, public embree
{
	public:
		splat_embree(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats)
		{
			build_bvh(pointclouds, {mat4::identity()});
		}

		vec3 closesthit_radiance(const vertex & org, const vertex & dir, const vertex * lights, const int & n_lights, const vertex & cam_pos, const bool & flat)
		{
			ray_hit r(org, dir);
			if(!intersect(r)) return {};

			eval_hit hit;
			splat_hit(hit, splats_pcs[r.hit.geomID], r.hit.primID, r.pos());
			
			return eval_li(	hit, lights, n_lights, cam_pos,
							[&](const vec3 & position, const vec3 & wi, const float & light_dist) -> bool
							{
								ray_hit ro(position, wi, 1e-3f, light_dist - 1e-3f);
								return occluded(ro);
							});
		}
};


} // namespace gproshan

#endif // RT_SPLAT_EMBREE_H

