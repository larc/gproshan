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

		bool closesthit_radiance(	vertex & color,
									vertex & attenuation,
									vertex & position,
									vertex & ray_dir,
									random<real_t> & rnd,
									const render_params & params,
									const bool & flat
									) const
		{
			ray_hit r(position, ray_dir);
			if(!intersect(r)) return false;

			eval_hit hit;
			splat_hit(hit, *g_meshes[r.hit.geomID], splats_pcs[r.hit.geomID], r.hit.primID, r.pos(), ray_dir, r.ray.tfar);

			color = eval_li(	hit, params.ambient, params.lights, params.n_lights, params.cam_pos,
								[&](const vec3 & position, const vec3 & wi, const float light_dist) -> bool
								{
									ray_hit ro((position - r.pos(), ray_dir) < 0 ? position : r.pos(), wi, 1e-3f, light_dist - 1e-3f);
									return occluded(ro);
								});

			return true;
		}
};


} // namespace gproshan

#endif // RT_SPLAT_EMBREE_H

