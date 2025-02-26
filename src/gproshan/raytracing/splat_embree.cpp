#include <gproshan/raytracing/splat_embree.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat_embree::splat_embree(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats)
{
	std::vector<const che *> pcs;
	pcs.reserve(size(pointclouds));
	for(che * pc: pointclouds)
		pcs.push_back(pc);

	build_bvh(pcs, {mat4::identity()});
}

bool splat_embree::closesthit_radiance(	vertex & color,
										vertex & attenuation,
										vertex & position,
										vertex & ray_dir,
										float & dist,
										random<float> & rnd,
										const render_params & params,
										const bool
										) const
{
	ray_hit r(position, ray_dir);
	if(!intersect(r)) return false;

	dist += r.ray.tfar;

	eval_hit hit;
	const float w = splat_hit(hit, *g_meshes[r.hit.geomID], splats_pcs[r.hit.geomID], r.hit.primID, r.pos(), ray_dir, dist);

	color = eval_li(	hit, params.ambient, params.lights, params.n_lights, params.cam_pos,
						[&](const vec3 & position, const vec3 & wi, const float light_dist) -> bool
						{
							ray_hit ro((position - r.pos(), ray_dir) < 0 ? position : r.pos(), wi, 1e-3f, light_dist - 1e-3f);
							return occluded(ro);
						});

	color *= attenuation;
	position = r.pos();

	if(w < 1e-3f)
	{
		color *= w;
		return true;
	}

	if(!hit.scatter_diffuse(ray_dir, rnd))
		attenuation = 0;

	attenuation /= 2;

	return true;
}


} // namespace gproshan

