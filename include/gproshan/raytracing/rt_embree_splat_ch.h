#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_SPLAT_CH_H
#define RT_EMBREE_SPLAT_CH_H

#include "raytracing/rt_embree.h"
#include "geometry/convex_hull.h"


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {

inline glm::vec3 glm_vec3(const vertex & v)
{
	return glm::vec3(v.x, v.y, v.z);
}

class embree_splat_ch : public embree
{
	static const size_t K = 32;

	struct splat
	{
		std::vector<index_t> points;

		operator std::vector<index_t> & ()
		{
			return points;
		}

		float shading(const rt_mesh & mesh, const glm::vec3 & p, glm::vec3 & normal, glm::vec3 & color)
		{
			normal = glm::vec3(0);
			color = glm::vec3(0);

			float w, sum_w = 0, sigma = pc_radius;

			for(index_t & v: points)
			{
				w = glm::length(p - glm_vec3(mesh->gt(v)));
				w = exp(-0.5 * w * w / (sigma * sigma));
				normal += w * glm_vec3(mesh->normal(v));
				color += w * glm_vec3(mesh->color(v));
				sum_w += w;
			}

			normal /= sum_w;
			color /= sum_w;

			return sum_w;
		}
	};

	std::vector<splat> vsplat;
	che * ch_mesh = nullptr;

	public:
		embree_splat_ch(const std::vector<che *> & meshes, const bool & pointcloud);

	private:
		index_t add_pointcloud(const che * mesh);
		float pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r);

		void init_splats(const che * mesh);
};


} // namespace gproshan

#endif // RT_EMBREE_SPLAT_CH_H

#endif // GPROSHAN_EMBREE

