#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_SPLAT_H
#define RT_EMBREE_SPLAT_H

#include "raytracing/rt_embree.h"


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


class embree_splat : public embree
{
	static const size_t K = 16;

	struct splat
	{
		glm::vec3 P[K];
		glm::vec3 N[K];
		glm::vec3 C[K];
		float radio;

		const glm::vec4 xyzr()
		{
			return glm::vec4(P[0], radio);
		}

		const glm::vec3 & normal()
		{
			return N[0];
		}

		const glm::vec3 & color()
		{
			return C[0];
		}

		void shading(const glm::vec3 & p, glm::vec3 & normal, glm::vec3 & color)
		{
			normal = glm::vec3(0);
			color = glm::vec3(0);

			float w, sum_w = 0, sigma = radio * 0.1;

			for(index_t i = 0; i < K; i++)
			{
				w = glm::length(p - P[i]);
				w = exp(-0.5 * w * w / (sigma * sigma)) / (sigma * sqrt(2 * M_PI));
				normal += w * N[i];
				color += w * C[i];
				sum_w += w;
			}
			
			normal /= sum_w;
			color /= sum_w;
		}
	};
	
	std::vector<splat> vsplat;

	public:
		embree_splat(const std::vector<che *> & meshes, const bool & pointcloud);
	
	private:
		index_t add_pointcloud(const che * mesh);
		float pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r);
		
		void init_splats(const che * mesh);
};


} // namespace gproshan

#endif // RT_EMBREE_SPLAT_H

#endif // GPROSHAN_EMBREE

