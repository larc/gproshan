#ifdef GPROSHAN_EMBREE

#ifndef RT_EMBREE_SPLAT_CH_H
#define RT_EMBREE_SPLAT_CH_H

#include "raytracing/rt_embree.h"
#include "geometry/convex_hull.h"

#include <algorithm>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


unsigned int expand_bits(unsigned int v);
unsigned int morton_2d(float x, float y);

class embree_splat_ch : public embree
{
	struct splat
	{
		struct ipoint_code	// index point and 2d morton code
		{
			index_t p, code;
			bool operator < (const ipoint_code & ipc) const
			{
				return code < ipc.code;
			}
		};

		std::vector<ipoint_code> ipoints;
		vertex c, t, b, n;				// center, tbn matrix

		index_t & operator [] (const index_t & i)
		{
			return ipoints[i].p;
		}

		index_t & code(const index_t & i)
		{
			return ipoints[i].code;
		}

		size_t size() const
		{
			return ipoints.size();
		}

		void push_back(const index_t & p)
		{
			ipoints.push_back({p, 0});
		}

		float shading(const rt_mesh & mesh, const glm::vec3 & p, glm::vec3 & normal, glm::vec3 & color)
		{
			normal = glm::vec3(0);
			color = glm::vec3(0);

			float w, sum_w = 0;
			float & sigma = pc_radius;

			vertex h(p.x, p.y, p.z);
			to2d(h);
			int k = std::lower_bound(ipoints.begin(), ipoints.end(),
										ipoint_code{0, morton_2d((h.x + 1) / 2, (h.y + 1) / 2)}) - ipoints.begin();

			const int nk = 8;
			int begin = std::max(k - nk, 0);
			int end = std::min(k + nk, (int) ipoints.size());
			for(int i = begin; i < end; ++i)
//			for(auto & i: ipoints)
			{
				const index_t & v = ipoints[i].p;
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

		void to2d(vertex & v) const
		{
			v -= c;
			v = {(t, v), (b, v), (n, b)};
		}

		void to3d(vertex & v) const
		{
			v = vertex{	(vertex{t.x, b.x, n.x}, v),
						(vertex{t.y, b.y, n.y}, v),
						(vertex{t.z, b.z, n.z}, v)
						} + c;
		}
	};

	public:
		static float r_threshold;
		static float n_threshold;
		static size_t max_neighbors;

	private:
		std::vector<splat> vsplat;
		std::vector<float> csplat;
		std::vector<index_t> primID_splat;

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

