#ifndef RT_EMBREE_SPLAT_CH_H
#define RT_EMBREE_SPLAT_CH_H

#include <gproshan/raytracing/rt_embree.h>
#include <gproshan/geometry/convex_hull.h>

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
		vertex c;	// center
		mat3 tbn;	// tbn matrix

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

		float shading(const rt_mesh & mesh, const vec3 & p, vec3 & normal, vec3 & color)
		{
			normal = vec3(0);
			color = vec3(0);

			vertex h = p;
			to2d(h);
			int k = std::lower_bound(ipoints.begin(), ipoints.end(), ipoint_code{0, morton_2d((h.x() + 1) / 2, (h.y() + 1) / 2)}) - ipoints.begin();

			float w, sum_w = 0;
			float sigma = k < ipoints.size() ? length(p - mesh->point(ipoints[k].p)) :
											length(p - mesh->point(ipoints[k - 1].p));
			int begin = std::max(k - k_neighbors, 0);
			int end = std::min(k + k_neighbors, (int) ipoints.size());

			sigma /= 2;
			for(int i = begin; i < end; ++i)
			{
				const index_t & v = ipoints[i].p;
				w = length(p - mesh->point(v));
				w = exp(-0.5 * w * w / (sigma * sigma));
				normal += w * mesh->normal(v);
				color += w * mesh->color(v);
				sum_w += w;
			}

			normal /= sum_w;
			color /= sum_w;

			return sum_w;
		}

		void to2d(vertex & v) const
		{
			v -= c;
			v = mat3::transpose(tbn) * v;
		}

		void to3d(vertex & v) const
		{
			v = tbn * v + c;
		}
	};

	public:
		static bool show_chsplats;
		static int k_neighbors;
		static float r_threshold;
		static float n_threshold;
		static size_t max_neighbors;

	private:
		std::vector<splat> vsplat;
		std::vector<float> csplat;
		std::vector<index_t> primID_splat;

	public:
		embree_splat_ch(const std::vector<che *> & meshes,
						const std::vector<mat4> & model_mats
						);

	private:
		index_t add_pointcloud(const che * mesh, const mat4 & model_mat);
		float pointcloud_hit(vec3 & position, vec3 & normal, vec3 & color, ray_hit r);

		void init_splats(const che * mesh);
};


} // namespace gproshan

#endif // RT_EMBREE_SPLAT_CH_H

