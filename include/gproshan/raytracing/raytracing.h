#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/render_params.h>

#include <vector>
#include <map>
#include <random>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


struct hit
{
	index_t idx = NIL;
	real_t dist = INFINITY;
	vertex color;
	vertex normal;
};


class raytracing
{
	protected:
		struct rt_mesh
		{
			che * mesh;
			bool pointcloud;

			che * operator -> () const
			{
				return mesh;
			}
		};

		std::map<index_t, rt_mesh> geomID_mesh;

		size_t n_samples = 0;

		static std::default_random_engine gen;
		static std::uniform_real_distribution<float> randf;

	public:
		raytracing() = default;
		virtual ~raytracing() = default;

		virtual void render(vec4 * img, const render_params & params, const bool & flat);

		virtual float * raycaster(	const uvec2 & windows_size,
									const mat4 & inv_proj_view,
									const vertex & cam_pos,
									const index_t & samples = 4
									);

		vertex ray_view_dir(	const index_t & x, const index_t & y,
								const uvec2 & windows_size,
								const mat4 & inv_proj_view,
								const vertex & cam_pos
								);

		virtual hit intersect(	const vertex &,	// org
								const vertex &	//dir
								)	{ return hit(); }

		virtual index_t closest_vertex(	const vertex &,	// org,
										const vertex &	// dir
										) { return NIL; };

	protected:
		virtual vec4 intersect_li(	const vertex &,	// org,
										const vertex &,	// dir,
										const vertex &,	// light,
										const bool &		// flat
										) { return vec4(0); };

		virtual float intersect_depth(	const vertex &,	// org,
										const vertex &	// dir
										) { return 0; };
};


} // namespace gproshan

#endif // RAYTRACING_H

