#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/render_params.h>

#include <vector>
#include <map>
#include <random>

#include <glm/glm.hpp>


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
	public:
		bool restart = false;

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

		virtual void render(glm::vec4 * img, const render_params & params, const bool & flat);

		virtual float * raycaster(	const glm::uvec2 & windows_size,
									const glm::mat4 & proj_view_mat,
									const glm::vec3 & cam_pos,
									const index_t & samples = 4
									);

		glm::vec3 ray_view_dir(	const index_t & x, const index_t & y,
								const glm::vec2 & windows_size,
								const glm::mat4 & inv_proj_view,
								const glm::vec3 & cam_pos
								);

		virtual hit intersect(	const glm::vec3 &,	// org
								const glm::vec3 &	//dir
								)	{ return hit(); }

		virtual index_t closest_vertex(	const glm::vec3 &,	// org,
										const glm::vec3 &	// dir
										) { return NIL; };

	protected:
		virtual glm::vec4 intersect_li(	const glm::vec3 &,	// org,
										const glm::vec3 &,	// dir,
										const glm::vec3 &,	// light,
										const bool &		// flat
										) { return glm::vec4(0); };

		virtual float intersect_depth(	const glm::vec3 &,	// org,
										const glm::vec3 &	// dir
										) { return 0; };
};


} // namespace gproshan

#endif // RAYTRACING_H

