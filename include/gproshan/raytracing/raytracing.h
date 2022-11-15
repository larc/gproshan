#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/render_params.h>
#include <gproshan/raytracing/rt_utils.h>

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class raytracing
{
	protected:
		size_t n_samples = 0;

	public:
		raytracing() = default;
		virtual ~raytracing() = default;

		virtual void render(vec4 * img, const render_params & params, const bool & flat);

		virtual float * raycaster(	const ivec2 & windows_size,
									const mat4 & inv_proj_view,
									const vertex & cam_pos,
									const index_t & samples = 4
									);

		virtual eval_hit intersect(	const vertex &,	// org
									const vertex &	//dir
									) { return {}; }

		virtual index_t closest_vertex(	const vertex &,	// org,
										const vertex &	// dir
										) { return NIL; };

	protected:
		virtual vec3 closesthit_radiance(	const vertex &, // org
											const vertex &, // dir
											const vertex *, // lights
											const int &, // n_lights
											const bool & // flat
											) { return {}; };

		virtual float intersect_depth(	const vertex &,	// org,
										const vertex &	// dir
										) { return 0; };
};


} // namespace gproshan

#endif // RAYTRACING_H

