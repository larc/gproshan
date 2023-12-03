#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/render_params.h>
#include <gproshan/raytracing/utils.h>

#include <vector>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class raytracing
{
	public:
		raytracing() = default;
		virtual ~raytracing() = default;

		virtual void render(vec4 * img, const render_params & params, const bool & flat);

		virtual std::vector<float> raycaster(	const uvec2 & windows_size,
												const mat4 & inv_proj_view,
												const vertex & cam_pos,
												const index_t & samples = 4
												) const;

		virtual eval_hit intersect(	const vertex &,	// org
									const vertex &,	//dir
									[[maybe_unused]] const bool & flat = true
									) const { return {}; }

		virtual index_t closest_vertex(	const vertex &,	// org,
										const vertex &	// dir
										) const { return NIL; }

	protected:
		virtual vec3 closesthit_radiance(	const vertex &, // org
											const vertex &, // dir
											const light &, // ambient light
											const light *, // lights
											const int &, // n_lights
											const vertex &, // cam_pos
											const bool & // flat
											) const { return {}; };

		virtual float intersect_depth(	const vertex &,	// org,
										const vertex &	// dir
										) const { return 0; };
};


} // namespace gproshan

#endif // RAYTRACING_H

