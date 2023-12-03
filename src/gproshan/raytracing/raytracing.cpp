#include <gproshan/raytracing/raytracing.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan::rt {


void raytracing::render(vec4 * img, const render_params & params, const bool & flat)
{
	uvec2 window_size = params.window_size;
	if(params.viewport_is_window)
		window_size = params.viewport_size;

	#pragma omp parallel for
	for(unsigned int i = 0; i < params.viewport_size.x(); ++i)
	for(unsigned int j = 0; j < params.viewport_size.y(); ++j)
	{
		const uvec2 & pos = params.viewport_size + uvec2{i, j};

		random<real_t> rnd(pos.x(), pos.y());
	
		//row major
		vec4 & color = img[j * params.viewport_size.x() + i];
		const vertex & dir = ray_view_dir(pos, window_size, params.inv_proj_view, params.cam_pos, rnd);

		const vec3 & li = closesthit_radiance(params.cam_pos, dir, params.ambient, params.lights, params.n_lights, params.cam_pos, flat);

		color = (color * params.n_frames + (li, 1)) / (params.n_frames + 1);
	}
}

std::vector<float> raytracing::raycaster(	const uvec2 & windows_size,
											const mat4 & inv_proj_view,
											const vertex & cam_pos,
											const index_t & samples
											) const
{
	std::vector<float> frame(windows_size.x() * windows_size.y());

	#pragma omp parallel for
	for(unsigned int i = 0; i < windows_size.x(); ++i)
	for(unsigned int j = 0; j < windows_size.y(); ++j)
	{
		random<real_t> rnd(i, j);
		
		//row major
		float & color = frame[(windows_size.y() - j - 1) * windows_size.x() + i] = 0;
		vertex dir = ray_view_dir({i, j}, windows_size, inv_proj_view, cam_pos, rnd);

		for(index_t s = 0; s < samples; ++s)
			color += intersect_depth(cam_pos, dir);

		color /= samples;
	}

	return frame;
}


} // namespace gproshan

