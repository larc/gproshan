#include <gproshan/raytracing/raytracing.h>

#include <cstring>


// geometry processing and shape analysis framework
namespace gproshan::rt {


void raytracing::render(vec4 * img, const render_params & params, const bool & flat)
{
	uvec2 window_size = params.window_size;
	if(params.viewport_is_window)
		window_size = params.viewport_size;

	vec3 color_acc, position, ray_dir, attenuation;

	#pragma omp parallel for private(color_acc, position, ray_dir, attenuation)
	for(unsigned int i = 0; i < params.viewport_size.x(); ++i)
	for(unsigned int j = 0; j < params.viewport_size.y(); ++j)
	{
		const uvec2 & pos = params.viewport_pos + uvec2{i, j};

		random<real_t> rnd(pos.x() + window_size.x() * pos.y(), params.n_frames);

		int depth	= params.depth;
		int samples	= params.n_samples;

		do
		{
			color_acc	= 0;
			attenuation = 1;
			position	= params.cam_pos;
			ray_dir		= ray_view_dir(pos, window_size, params.inv_proj_view, params.cam_pos, rnd);

			depth = params.depth;
			do
			{
				const vec3 & color = closesthit_radiance(position, ray_dir, params.ambient, params.lights, params.n_lights, params.cam_pos, flat);
				color_acc += color * attenuation;
			}
			while(--depth);
		}
		while(--samples);

		color_acc /= params.n_samples;

		vec4 & pixel_color = img[i + j * params.viewport_size.x()];
		pixel_color = (pixel_color * params.n_frames + (color_acc, 1)) / (params.n_frames + 1);
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

