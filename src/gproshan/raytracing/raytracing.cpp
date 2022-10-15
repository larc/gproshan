#include <gproshan/raytracing/raytracing.h>

#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


void raytracing::render(vec4 * img, const render_params & params, const bool & flat)
{
	if(params.restart) n_samples = 0;

	int window_width = params.window_width;
	int window_height = params.window_height;
	if(params.viewport_is_window)
	{
		window_width = params.viewport_width;
		window_height = params.viewport_height;
	}

	vec4 li;

	#pragma omp parallel for private(li)
	for(int i = 0; i < params.viewport_width; ++i)
	for(int j = 0; j < params.viewport_height; ++j)
	{
		//row major
		vec4 & color = img[j * params.viewport_width + i];
		const vertex & dir = ray_view_dir(	{i + params.viewport_x, j + params.viewport_y},
											{window_width, window_height},
											params.inv_proj_view,
											params.cam_pos
											);

		li = vec4(0);
		for(int l = 0; l < params.n_lights; ++l)
			li += intersect_li(params.cam_pos, dir, params.lights[l], flat);

		color = (color * n_samples + li / params.n_lights) / (n_samples + 1);
	}

	++n_samples;
}

float * raytracing::raycaster(	const ivec2 & windows_size,
								const mat4 & inv_proj_view,
								const vertex & cam_pos,
								const index_t & samples	)
{
	float * frame = new float[windows_size.x() * windows_size.y()];

	#pragma omp parallel for
	for(int i = 0; i < windows_size.x(); ++i)
	for(int j = 0; j < windows_size.y(); ++j)
	{
		//row major
		float & color = frame[(windows_size.y() - j - 1) * windows_size.x() + i] = 0;
		vertex dir = ray_view_dir({i, j}, windows_size, inv_proj_view, cam_pos);

		for(index_t s = 0; s < samples; ++s)
			color += intersect_depth(cam_pos, dir);

		color /= samples;
	}

	return frame;
}


} // namespace gproshan

