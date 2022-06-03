#include <gproshan/raytracing/raytracing.h>

#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


std::default_random_engine raytracing::gen;
std::uniform_real_distribution<float> raytracing::randf;

void raytracing::render(glm::vec4 * img, const render_params & params, const bool & flat)
{
	if(restart) n_samples = 0;

	glm::mat4 inv_proj_view = glm::inverse(params.proj_view_mat);

	int window_width = params.window_width;
	int window_height = params.window_height;
	if(params.viewport_is_window)
	{
		window_width = params.viewport_width;
		window_height = params.viewport_height;
	}

	glm::vec4 li;

	#pragma omp parallel for private(li)
	for(int i = 0; i < params.viewport_width; ++i)
	for(int j = 0; j < params.viewport_height; ++j)
	{
		//row major
		glm::vec4 & color = img[j * params.viewport_width + i];
		const glm::vec3 & dir = ray_view_dir(	i + params.viewport_x,
												j + params.viewport_y,
												{window_width, window_height},
												inv_proj_view,
												params.cam_pos
												);

		li = glm::vec4(0);
		for(auto & l: params.lights)
			li += intersect_li(params.cam_pos, dir, l, flat);

		color = (color * float(n_samples) + li / float(params.lights.size())) / float(n_samples + 1);
	}

	++n_samples;
}

float * raytracing::raycaster(	const glm::uvec2 & windows_size,
								const glm::mat4 & proj_view_mat,
								const glm::vec3 & cam_pos,
								const index_t & samples	)
{
	float * frame = new float[windows_size.x * windows_size.y];

	glm::mat4 inv_proj_view = glm::inverse(proj_view_mat);

	#pragma omp parallel for
	for(index_t i = 0; i < windows_size.x; ++i)
	for(index_t j = 0; j < windows_size.y; ++j)
	{
		//row major
		float & color = frame[(windows_size.y - j - 1) * windows_size.x + i] = 0;
		glm::vec3 dir = ray_view_dir(i, j, windows_size, inv_proj_view, cam_pos);

		for(index_t s = 0; s < samples; ++s)
			color += intersect_depth(cam_pos, dir);

		color /= samples;
	}

	return frame;
}

glm::vec3 raytracing::ray_view_dir(const index_t & x, const index_t & y, const glm::vec2 & windows_size, const glm::mat4 & inv_proj_view, const glm::vec3 & cam_pos)
{
	glm::vec2 screen = glm::vec2((float(x) + randf(gen)) / windows_size.x, (float(y) + randf(gen)) / windows_size.y);
	glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
	glm::vec4 q = inv_proj_view * view;
	glm::vec3 p = glm::vec3(q * (1.f / q.w));

	return glm::normalize(p - cam_pos);
}


} // namespace gproshan
