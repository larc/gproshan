#include "raytracing/raytracing.h"

#include <cstring>
#include <random>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


raytracing::raytracing()
{
	width = height = n_samples = 0;
	img = nullptr;
}

raytracing::~raytracing()
{
	if(img) delete [] img;
}

bool raytracing::rt_restart(const size_t & w, const size_t & h)
{
	if(width * height < w * h)
	{
		width = w;
		height = h;

		delete [] img;
		img = new glm::vec4[width * height];

		return true;
	}

	if(width != w || height != h)
	{
		width = w;
		height = h;

		return true;
	}

	return false;
}

void raytracing::render(	const glm::uvec2 & windows_size,
								const glm::mat4 & view_mat,
								const glm::mat4 & proj_mat,
								const std::vector<glm::vec3> & light,
								const bool & flat,
								const bool & restart )
{
	if(rt_restart(windows_size.x, windows_size.y) || restart)
		n_samples = 0;

	std::default_random_engine gen;
	std::uniform_real_distribution<float> randf(0, 1);

	glm::vec3 cam_pos = glm::vec3(glm::inverse(view_mat) * glm::vec4(0, 0, 0, 1));
	glm::mat4 inv_proj_view = glm::inverse(proj_mat * view_mat);

	glm::vec4 li;

	if(!n_samples)
	{
		#pragma omp parallel for
		for(index_t i = 0; i < width; ++i)
		for(index_t j = 0; j < height; ++j)
			img[j * width + i] = glm::vec4(0);
	}

	#pragma omp parallel for private(li)
	for(index_t i = 0; i < width; ++i)
	for(index_t j = 0; j < height; ++j)
	{
		//row major
		glm::vec4 & color = img[j * width + i];

		glm::vec2 screen = glm::vec2(	(float(i) + randf(gen)) / width,
										(float(j) + randf(gen)) / height
										);

		glm::vec4 view = glm::vec4(screen.x * 2 - 1, screen.y * 2 - 1, 1, 1);
		glm::vec4 q = inv_proj_view * view;
		glm::vec3 p = glm::vec3(q * (1.f / q.w));

		li = glm::vec4(0);
		for(auto & l: light)
			li += intersect_li(cam_pos, glm::normalize(p - cam_pos), l, flat);

		color = (color * float(n_samples) + li / float(light.size())) / float(n_samples + 1);
	}

	++n_samples;
}

float * raytracing::raycaster(	const glm::uvec2 & windows_size,
								const glm::mat4 & view_mat,
								const glm::mat4 & proj_mat,
								const index_t & samples	)
{
	float * frame = new float[windows_size.x * windows_size.y];

	std::default_random_engine gen;
	std::uniform_real_distribution<float> randf(0.f, 1.f);

	glm::vec3 cam_pos = glm::vec3(glm::inverse(view_mat) * glm::vec4(0.f, 0.f, 0.f, 1.f));
	glm::mat4 inv_proj_view = glm::inverse(proj_mat * view_mat);

	#pragma omp parallel for
	for(index_t i = 0; i < windows_size.x; ++i)
	for(index_t j = 0; j < windows_size.y; ++j)
	{
		//row major
		float & color = frame[(windows_size.y - j - 1) * windows_size.x + i] = 0;

		for(index_t s = 0; s < samples; ++s)
		{
			glm::vec2 screen = glm::vec2(	(float(i) + randf(gen)) / windows_size.x,
											(float(j) + randf(gen)) / windows_size.y
											);

			glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
			glm::vec4 q = inv_proj_view * view;
			glm::vec3 p = glm::vec3(q * (1.f / q.w));

			color += intersect_depth(cam_pos, glm::normalize(p - cam_pos));
		}

		color /= samples;
	}

	return frame;
}


} // namespace gproshan

