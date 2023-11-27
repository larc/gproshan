#include <gproshan/mesh/che_sphere.h>

#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>


// geometry processing and shape analysis framework
namespace gproshan {


che_sphere::che_sphere(const real_t & r, const size_t & n)
{
	filename = "sphere";

	std::vector<vertex> vertices;
	std::vector<index_t> trigs;

	const real_t delta = M_PI / n;

	for(real_t phi = 0; phi < 2 * M_PI - 0.5 * delta; phi += delta)
	for(real_t theta = delta; theta < M_PI - 0.5 * delta; theta += delta)
		vertices.push_back({r * std::sin(theta) * std::cos(phi), r * std::sin(theta) * std::sin(phi), r * std::cos(theta)});

	vertices.push_back({0, 0, r});
	vertices.push_back({0, 0, -r});

	size_t v, cols = n - 1;

	for(index_t i = 0; i < 2 * n - 1; ++i)
	{
		for(index_t j = 0; j < cols - 1; ++j)
		{
			v = i * cols + j;

			trigs.push_back(v);
			trigs.push_back(v + 1);
			trigs.push_back(v + cols);

			trigs.push_back(v + cols);
			trigs.push_back(v + 1);
			trigs.push_back(v + cols + 1);
		}

		v = i * cols;
		trigs.push_back(size(vertices) - 2);
		trigs.push_back(v);
		trigs.push_back(v + cols);

		v = (i + 1) * cols - 1;
		trigs.push_back(size(vertices) - 1);
		trigs.push_back(v + cols);
		trigs.push_back(v);
	}

	for(index_t j = 0; j < cols - 1; ++j)
	{
		v = (2 * n - 1) * cols + j;

		trigs.push_back(v + 1);
		trigs.push_back(j);
		trigs.push_back(v);

		trigs.push_back(j + 1);
		trigs.push_back(j);
		trigs.push_back(v + 1);
	}

	v = (2 * n - 1) * cols;
	trigs.push_back(size(vertices) - 2);
	trigs.push_back(v);
	trigs.push_back(0);

	v = (2 * n) * cols - 1;
	trigs.push_back(size(vertices) - 1);
	trigs.push_back(cols - 1);
	trigs.push_back(v);

	init(vertices.data(), size(vertices), trigs.data(), size(trigs) / 3);
}


} // namespace gproshan

