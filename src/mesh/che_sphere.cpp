#include "mesh/che_sphere.h"

#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_sphere::che_sphere(const real_t & r, const size_t & n)
{
	filename = "sphere";

	std::vector<vertex> vertices;
	std::vector<index_t> faces;

	const real_t delta = M_PI / n;

	for(real_t phi = 0; phi < 2 * M_PI - 0.5 * delta; phi += delta)
	for(real_t theta = delta; theta < M_PI - 0.5 * delta; theta += delta)
		vertices.push_back({r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)});

	vertices.push_back({0, 0, r});
	vertices.push_back({0, 0, -r});

	size_t v, cols = n - 1;

	for(index_t i = 0; i < 2 * n - 1; ++i)
	{
		for(index_t j = 0; j < cols - 1; ++j)
		{
			v = i * cols + j;

			faces.push_back(v);
			faces.push_back(v + 1);
			faces.push_back(v + cols);

			faces.push_back(v + cols);
			faces.push_back(v + 1);
			faces.push_back(v + cols + 1);
		}

		v = i * cols;
		faces.push_back(vertices.size() - 2);
		faces.push_back(v);
		faces.push_back(v + cols);

		v = (i + 1) * cols - 1;
		faces.push_back(vertices.size() - 1);
		faces.push_back(v + cols);
		faces.push_back(v);
	}

	for(index_t j = 0; j < cols - 1; ++j)
	{
		v = (2 * n - 1) * cols + j;

		faces.push_back(v + 1);
		faces.push_back(j);
		faces.push_back(v);

		faces.push_back(j + 1);
		faces.push_back(j);
		faces.push_back(v + 1);
	}

	v = (2 * n - 1) * cols;
	faces.push_back(vertices.size() - 2);
	faces.push_back(v);
	faces.push_back(0);

	v = (2 * n) * cols - 1;
	faces.push_back(vertices.size() - 1);
	faces.push_back(cols - 1);
	faces.push_back(v);

	init(vertices.data(), vertices.size(), faces.data(), faces.size() / 3);
}


} // namespace gproshan

