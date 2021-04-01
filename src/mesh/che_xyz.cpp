#include "mesh/che_xyz.h"

#include <cstring>
#include <cstdio>
#include <cassert>
#include <fstream>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_xyz::che_xyz(const string & file)
{
	init(file);
}

void che_xyz::read_file(const string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	char line[256], c;
	float x, y, z, r, g, b;
	size_t n;

	vector<vertex> vertices;
	vector<vertex> vertices_color;

	while(fgets(line, sizeof(line), fp))
	{
		sscanf(line, "%c", &c);
		if(c == '#') continue;

		n = sscanf(line, "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
		vertices.push_back({x, y, z});
		vertices_color.push_back(n == 6 ? vertex{r, g, b} / 255 : vcolor);
	}

	fclose(fp);

	alloc(vertices.size(), 0);
	memcpy(GT, vertices.data(), n_vertices * sizeof(vertex));
	memcpy(VC, vertices_color.data(), n_vertices * sizeof(vertex));
}

void che_xyz::write_file(const che * mesh, const string & file)
{
}


} // namespace gproshan

