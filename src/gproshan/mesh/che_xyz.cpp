#include <gproshan/mesh/che_xyz.h>

#include <cstring>
#include <cstdio>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


che_xyz::che_xyz(const std::string & file)
{
	init(file);
}

void che_xyz::read_file(const std::string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	char line[256], c;
	float x, y, z;
	unsigned char r, g, b;
	size_t n;

	std::vector<vertex> vertices;
	std::vector<rgb_t> vertices_color;

	while(fgets(line, sizeof(line), fp))
	{
		sscanf(line, "%c", &c);
		if(c == '#') continue;

		n = sscanf(line, "%f %f %f %hhu %hhu %hhu", &x, &y, &z, &r, &g, &b);
		vertices.push_back({x, y, z});
		vertices_color.push_back(n == 6 ? rgb_t{r, g, b} : rgb_t());
	}

	fclose(fp);

	alloc(size(vertices), 0);
	memcpy(GT, vertices.data(), n_vertices * sizeof(vertex));
	memcpy(VC, vertices_color.data(), n_vertices * sizeof(rgb_t));
}

void che_xyz::write_file(const che * mesh, const std::string & file, const bool & color)
{
	FILE * fp = fopen((file + ".xyz").c_str(), "w");
	assert(fp);

	fprintf(fp, "# XYZ generated by gproshan\n");
	fprintf(fp, "# vertices %lu\n", mesh->n_vertices);

	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		const vertex & v = mesh->point(i);
		fprintf(fp, "%f %f %f", (float) v.x(), (float) v.y(), (float) v.z());
		if(color)
		{
			const rgb_t & c = mesh->rgb(i);
			fprintf(fp, " %hhu %hhu %hhu", c.r, c.g, c.b);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


} // namespace gproshan

