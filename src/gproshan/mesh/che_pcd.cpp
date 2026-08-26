#include <gproshan/mesh/che_pcd.h>

#include <cstring>
#include <cstdio>
#include <cassert>
#include <unordered_map>


// geometry processing and shape analysis framework
namespace gproshan {


che_pcd::che_pcd(const std::string & file)
{
	init(file);
}

void che_pcd::read_file(const std::string & file)
{
	std::unordered_map<std::string, size_t> bytes = {
														{"char", 1},
														{"uchar", 1},
														{"short", 2},
														{"ushort", 2},
														{"int", 4},
														{"uint", 4},
														{"float", 4},
														{"float32", 4},
														{"float64", 8},
														{"double", 8}
													};

	FILE * fp = fopen(file.c_str(), "rb");
	assert(fp);

	size_t n_points = 0;
	char line[512], str[32], format[32];

	while(fgets(line, sizeof(line), fp) && str[0] != 'D')	// DATA is the end header
	{
		sscanf(line, "%s", str);

		switch(str[0])
		{
			case 'V':	// VERSION / VIEWPOINT
				break;
			case 'F':	// FIELDS
				break;
			case 'S':	// SIZE
				break;
			case 'T':	// TYPE
				break;
			case 'C':	// COUNT
				break;
			case 'W':	// WIDTH
				break;
			case 'H':	// HEIGHT
				break;
			case 'P':	// POINTS
				sscanf(line, "%*s %lu", &n_points);
				break;
			case 'D':	// DATA
				sscanf(line, "%*s %s", format);
				break;
		}
	}

	alloc(n_points, 0);


	if(format[0] == 'a')	// ascii
	{
		float x, y, z;
		for(index_t v = 0; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);
			sscanf(line, "%f %f %f", &x, &y, &z);
			GT[v] = {x, y, z};
		}
	}

	fclose(fp);
}

void che_pcd::write_file(const che * mesh, const std::string & file)
{
	FILE * fp = fopen((file + ".pcd").c_str(), "wb");
	assert(fp);

	fprintf(fp, "# .PCD v.7 - Point Cloud Data file format\n");
	fprintf(fp, "# written by gproshan\n");
	fprintf(fp, "VERSION .7\n");
	fprintf(fp, "FIELDS x y z\n");
	fprintf(fp, "SIZE 4 4 4\n");
	fprintf(fp, "TYPE F F F\n");
	fprintf(fp, "COUNT 1 1 1\n");
	fprintf(fp, "WIDTH %lu\n", mesh->n_vertices);
	fprintf(fp, "HEIGHT 1\n");
	fprintf(fp, "VIEWPOINT 0 0 0 1 0 0 0\n");
	fprintf(fp, "POINTS %lu\n", mesh->n_vertices);
	fprintf(fp, "DATA ascii\n");

	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		const vertex & p = mesh->point(v);
		fprintf(fp, "%f %f %f\n", p.x(), p.y(), p.z());
	}

	fclose(fp);
}


} // namespace gproshan

