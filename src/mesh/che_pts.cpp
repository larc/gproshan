#include "mesh/che_pts.h"

#include <cstring>
#include <cstdio>
#include <cassert>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_pts::che_pts(const string & file)
{
	init(file);
}

void che_pts::read_file(const string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	float x, y, z;
	int intensity;
	unsigned char r, g, b;
	size_t n;

	fscanf(fp, "%lu", &n);

	alloc(n, 0);
	
	for(index_t v = 0; v < n_vertices; ++v)
	{
		n = fscanf(fp, "%f %f %f %d %hhu %hhu %hhu", &x, &y, &z, &intensity, &r, &g, &b);

		GT[v] = {x, y, z};
		if(n == 7)
		{
			VC[v] = {r, g, b};
			VHC[v] = float(intensity + 2048) / 4095;
		}
	}

	fclose(fp);
}

void che_pts::write_file(const che * mesh, const string & file)
{
	FILE * fp = fopen((file + ".pts").c_str(), "w");
	assert(fp);

	fprintf(fp, "%lu\n", mesh->n_vertices);
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		const vertex & v = mesh->gt(i);
		const rgb_t & c = mesh->rgb(i);
		fprintf(fp, "%f %f %f", (float) v.x, (float) v.y, (float) v.z);
		fprintf(fp, " %d ", int(mesh->heatmap(i) * 4095) - 2048);
		fprintf(fp, "%hhu %hhu %hhu\n", c.r, c.g, c.b);
	}
	fclose(fp);
}


} // namespace gproshan

