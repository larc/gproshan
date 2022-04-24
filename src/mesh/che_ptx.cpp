#include "mesh/che_ptx.h"

#include <cstring>
#include <cassert>
#include <cstdio>
#include <thread>

#include <CImg.h>

using namespace cimg_library;
using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_ptx::che_ptx(const string & file)
{
	init(file);
}

void che_ptx::read_file(const string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	size_t n_rows, n_cols;
	float T[12], R[12], tr[4];

	fscanf(fp, "%lu %lu", &n_rows, &n_cols);

	for(index_t i = 0; i < 12; ++i)
		fscanf(fp, "%f", T + i);

	for(index_t i = 0; i < 12; ++i)
		fscanf(fp, "%f", R + i);

	for(index_t i = 0; i < 4; ++i)
		fscanf(fp, "%f", tr + i);


	alloc(n_rows * n_cols, 2 * (n_rows - 1) * (n_cols - 1));

	float x, y, z, intensity;
	unsigned char r, g, b;
	char line[128];

	bool rgb = false;

	// vertex 0: x y z a or x y z a r g b
	fgets(line, sizeof(line), fp);
	fgets(line, sizeof(line), fp);
	rgb = sscanf(line, "%f %f %f %f %hhu %hhu %hhu", &x, &y, &z, &intensity, &r, &g, &b) == 7;

	if(rgb)
	{
		GT[0] = { x, y, z };
		VC[0] = { r, g, b };
		VHC[0] = intensity;

		for(index_t v = 1; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);
			sscanf(line, "%f %f %f %f %hhu %hhu %hhu", &x, &y, &z, &intensity, &r, &g, &b);
			GT[v] = { x, y, z };
			VC[v] = { r, g, b };
			VHC[v] = intensity;
		}

		CImg<unsigned char> img((unsigned char *) VC, 3, n_cols, n_rows);
		img.permute_axes("zycx");
		img.save((file + ".jpg").c_str());

		thread([](CImg<real_t> img) { img.mirror("y").display(); }, img).detach();
	}
	else
	{
		GT[0] = { x, y, z };
		VHC[0] = intensity;

		for(index_t v = 1; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);
			sscanf(line, "%f %f %f %f", &x, &y, &z, &intensity);
			GT[v] = { x, y, z };
			VHC[v] = intensity;
		}
	}

	fclose(fp);


	index_t he = 0;
	auto add_trig = [&](const index_t & i, const index_t & j, const index_t & k)
	{
		if(GT[i].is_zero() || GT[j].is_zero() || GT[k].is_zero())
			return;

		VT[he++] = i;
		VT[he++] = j;
		VT[he++] = k;

		if(pdetriq(trig(he - 1)) < 0.1)
			he -= 3;
	};

	for(index_t r = 0; r < n_rows - 1; ++r)
	for(index_t c = 0; c < n_cols - 1; ++c)
	{
		add_trig((c    ) + (r    ) * n_cols,
				 (c    ) + (r + 1) * n_cols,
				 (c + 1) + (r    ) * n_cols);

		add_trig((c + 1) + (r + 1) * n_cols,
				 (c + 1) + (r    ) * n_cols,
				 (c    ) + (r + 1) * n_cols);
	}

	rw(n_half_edges)	= he;
	rw(n_faces)			= he / che::mtrig;
}


} // namespace gproshan

