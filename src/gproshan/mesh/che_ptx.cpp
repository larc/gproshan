#include <gproshan/mesh/che_ptx.h>

#include <cstring>
#include <cassert>
#include <cstdio>
#include <thread>

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


che_ptx::che_ptx(const std::string & file)
{
	init(file);
}

void che_ptx::read_file(const std::string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	size_t n_rows, n_cols;
	vertex p;	// scanner position
	mat3 A;		// scanner axis
	mat4 T;		// transformation matrix

	fscanf(fp, "%lu %lu", &n_rows, &n_cols);
	fscanf(fp, "%f %f %f", &p.x(), &p.y(), &p.z());

	for(index_t i = 0; i < 3; ++i)
	for(index_t j = 0; j < 3; ++j)
		fscanf(fp, "%f", &A(i, j));

	for(index_t i = 0; i < 4; ++i)
	for(index_t j = 0; j < 4; ++j)
		fscanf(fp, "%f", &T(i, j));

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
		GT[0] = {x, y, z};
		VC[0] = {r, g, b};
		VHC[0] = intensity;

		for(index_t v = 1; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);
			sscanf(line, "%f %f %f %f %hhu %hhu %hhu", &x, &y, &z, &intensity, &r, &g, &b);
			GT[v] = {x, y, z};
			VC[v] = {r, g, b};
			VHC[v] = intensity;
		}

	#ifndef NDEBUG
		CImg<unsigned char> img((unsigned char *) VC, 3, n_cols, n_rows);
		img.permute_axes("zycx");
		img.save((file + ".jpg").c_str());

		std::thread([](CImg<unsigned char> img) { img.mirror("y").display(); }, img).detach();
	#endif // NDEBUG
	}
	else
	{
		GT[0] = {x, y, z};
		VHC[0] = intensity;

		for(index_t v = 1; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);
			sscanf(line, "%f %f %f %f", &x, &y, &z, &intensity);
			GT[v] = {x, y, z};
			VHC[v] = intensity;
		}
	}

	fclose(fp);


	index_t he = 0;
	auto add_trig = [&](const index_t i, const index_t j, const index_t k)
	{
		VT[he++] = i;
		VT[he++] = j;
		VT[he++] = k;

		if(pdetriq(he_trig(he - 1)) < 0.1)
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
	rw(n_trigs)			= he / che::mtrig;
}

void che_ptx::write_file(const che * mesh, const std::string & file, const size_t & n_rows, const size_t & n_cols)
{
	FILE * fp = fopen((file + ".ptx").c_str(), "wb");
	assert(fp);

	fprintf(fp, "%lu\n", n_rows);
	fprintf(fp, "%lu\n", n_cols);
	fprintf(fp, "0 0 0\n");
	fprintf(fp, "1 0 0\n");
	fprintf(fp, "0 1 0\n");
	fprintf(fp, "0 0 1\n");
	fprintf(fp, "1 0 0 0\n");
	fprintf(fp, "0 1 0 0\n");
	fprintf(fp, "0 0 1 0\n");
	fprintf(fp, "0 0 0 1\n");

	for(size_t i = 0; i < mesh->n_vertices; ++i)
	{
		const vertex & v = mesh->point(i);
		const rgb_t & c = mesh->rgb(i);
		const real_t & h = mesh->heatmap(i);

		fprintf(fp, "%f %f %f %f %hhu %hhu %hhu\n", (float) v.x(), (float) v.y(), (float) v.z(), (float) h, c.r, c.g, c.b );
	}

	fclose(fp);
}


} // namespace gproshan

