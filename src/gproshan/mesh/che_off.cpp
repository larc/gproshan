#include <gproshan/mesh/che_off.h>

#include <cstring>
#include <cstdio>
#include <cassert>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_off::che_off(const string & file)
{
	init(file);
}

void che_off::read_file(const string & file)
{
	char soff[32];
	size_t nv, nf, n;

	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	fgets(soff, sizeof(soff), fp);
	fscanf(fp, "%lu %lu %lu", &nv, &nf, &n);

	alloc(nv, nf);

	float x, y, z;
	unsigned char r, g, b, a;
	for(index_t v = 0; v < n_vertices; ++v)
	{
		fscanf(fp, "%f %f %f", &x, &y, &z);
		GT[v] = { x, y, z };

		if(soff[0] == 'C' || soff[1] == 'C')
		{
			fscanf(fp, "%hhu %hhu %hhu %hhu", &r, &g, &b, &a);
			VC[v] = { r, g, b };
		}

		if(soff[0] == 'N')
		{
			fscanf(fp, "%f %f %f", &x, &y, &z);
			VN[v] = { x, y, z };
		}
	}

	vector<index_t> faces;
	faces.reserve(che::mtrig * n_faces);

	index_t P[32];
	while(nf--)
	{
		fscanf(fp, "%lu", &n);
		for(index_t i = 0; i < n; ++i)
			fscanf(fp, "%u", P + i);

		for(const index_t & v: trig_convex_polygon(P, n))
			faces.push_back(v);
	}

	fclose(fp);


	if(faces.size() != che::mtrig * n_faces)
	{
		vertex * tGT = GT; GT = nullptr;
		vertex * tVN = VN; VN = nullptr;
		rgb_t * tVC = VC; VC = nullptr;

		free();
		alloc(nv, faces.size() / che::mtrig);

		GT = tGT;
		VN = tVN;
		VC = tVC;
	}

	memcpy(VT, faces.data(), faces.size() * sizeof(index_t));
}

void che_off::write_file(const che * mesh, const string & file, const che_off::type & off, const bool & pointcloud)
{
	static const char * str_off[] = {"OFF", "NOFF", "COFF", "NCOFF"};

	FILE * fp = fopen((file + ".off").c_str(), "w");
	assert(fp);

	fprintf(fp, "%s\n", str_off[off]);
	fprintf(fp, "%lu %lu 0\n", mesh->n_vertices, pointcloud ? 0 : mesh->n_faces);

	for(size_t i = 0; i < mesh->n_vertices; ++i)
	{
		const vertex & v = mesh->point(i);
		fprintf(fp, "%f %f %f", (float) v.x(), (float) v.y(), (float) v.z());

		if(off == COFF || off == NCOFF)
		{
			const rgb_t & c = mesh->rgb(i);
			fprintf(fp, " %hhu %hhu %hhu 1", c.r, c.g, c.b);
		}

		if(off == NOFF || off == NCOFF)
		{
			const vertex & n = mesh->normal(i);
			fprintf(fp, " %f %f %f", (float) n.x(), (float) n.y(), (float) n.z());
		}

		fprintf(fp, "\n");
	}

	if(!pointcloud)
	{
		for(index_t he = 0; he < mesh->n_half_edges; )
		{
			fprintf(fp, "%u", che::mtrig);
			for(index_t i = 0; i < che::mtrig; ++i)
				fprintf(fp, " %u", mesh->halfedge(he++));
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}


} // namespace gproshan

