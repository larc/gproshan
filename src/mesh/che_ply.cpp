#include "mesh/che_ply.h"

#include <cstring>
#include <cstdio>
#include <cassert>
#include <map>


using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_ply::che_ply(const string & file)
{
	init(file);
}

void che_ply::read_file(const string & file)
{
	map<string, size_t> bytes = {
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

	size_t nv = 0, nf = 0;
	size_t nb, xyz = 0, rgb = 0, vbytes = 0;
	index_t ixyz = 0, irgb = 0;
	size_t fn = 0, fbytes = 0;
	index_t P[32];

	char line[512], type[32], str[32], format[32], element[32];

	while(fgets(line, sizeof(line), fp) && line[1] != 'n')	// end_header
	{
		sscanf(line, "%s", str);

		if(str[0] == 'f')	// format
			sscanf(line, "%*s %s", format);

		if(str[0] == 'e')	// element
		{
			sscanf(line, "%*s %s", element);
			if(element[0] == 'v')	// vertex
				sscanf(line, "%*s %*s %lu", &nv);
			if(element[0] == 'f')	// face
				sscanf(line, "%*s %*s %lu", &nf);
		}

		if(str[0] == 'p' && element[0] == 'v')	// property vertex
		{
			sscanf(line, "%*s %s %s", type, str);
			nb = bytes[type];

			if(str[0] == 'x')
			{
				xyz = nb;
				ixyz = vbytes;
			}
			if(str[0] == 'r')
			{
				rgb = nb;
				irgb = vbytes;
			}

			vbytes += nb;
		}

		if(str[0] == 'p' && element[0] == 'f')	// property face
		{
			sscanf(line, "%*s %s", str);
			if(str[0] == 'l')	// list
			{
				sscanf(line, "%*s %*s %s", type);
				fn = bytes[type];
				sscanf(line, "%*s %*s %*s %s", type);
				fbytes = bytes[type];
			}
		}
	}

	alloc(nv, nf);

	vector<index_t> faces;
	faces.reserve(che::mtrig * n_faces);

	if(format[0] == 'a')	// ascii
	{
		float x, y, z;
		for(index_t v = 0; v < n_vertices; ++v)
		{
			fscanf(fp, "%f %f %f", &x, &y, &z);
			GT[v] = {x, y, z};
		}

		while(nf--)
		{
			fscanf(fp, "%lu", &nv);
			for(index_t i = 0; i < nv; ++i)
				fscanf(fp, "%u", P + i);

			for(const index_t & v: trig_convex_polygon(P, nv))
				faces.push_back(v);
		}
	}
	else // binary_little_endian or binary_big_endian
	{
		bool big_endian = format[7] == 'b';
		auto big_to_little = [](char * buffer, const index_t & n)
		{
			for(index_t i = 0, j = n - 1; i < j; ++i, --j)
				swap(buffer[i], buffer[j]);
		};

		char * buffer = vbytes == sizeof(vertex) ? (char *) GT : new char[vbytes * n_vertices];

		fread(buffer, vbytes, n_vertices, fp);
		if(big_endian)
		{
			char * pb;

			#pragma omp parallel for private(pb)
			for(index_t v = 0; v < n_vertices; ++v)
			{
				pb = buffer + v * vbytes;

				for(index_t i = 0; i < 3; ++i)
					big_to_little(pb + ixyz + i * xyz, xyz);

				if(rgb)
				{
					for(index_t i = 0; i < 3; ++i)
						big_to_little(pb + irgb + i * rgb, rgb);
				}
			}
		}

		if(vbytes != sizeof(vertex))
		{
			char * pb;

			#pragma omp parallel for private(pb)
			for(index_t v = 0; v < n_vertices; ++v)
			{
				pb = buffer + v * vbytes;

				for(index_t i = 0; i < 3; ++i)
				{
					if(xyz == 4)
						GT[v][i] = (real_t) *(float *) (pb + ixyz + i * xyz);
					else
						GT[v][i] = (real_t) *(double *) (pb + ixyz + i * xyz);
				}

				if(rgb)
				{
					for(index_t i = 0; i < 3; ++i)
						if(rgb == 1)
							VC[v][i] = *(unsigned char *) (pb + irgb + i * rgb);
						else
							GT[v][i] = (unsigned char) (*(float *) (pb + irgb + i * rgb)) * 255;
				}
			}
		}

		if(buffer != (char *) GT) delete [] buffer;


		while(nf--)
		{
			char buffer[8];

			fread(buffer, 1, fn, fp);
			if(big_endian) big_to_little(buffer, fn);

			if(fn == 1) nv = *((char *) buffer);
			if(fn == 2) nv = *((short *) buffer);
			if(fn == 4) nv = *((int *) buffer);

			for(index_t i = 0; i < nv; ++i)
			{
				fread(buffer, 1, fbytes, fp);
				if(big_endian) big_to_little(buffer, fbytes);

				if(fbytes == 1) P[i] = *((char *) buffer);
				if(fbytes == 2) P[i] = *((short *) buffer);
				if(fbytes == 4) P[i] = *((int *) buffer);
			}

			for(const index_t & v: trig_convex_polygon(P, nv))
				faces.push_back(v);
		}
	}

	fclose(fp);


	if(faces.size() != che::mtrig * n_faces)
	{
		vertex * tGT = GT; GT = nullptr;
		rgb_t * tVC = VC; VC = nullptr;

		free();
		alloc(nv, faces.size() / che::mtrig);

		GT = tGT;
		VC = tVC;
	}

	memcpy(VT, faces.data(), faces.size() * sizeof(index_t));
}

void che_ply::write_file(const che * mesh, const string & file, const bool & color)
{
	FILE * fp = fopen((file + ".ply").c_str(), "wb");
	assert(fp);

	const char * type = sizeof(real_t) == 4 ? "float" : "double";

	fprintf(fp, "ply\n");
	fprintf(fp, "format binary_little_endian 1.0\n");
	fprintf(fp, "comment PLY generated by gproshan\n");
	fprintf(fp, "element vertex %lu\n", mesh->n_vertices);
	fprintf(fp, "property %s x\n", type);
	fprintf(fp, "property %s y\n", type);
	fprintf(fp, "property %s z\n", type);
	if(color)
	{
		fprintf(fp, "property uchar red\n");
		fprintf(fp, "property uchar green\n");
		fprintf(fp, "property uchar blue\n");
	}
	fprintf(fp, "element face %lu\n", mesh->n_faces);
	fprintf(fp, "property list uchar uint vertex_index\n");
	fprintf(fp, "end_header\n");

	if(color)
	{
		for(index_t v = 0; v < mesh->n_vertices; ++v)
		{
			fwrite(&mesh->gt(v), sizeof(vertex), 1, fp);
			fwrite(&mesh->rgb(v), sizeof(rgb_t), 1, fp);
		}
	}
	else fwrite(&mesh->gt(0), sizeof(vertex), mesh->n_vertices, fp);

	unsigned char mtrig = che::mtrig;
	for(index_t he = 0; he < mesh->n_half_edges; he += che::mtrig)
	{
		fwrite(&mtrig, 1, 1, fp);
		fwrite(&mesh->vt(he), sizeof(index_t), che::mtrig, fp);
	}

	fclose(fp);
}


} // namespace gproshan
