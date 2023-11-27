#include <gproshan/mesh/che_ply.h>

#include <cstring>
#include <cstdio>
#include <cassert>
#include <unordered_map>


// geometry processing and shape analysis framework
namespace gproshan {


che_ply::che_ply(const std::string & file)
{
	init(file);
}

void che_ply::read_file(const std::string & file)
{
	std::unordered_map<std::string, size_t> bytes = {
														{"char"		, 1},
														{"uchar"	, 1},
														{"short"	, 2},
														{"ushort"	, 2},
														{"int"		, 4},
														{"uint"		, 4},
														{"float"	, 4},
														{"float32"	, 4},
														{"float64"	, 8},
														{"double"	, 8}
													};

	FILE * fp = fopen(file.c_str(), "rb");
	assert(fp);

	size_t nv = 0;
	size_t nf = 0;
	size_t xyz = 0;
	size_t rgb = 0;
	size_t normal = 0;
	size_t vbytes = 0;
	index_t ixyz = 0, irgb = 0, inormal = 0;
	size_t fn = 0, fbytes = 0;
	index_t P[32];

	char line[512], type[32], str[32], format[32], element[32];

	auto add_vproperty = [&vbytes](size_t & p, index_t & i, const size_t & nb)
	{
		p = nb;
		i = vbytes;
	};

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
			const size_t & nb = bytes[type];

			switch(str[0])
			{
				case 'x': add_vproperty(xyz, ixyz, nb);
					break;
				case 'r': add_vproperty(rgb, irgb, nb);
					break;
				case 'n':
					if(str[1] == 'x')
						add_vproperty(normal, inormal, nb);
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

	std::vector<index_t> trigs;
	trigs.reserve(che::mtrig * n_trigs);

	if(format[0] == 'a')	// ascii
	{
		float x, y, z, nx, ny, nz;
		unsigned char r, g, b;
		for(index_t v = 0; v < n_vertices; ++v)
		{
			fgets(line, sizeof(line), fp);

			rgb ? sscanf(line, "%f %f %f %hhu %hhu %hhu %f %f %f", &x, &y, &z, &r, &g, &b, &nx, &ny, &nz)
				: sscanf(line, "%f %f %f %f %f %f", &x, &y, &z, &nx, &ny, &nz);

			GT[v] = {x, y, z};
			if(rgb) VC[v] = {r, g, b};
			if(normal) VN[v] = {nx, ny, nz};
		}

		while(nf--)
		{
			fscanf(fp, "%lu", &nv);
			for(index_t i = 0; i < nv; ++i)
				fscanf(fp, "%u", P + i);

			for(const index_t & v: trig_convex_polygon(P, nv))
				trigs.push_back(v);
		}
	}
	else // binary_little_endian or binary_big_endian
	{
		bool big_endian = format[7] == 'b';
		auto big_to_little = [](char * buffer, const index_t & n)
		{
			for(index_t i = 0, j = n - 1; i < j; ++i, --j)
				std::swap(buffer[i], buffer[j]);
		};

		char * buffer = vbytes == sizeof(vertex) ? (char *) GT : new char[vbytes * n_vertices];

		fread(buffer, vbytes, n_vertices, fp);
		if(big_endian)
		{
			#pragma omp parallel for
			for(index_t v = 0; v < n_vertices; ++v)
			{
				char * pb = buffer + v * vbytes;

				for(index_t i = 0; i < 3; ++i)
					big_to_little(pb + ixyz + i * xyz, xyz);

				if(rgb)
				{
					for(index_t i = 0; i < 3; ++i)
						big_to_little(pb + irgb + i * rgb, rgb);
				}

				if(normal)
				{
					for(index_t i = 0; i < 3; ++i)
						big_to_little(pb + inormal + i * normal, normal);
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
							VC[v][i] = (unsigned char) (*(float *) (pb + irgb + i * rgb)) * 255;
				}

				if(normal)
				{
					for(index_t i = 0; i < 3; ++i)
					{
						if(normal == 4)
							VN[v][i] = (real_t) *(float *) (pb + inormal + i * normal);
						else
							VN[v][i] = (real_t) *(double *) (pb + inormal + i * normal);
					}
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
				trigs.push_back(v);
		}
	}

	fclose(fp);


	if(size(trigs) != che::mtrig * n_trigs)
	{
		vertex * tGT = GT; GT = nullptr;
		rgb_t * tVC = VC; VC = nullptr;

		free();
		alloc(nv, trigs.size() / che::mtrig);

		GT = tGT;
		VC = tVC;
	}

	memcpy(VT, trigs.data(), trigs.size() * sizeof(index_t));
}

void che_ply::write_file(const che * mesh, const std::string & file, const bool & color)
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
	fprintf(fp, "property %s nx\n", type);
	fprintf(fp, "property %s ny\n", type);
	fprintf(fp, "property %s nz\n", type);
	if(color)
	{
		fprintf(fp, "property uchar red\n");
		fprintf(fp, "property uchar green\n");
		fprintf(fp, "property uchar blue\n");
	}
	fprintf(fp, "element face %lu\n", mesh->n_trigs);
	fprintf(fp, "property list uchar uint vertex_index\n");
	fprintf(fp, "end_header\n");

	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		fwrite(&mesh->point(v), sizeof(vertex), 1, fp);
		fwrite(&mesh->normal(v), sizeof(vertex), 1, fp);
		if(color) fwrite(&mesh->rgb(v), sizeof(rgb_t), 1, fp);
	}

	unsigned char mtrig = che::mtrig;
	for(index_t he = 0; he < mesh->n_half_edges; he += che::mtrig)
	{
		fwrite(&mtrig, 1, 1, fp);
		fwrite(&mesh->halfedge(he), sizeof(index_t), che::mtrig, fp);
	}

	fclose(fp);
}


} // namespace gproshan

