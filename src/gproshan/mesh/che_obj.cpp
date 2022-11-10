#include <gproshan/mesh/che_obj.h>

#include <cstring>
#include <cstdio>
#include <cassert>


// geometry processing and shape analysis framework
namespace gproshan {


che_obj::che_obj(const string & file)
{
	init(file);
}

void che_obj::read_file(const string & file)
{
	FILE * fp = fopen(file.c_str(), "r");
	assert(fp);

	float x, y, z, r, g, b;
	index_t P[32], n;

	vector<vertex> vertices;
	vector<rgb_t> vertices_color;
	vector<index_t> faces;

	char line[256], str[64];
	char * line_ptr;
	index_t offset;

	while(fgets(line, sizeof(line), fp))
	{
		str[0] = 0;
		line_ptr = line;

		sscanf(line_ptr, "%s%n", str, &offset);
		line_ptr += offset;

		if(str[0] == 'v' && !str[1])	// v x y z
		{
			n = sscanf(line_ptr, "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
			vertices.push_back({x, y, z});
			vertices_color.push_back(n == 6 ? rgb_t{(unsigned char) (r * 255), (unsigned char) (g * 255), (unsigned char) (b * 255)} : rgb_t());
		}

		if(str[0] == 'f')				// f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
		{
			n = 0;
			while(sscanf(line_ptr, "%s%n", str, &offset) > 0)
			{
				line_ptr += offset;
				sscanf(str, "%d%*s", P + n);
				P[n] += P[n] > vertices.size() ? vertices.size() : -1;
				++n;
			}

			for(const index_t & v: trig_convex_polygon(P, n))
				faces.push_back(v);
		}
	}

	fclose(fp);


	alloc(vertices.size(), faces.size() / che::mtrig);
	memcpy(GT, vertices.data(), vertices.size() * sizeof(vertex));
	memcpy(VC, vertices_color.data(), vertices_color.size() * sizeof(rgb_t));
	memcpy(VT, faces.data(), faces.size() * sizeof(index_t));
}

void che_obj::write_file(const che * mesh, const string & file, const bool & color, const bool & pointcloud)
{
	FILE * fp = fopen((file + ".obj").c_str(), "w");
	assert(fp);

	fprintf(fp, "# OBJ generated by gproshan\n");
	fprintf(fp, "# vertices %lu\n", mesh->n_vertices);
	fprintf(fp, "# faces %lu\n", mesh->n_faces);

	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		const vertex & v = mesh->point(i);
		fprintf(fp, "v %f %f %f", (float) v.x(), (float) v.y(), (float) v.z());
		if(color)
		{
			const vertex & c = mesh->color(i);
			fprintf(fp, " %f %f %f", (float) c.x(), (float) c.y(), (float) c.z());
		}
		fprintf(fp, "\n");
	}

	if(!pointcloud)
	{
		for(index_t he = 0; he < mesh->n_half_edges; )
		{
			fprintf(fp, "f");
			for(index_t i = 0; i < che::mtrig; ++i)
				fprintf(fp, " %u", mesh->halfedge(he++) + 1);
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}


} // namespace gproshan

