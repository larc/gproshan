#include "raytracing/rt_embree_splat_ch.h"

#ifdef GPROSHAN_EMBREE

#include "mesh/che_off.h"

#include <cstring>
#include <set>
#include <queue>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


float embree_splat_ch::r_threshold = 0.95;
float embree_splat_ch::n_threshold = 0.45;
size_t embree_splat_ch::max_neigs = 64;

embree_splat_ch::embree_splat_ch(const std::vector<che *> & meshes, const bool & pointcloud)
{
	build_bvh(meshes, pointcloud);
}

index_t embree_splat_ch::add_pointcloud(const che * mesh)
{
	init_splats(mesh);

	std::vector<index_t> vstart(vsplat.size() + 1);
	std::vector<convex_hull *> vch(vsplat.size());

	vstart[0] = 0;
	for(index_t i = 1; i < vstart.size(); ++i)
		vstart[i] = vstart[i - 1] + vsplat[i - 1].points.size();

	std::vector<vertex> vertices(vstart.back());
	std::vector<index_t> faces;

	#pragma omp parallel for
	for(index_t i = 0; i < vsplat.size(); ++i)
	{
		splat & is = vsplat[i];
		const index_t & begin = vstart[i];
		const index_t & end = vstart[i + 1];

		std::vector<index_t> & points = is;

		is.n = 0;
		for(index_t j = 0; j < points.size(); ++j)
		{
			index_t & v = points[j];
			vertices[j + begin] = mesh->gt(v);
			is.n += mesh->normal(v);
		}

		is.n = is.n.unit();
		is.c = vertices[begin];
		is.t = vertices[end - 1] - is.c;
		is.t = (is.t - ((is.t, is.n) * is.n)).unit();
		is.b = (is.n * is.t).unit();

		for(index_t j = begin; j < end; ++j)
			is.to2d(vertices[j]);

		if(points.size() >= 3)
			vch[i] = new convex_hull(vertices.data() + begin, points.size());
		else
			vch[i] = nullptr;

		for(index_t j = begin; j < end; ++j)
			is.to3d(vertices[j]);
	}

	for(index_t i = 0; i < vsplat.size(); ++i)
	{
		if(!vch[i]) continue;

		const std::vector<index_t> & ch = *vch[i];

		index_t f = 0;
		for(index_t & v: che::trig_convex_polygon(ch.data(), ch.size()))
		{
			faces.push_back(v + vstart[i]);
			if(!(f % 3)) primID_splat.push_back(i);
			++f;
		}
	}

	for(convex_hull * ch: vch)
		delete ch;

	che ch_mesh(vertices.data(), vertices.size(), faces.data(), faces.size() / 3);
	che_off::write_file(&ch_mesh, "ch_splats");
	return add_mesh(&ch_mesh);
}

float embree_splat_ch::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	float w = vsplat[primID_splat[r.hit.primID]].shading(geomID_mesh[r.hit.geomID], position, normal, color);

	return 1e-2;
}

void embree_splat_ch::init_splats(const che * mesh)
{
	vsplat.reserve(mesh->n_vertices);

	std::vector<bool> visited;
	visited.assign(mesh->n_vertices, 0);


	real_t radio;
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(visited[v]) continue;

		std::set<index_t> points;
		std::queue<index_t> q;

		q.push(v);
		points.insert(v);

		index_t u;
		while(!q.empty() && points.size() < max_neigs)
		{
			for_star(he, mesh, q.front())
			{
				u = mesh->vt(prev(he));
				if(points.find(u) == points.end())
				{
					points.insert(u);
					q.push(u);
				}
			}

			q.pop();
		}

		const vertex & n = mesh->normal(v);
		const vertex & c = mesh->gt(v);

		radio = 0;

		vsplat.push_back(splat());
		std::vector<index_t> & splat_points = vsplat.back();
		for(const index_t & p: points)
			if((n, mesh->normal(p)) > n_threshold)
			{
				splat_points.push_back(p);
				radio = std::max(radio, *(mesh->gt(p) - c));
			}
			else break;

		radio *= r_threshold;
		for(const index_t & p: splat_points)
			visited[p] = *(mesh->gt(p) - c) < radio;

		if(splat_points.size() < 3)
			vsplat.pop_back();
	}

	gproshan_error_var(float(vsplat.size()) / mesh->n_vertices);
}


// FROM: https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int expand_bits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
// UPDATED ONLY 2D
unsigned int morton_2d(float x, float y)
{
    x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
    y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx = expand_bits(x);
    unsigned int yy = expand_bits(y);
    return (xx >> 2) + (yy >> 1);
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

