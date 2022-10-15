#include "raytracing/rt_embree_splat_ch.h"


#include "mesh/che_off.h"

#include <cstring>
#include <set>
#include <queue>
#include <algorithm>
#include <numeric>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


bool embree_splat_ch::show_chsplats = true;
int embree_splat_ch::k_neighbors = 4;
float embree_splat_ch::r_threshold = 0.50;	// cos overlapping radius
float embree_splat_ch::n_threshold = 0.81;	// 30 degrees angle normals
size_t embree_splat_ch::max_neighbors = 256;	// max neighbors per splat


glm::vec3 colormap(const float & x)
{
	float r = x < 0.75 ? 1012.0 * x - 389.0 : -1.11322769567548E+03 * x + 1.24461193212872E+03;
	float g = x < 0.50 ? 1012.0 * x - 129.0 : -1012.0 * x + 899.0;
	float b = x < 0.25 ? 1012.0 * x + 131.0 : -1012.0 * x + 643.0;
	r = std::min(std::max(r / 255.0, 0.0), 1.0);
	g = std::min(std::max(g / 255.0, 0.0), 1.0);
	b = std::min(std::max(b / 255.0, 0.0), 1.0);
	return glm::vec3(r, g, b);
}


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
		vstart[i] = vstart[i - 1] + vsplat[i - 1].size();

	std::vector<vertex> vertices(vstart.back());
	std::vector<index_t> faces;

	#pragma omp parallel for
	for(index_t i = 0; i < vsplat.size(); ++i)
	{
		splat & is = vsplat[i];

		vch[i] = nullptr;
		if(is.size() < 3) continue;

		const index_t & begin = vstart[i];
		const index_t & end = vstart[i + 1];

		is.c = is.n = 0;
		for(index_t j = 0; j < is.size(); ++j)
		{
			index_t & v = is[j];
			vertices[j + begin] = mesh->gt(v);
			is.n += mesh->normal(v);
			is.c += mesh->gt(v);
		}

		is.c /= is.size();
		is.n = is.n.unit();
		is.t = vertices[end - 1] - is.c;
		is.t = (is.t - ((is.t, is.n) * is.n)).unit();
		is.b = (is.n * is.t).unit();

		for(index_t j = begin; j < end; ++j)
		{
			vertex & v = vertices[j];
			is.to2d(vertices[j]);
			is.code(j - begin) = morton_2d((v.x + 1) / 2, (v.y + 1) / 2);
		}

		// sorting point by its morton code
		std::sort(is.ipoints.begin(), is.ipoints.end());

		vch[i] = new convex_hull(vertices.data() + begin, is.size());

		for(index_t j = begin; j < end; ++j)
			is.to3d(vertices[j]);
	}

	csplat.resize(vsplat.size());
	for(index_t i = 0; i < vsplat.size(); ++i)
	{
		csplat[i] = float(i) / (vsplat.size() - 1);
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

	std::random_shuffle(csplat.begin(), csplat.end());

	for(convex_hull * ch: vch)
		delete ch;

	che ch_mesh(vertices.data(), vertices.size(), faces.data(), faces.size() / 3);
	che_off::write_file(&ch_mesh, "ch_splats");

	return add_mesh(&ch_mesh);
}

float embree_splat_ch::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	if(show_chsplats)
	{
		normal = glm::normalize(glm::vec3(r.hit.Ng_x, r.hit.Ng_y, r.hit.Ng_z));
		color = colormap(csplat[primID_splat[r.hit.primID]]);
	}
	else
	{
		vsplat[primID_splat[r.hit.primID]].shading(geomID_mesh[r.hit.geomID], position, normal, color);
	}

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

		std::set<index_t> neigs;
		std::queue<index_t> q;

		q.push(v);
		neigs.insert(v);
		vertex n = mesh->normal(v);

		while(!q.empty() && neigs.size() < max_neighbors)
		{
			for_star(he, mesh, q.front())
			{
				const index_t & u = mesh->vt(prev(he));
				if(!visited[u] && neigs.find(u) == neigs.end())
				{
					if((n.unit(), mesh->normal(u)) > n_threshold)
					{
						q.push(u);
						n += mesh->normal(u);
						visited[u] = true;
					}
					neigs.insert(u);
				}
			}
			q.pop();
		}

		if(neigs.size() < 3) continue;

		vsplat.push_back(splat());
		splat & s = vsplat.back();

		for(const index_t & p: neigs)
			s.push_back(p);
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
    return (xx >> 1) + yy;
}


} // namespace gproshan

