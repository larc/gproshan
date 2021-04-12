#include "raytracing/rt_embree_splat_ch.h"

#ifdef GPROSHAN_EMBREE

#include "mesh/che_off.h"

#include <cstring>
#include <set>
#include <queue>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


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
		const index_t & begin = vstart[i];
		const index_t & end = vstart[i + 1];

		std::vector<index_t> & points = vsplat[i];

		vertex c, t, b, n = 0;

		for(index_t j = 0; j < points.size(); ++j)
		{
			index_t & v = points[j];
			vertices[j + begin] = mesh->gt(v);
			n += mesh->normal(v);
		}

		n = n.unit();
		c = vertices[begin];
		t = vertices[end - 1] - c;
		t = (t - ((t, n) * n)).unit();
		b = (n * t).unit();

		for(index_t j = begin; j < end; ++j)
		{
			vertex & v = vertices[j];
			v -= c;
			v = {(t, v), (b, v), (n, b)};
		}

		if(points.size() >= 3)
			vch[i] = new convex_hull(vertices.data() + begin, points.size());
		else
			vch[i] = nullptr;

		for(index_t j = begin; j < end; ++j)
		{
			vertex & v = vertices[j];
			v = vertex{	(vertex{t.x, b.x, n.x}, v),
						(vertex{t.y, b.y, n.y}, v),
						(vertex{t.z, b.z, n.z}, v)
						} + c;
		}
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

	return add_mesh(&ch_mesh);
}

float embree_splat_ch::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	float w = vsplat[primID_splat[r.hit.primID]].shading(geomID_mesh[r.hit.geomID], position, normal, color);

	if(w < 1e-5f)
	{
		r = ray_hit(r.position(), r.dir());
		if(intersect(r))
			return pointcloud_hit(position, normal, color, r);
	}

	return 1e-2f;
}

void embree_splat_ch::init_splats(const che * mesh)
{
	const size_t n = 100;
	vsplat.resize((mesh->n_vertices + n - 1) / n);

	gproshan_log_var(vsplat.size());

	#pragma omp parallel for
	for(index_t i = 0; i < vsplat.size(); ++i)
	{
		const index_t v = n * i;	// random, feature aware index

		std::set<index_t> points;
		std::queue<index_t> q;

		q.push(v);
		points.insert(v);

		index_t u;
		while(!q.empty() && points.size() < K)
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

		const vertex & n = mesh->normal(*points.begin());

		std::vector<index_t> & s = vsplat[i];
		for(const index_t & p: points)
			if((n, mesh->normal(p)) > 0.85)
				s.push_back(p);
			else break;
	}
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

