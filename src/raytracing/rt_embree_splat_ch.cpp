#include "raytracing/rt_embree_splat_ch.h"

#ifdef GPROSHAN_EMBREE


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

	return add_mesh(mesh);
}

float embree_splat_ch::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	float w = vsplat[r.hit.primID].shading(geomID_mesh[r.hit.geomID], position, normal, color);
	// normal = vsplat[r.hit.primID].normal();
	// color = vsplat[r.hit.primID].color();
/*	if(w < 1e-2f)
	{
		normal = glm::vec3(0);
		color = glm::vec3(0);
	}
*/
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
	const size_t n = 10;
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
		while(!q.empty() && points.size() < 2 * K)
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

		std::vector<index_t> & s = vsplat[i];
		for(const index_t & p: points)
			s.push_back(p);
	}
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

