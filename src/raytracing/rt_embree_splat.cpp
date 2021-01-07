#include "raytracing/rt_embree_splat.h"

#ifdef GPROSHAN_EMBREE


#include <set>
#include <queue>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


embree_splat::embree_splat(const std::vector<che *> & meshes, const bool & pointcloud)
{
	build_bvh(meshes, pointcloud);
}

index_t embree_splat::add_pointcloud(const che * mesh)
{
	init_splats(mesh);

	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

	glm::vec4 * pxyzr = (glm::vec4 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_VERTEX, 0,
																RTC_FORMAT_FLOAT4,
																4 * sizeof(float),
																vsplat.size()
																);
	
	glm::vec3 * normal = (glm::vec3 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_NORMAL, 0,
																RTC_FORMAT_FLOAT3,
																3 * sizeof(float),
																vsplat.size()
																);
	
	#pragma omp parallel for
	for(index_t i = 0; i < vsplat.size(); i++)
	{
		pxyzr[i] = vsplat[i].xyzr();
		normal[i] = vsplat[i].normal();
	}

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
	
	return geom_id;
}

float embree_splat::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	vsplat[r.hit.primID].shading(position, normal, color);
	// normal = vsplat[r.hit.primID].normal();
	// color = vsplat[r.hit.primID].color();

	return 1e-2f;
}

void embree_splat::init_splats(const che * mesh)
{
	const size_t n = 5;
	vsplat.resize((mesh->n_vertices() + n - 1) / n);

	gproshan_log_var(vsplat.size());
	
	#pragma omp parallel for
	for(index_t i = 0; i < vsplat.size(); i++)
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
		
		real_t dist, d;
		const vertex & c = mesh->gt(v);
		
		splat & s = vsplat[i];
		for(index_t j = 0; j < K; j++)
		{
			dist = INFINITY;

			for(const index_t & p: points)
			{
				d = *(mesh->gt(p) - c);
				if(d < dist)
				{
					dist = d;
					u = p;
				}
			}

			points.erase(u);

			s.P[j] = glm::vec3(mesh->gt(u).x, mesh->gt(u).y, mesh->gt(u).z);
			s.N[j] = glm::vec3(mesh->normal(u).x, mesh->normal(u).y, mesh->normal(u).z);
			s.C[j] = glm::vec3(mesh->color(u).x, mesh->color(u).y, mesh->color(u).z);
		}

		s.radio = glm::length(s.P[K - 1] - s.P[0]);
	}
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

