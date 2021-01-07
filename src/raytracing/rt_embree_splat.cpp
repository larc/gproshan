#include "raytracing/rt_embree_splat.h"

#ifdef GPROSHAN_EMBREE


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
		pxyzr[i] = vsplat[i].xyzr;
		normal[i] = vsplat[i].normal;
	}

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
	
	return geom_id;
}

float embree_splat::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	position = r.position();
	normal = vsplat[r.hit.primID].normal;
	color = vsplat[r.hit.primID].color;

	return 1e-5f;
}

void embree_splat::init_splats(const che * mesh)
{
	const size_t s = 10;
	vsplat.reserve(vsplat.size() + (mesh->n_vertices() + s - 1) / s);

	for(index_t i = 0; i < mesh->n_vertices(); i += s)
		vsplat.push_back({	glm::vec4(mesh->gt(i).x, mesh->gt(i).y, mesh->gt(i).z, pc_radius),
							glm::vec3(mesh->normal(i).x, mesh->normal(i).y, mesh->normal(i).z),
							glm::vec3(mesh->color(i).x, mesh->color(i).y, mesh->color(i).z)
							});
						
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

