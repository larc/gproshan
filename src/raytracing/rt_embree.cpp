#include "rt_embree.h"

#ifdef GPROSHAN_EMBREE

#include <iostream>
#include <random>
#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


void compute_normals(glm::vec3 * normals, const glm::vec3 * vertices, const size_t & n_vertices)
{

}

void embree_error(void * ptr, RTCError error, const char * str)
{
	fprintf(stderr, "EMBREE ERROR: %s\n", str);
}

embree::embree(const std::vector<che *> & meshes)
{
	device = rtcNewDevice(NULL);
	rtcSetDeviceErrorFunction(device, embree_error, NULL);

	scene = rtcNewScene(device);
	
	rtcInitIntersectContext(&intersect_context);

	build_bvh(meshes);
}

embree::~embree()
{
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);
}

void embree::build_bvh(const std::vector<che *> & meshes)
{
	for(auto & m: meshes)
		if(m->n_faces()) geomID_mesh[add_mesh(m)] = m;
		else geomID_mesh[add_point_cloud(m)] = m;

	rtcCommitScene(scene);
}

index_t embree::add_sphere(const glm::vec4 & xyzr)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

	glm::vec4 * pxyzr = (glm::vec4 *) rtcSetNewGeometryBuffer(geom,
															RTC_BUFFER_TYPE_VERTEX, 0,
															RTC_FORMAT_FLOAT4, 4 * sizeof(float), 1);
	*pxyzr = xyzr;

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
	
	return geom_id;
}

index_t embree::add_mesh(const che * mesh)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

	glm::vec3 * vertices = (glm::vec3 *) rtcSetNewGeometryBuffer(geom,
																RTC_BUFFER_TYPE_VERTEX, 0,
																RTC_FORMAT_FLOAT3, 3 * sizeof(float),
																mesh->n_vertices()
																);

	index_t * tri_idxs = (index_t *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_INDEX, 0,
																RTC_FORMAT_UINT3, 3 * sizeof(index_t),
																mesh->n_faces()
																);
#ifdef SINGLE_P
	memcpy(vertices, &mesh->gt(0), mesh->n_vertices() * sizeof(vertex));
#else
	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
		vertices[i] = glm::vec3(mesh->gt(i).x, mesh->gt(i).y, mesh->gt(i).z);
#endif // SINGLE_P
	
	memcpy(tri_idxs, &mesh->vt(0), mesh->n_half_edges() * sizeof(index_t));

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
	
	return geom_id;
}

index_t embree::add_point_cloud(const che * mesh)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

	glm::vec4 * pxyzr = (glm::vec4 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_VERTEX, 0,
																RTC_FORMAT_FLOAT4,
																4 * sizeof(float),
																mesh->n_vertices()
																);
	
	glm::vec3 * normal = (glm::vec3 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_NORMAL, 0,
																RTC_FORMAT_FLOAT3,
																3 * sizeof(float),
																mesh->n_vertices()
																);
	
	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		pxyzr[i] = glm::vec4(mesh->gt(i).x, mesh->gt(i).y, mesh->gt(i).z, 0.001);
		normal[i] = glm::vec3(mesh->normal(i).x, mesh->normal(i).y, mesh->normal(i).z);
	}

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

void embree::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, ray_hit & r)
{
	float dist, sum_dist = 0;
	position = glm::vec3(0);
	normal = glm::vec3(0);

	do
	{
		glm::vec4 * xyzr = (glm::vec4 *) rtcGetGeometryBufferData(rtcGetGeometry(scene, r.hit.geomID), RTC_BUFFER_TYPE_VERTEX, 0);
		
		sum_dist += dist = glm::length(r.position() - glm::vec3(xyzr[r.hit.primID]));
		position += dist * r.position();
		normal += dist * r.normal(geomID_mesh[r.hit.geomID], true);

		r = ray_hit(r.position(), r.dir());
	}
	while(intersect(r) && r.ray.tfar < 0.001);

	position /= sum_dist;
	normal /= sum_dist;
}

glm::vec4 embree::li(const glm::vec3 & light, const glm::vec3 & position, const glm::vec3 & normal)
{
	const glm::vec3 color(0.6, 0.8, 1.0);
	
	const glm::vec3 wi = normalize(light - position);
	const float dist_light = glm::length(light - position);
	const float falloff = 4.f / (dist_light * dist_light);	// intensity multiplier / falloff
	const float dot_wi_normal = glm::dot(wi, normal);

	const glm::vec4 L = glm::vec4(falloff * (dot_wi_normal < 0 ? -dot_wi_normal : dot_wi_normal) * color, 1);

	ray_hit r(position, wi);
	
	return occluded(r) ? 0.5f * L : L;
}

glm::vec4 embree::li(ray_hit r, const glm::vec3 & light, const bool & flat)
{
	const float max_tfar = 8;

	float total_tfar = 0;
	float tfar = r.ray.tfar;
	glm::vec4 L(0);
	
	glm::vec3 position, normal;
	while(tfar < max_tfar)
	{
		total_tfar += tfar;
		
		position = r.position();
		normal = r.normal(geomID_mesh[r.hit.geomID], flat);
		
		if(geomID_mesh[r.hit.geomID]->is_pointcloud())
			pointcloud_hit(position, normal, r);
		
		L += li(light, position, normal);

		/*
		r = ray_hit(r.position(), r.dir());
		if(!intersect(r)) break;

		tfar = r.ray.tfar + total_tfar;
		*/
		
		tfar = 10;
	}

	return L / total_tfar;
}

bool embree::intersect(ray_hit & r)
{
	rtcIntersect1(scene, &intersect_context, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}

bool embree::occluded(ray_hit & r)
{
	rtcIntersect1(scene, &intersect_context, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}


glm::vec4 embree::intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light,const bool & flat)
{
	ray_hit r(org, dir);
	return intersect(r) ? li(r, light, flat) : glm::vec4(0.f);
}

float embree::intersect_depth(const glm::vec3 & org, const glm::vec3 & dir)
{
	ray_hit r(org, dir);
	return intersect(r) ? r.ray.tfar : 0.f;
}


} // namespace gproshan

#endif // GPROSHAN_EMBREE

