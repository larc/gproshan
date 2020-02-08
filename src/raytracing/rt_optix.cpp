#include "rt_optix.h"

#ifdef GPROSHAN_OPTIX

#include <iostream>
#include <random>
#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


void optix_error(void * ptr, RTCError error, const char * str)
{
	fprintf(stderr, "OPTIX ERROR: %s\n", str);
}

optix::optix()
{
	device = rtcNewDevice(NULL);
	rtcSetDeviceErrorFunction(device, optix_error, NULL);

	scene = rtcNewScene(device);
	
	rtcInitIntersectContext(&intersect_context);
}

optix::~optix()
{
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);
}

void optix::build_bvh()
{
	rtcCommitScene(scene);
}

index_t optix::add_sphere(const glm::vec4 & xyzr)
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

index_t optix::add_mesh(const che * mesh, const glm::mat4 & model_matrix)
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
	
	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices(); i++)
	{
		glm::vec4 v(mesh->gt(i).x, mesh->gt(i).y, mesh->gt(i).z, 1.f);
		vertices[i] = glm::vec3(model_matrix * v);
	}
	
	memcpy(tri_idxs, &mesh->vt(0), mesh->n_half_edges() * sizeof(index_t));

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
	
	return geom_id;
}

index_t optix::add_point_cloud(const che * mesh, const glm::mat4 & model_matrix)
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
		pxyzr[i] = glm::vec4(mesh->gt(i).x, mesh->gt(i).y, mesh->gt(i).z, 0.001f);
		
		vertex n = mesh->normal(i);
		normal[i] = glm::vec3(n.x, n.y, n.z);
	}

	rtcCommitGeometry(geom);
	
	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

glm::vec4 optix::li(const ray_hit & r, const glm::vec3 & light)
{
	glm::vec3 color(.6f, .8f, 1.f);
	
	float dist_light = glm::length(light - r.position());
	float falloff = 4.f / (dist_light * dist_light);	// intensity multiplier / falloff
	
	glm::vec3 wi = normalize(light - r.position());
	
	float dot_wi_normal = glm::dot(wi, r.normal());
	if(dot_wi_normal < 0)
		dot_wi_normal = -dot_wi_normal;

	ray_hit ro(r.position() + 1e-5f * wi, wi);

	if(occluded(ro))
		return .5f * glm::vec4(color * falloff * dot_wi_normal, 1.f);
	
	return glm::vec4(color * falloff * dot_wi_normal, 1.f);
}

bool optix::intersect(ray_hit & r)
{
	rtcIntersect1(scene, &intersect_context, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}

bool optix::occluded(ray_hit & r)
{
	rtcIntersect1(scene, &intersect_context, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}


const glm::vec4 optix::intersect_li(const glm::vec3 & org, const glm::vec3 & dir, const glm::vec3 & light)
{
	ray_hit r(org, dir);
	return intersect(r) ? li(r, light) : glm::vec4(0.f);
}

const float optix::intersect_depth(const glm::vec3 & org, const glm::vec3 & dir)
{
	ray_hit r(org, dir);
	return intersect(r) ? r.ray.tfar : 0.f;
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

