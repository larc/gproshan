#include "embree.h"

#ifdef GPROSHAN_EMBREE

#include <iostream>
#include <random>
#include <cstring>

void embree_error(void * ptr, RTCError error, const char * str)
{
	fprintf(stderr, "EMBREE ERROR: %s\n", str);
}

embree::embree()
{
	device = rtcNewDevice(NULL);
	rtcSetDeviceErrorFunction(device, embree_error, NULL);

	scene = rtcNewScene(device);
	
	rtcInitIntersectContext(&intersect_context);

	width = height = n_samples = 0;
	img = nullptr;
}

embree::~embree()
{
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);

	if(img) delete [] img;
}

void embree::build_bvh()
{
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

index_t embree::add_mesh(const che * mesh, const glm::mat4 & model_matrix)
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

index_t embree::add_point_cloud(const che * mesh, const glm::mat4 & model_matrix)
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

glm::vec4 embree::Li(const ray_hit & r, const glm::vec3 & light)
{
	glm::vec3 color(.6f, .8f, 1.f);
	
	float dist_light = glm::length(light - r.position());
	float falloff = 4.f / (dist_light * dist_light);	// intensity multiplier / falloff
	
	glm::vec3 wi = normalize(light - r.position());
	
	ray_hit ro(r.position() + 1e-5f * wi, wi);
	if(occluded(ro))
		return glm::vec4(.2f * color * falloff * std::max(0.f, glm::dot(wi, r.normal())), 1.f);

	return glm::vec4(color * falloff * std::max(0.f, glm::dot(wi, r.normal())), 1.f);
}

bool embree::rt_restart(const size_t & w, const size_t & h)
{
	if(width * height < w * h)
	{
		if(img) delete [] img;

		width = w;
		height = h;
		img = new glm::vec4[width * height];

		return true;
	}
	
	if(width != w || height != h)
	{
		width = w;
		height = h;

		return true;
	}

	return false;
}

void embree::raytracing(	const glm::uvec2 & windows_size,
							const glm::mat4 & view_mat,
							const glm::mat4 & proj_mat,
							const glm::vec3 & light, const bool & restart )
{
	if(rt_restart(windows_size.x, windows_size.y) || restart)
		n_samples = 0;

	std::default_random_engine gen;
	std::uniform_real_distribution<float> randf(0.f, 1.f);

	glm::vec3 cam_pos = glm::vec3(glm::inverse(view_mat) * glm::vec4(0.f, 0.f, 0.f, 1.f));
	glm::mat4 inv_proj_view = glm::inverse(proj_mat * view_mat);
	
	#pragma omp parallel for
	for(index_t i = 0; i < width; i++)
	for(index_t j = 0; j < height; j++)
	{
		//row major
		glm::vec4 & color = img[j * windows_size.x + i];
		
			glm::vec2 screen = glm::vec2(	(float(i) + randf(gen)) / width, 
											(float(j) + randf(gen)) / height
											);

			glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
			glm::vec4 q = inv_proj_view * view;
			glm::vec3 p = glm::vec3(q * (1.f / q.w));

			ray_hit r(cam_pos, glm::normalize(p - cam_pos));
		
			color *= float(n_samples);

			if(intersect(r)) color += Li(r, light);

			color /= float(n_samples + 1);
	}

	n_samples++;
}

float * embree::raycaster(	const glm::uvec2 & windows_size,
							const glm::mat4 & view_mat,
							const glm::mat4 & proj_mat,
							const index_t & samples	)
{
	float * frame = new float[windows_size.x * windows_size.y];
	
	std::default_random_engine gen;
	std::uniform_real_distribution<float> randf(0.f, 1.f);

	glm::vec3 cam_pos = glm::vec3(glm::inverse(view_mat) * glm::vec4(0.f, 0.f, 0.f, 1.f));
	glm::mat4 inv_proj_view = glm::inverse(proj_mat * view_mat);

	#pragma omp parallel for
	for(index_t i = 0; i < windows_size.x; i++)
	for(index_t j = 0; j < windows_size.y; j++)
	{
		//row major
		float & color = frame[(windows_size.y - j - 1) * windows_size.x + i] = 0;
		
		for(index_t s = 0; s < samples; s++)
		{
			glm::vec2 screen = glm::vec2(	(float(i) + randf(gen)) / windows_size.x, 
											(float(j) + randf(gen)) / windows_size.y
											);

			glm::vec4 view = glm::vec4(screen.x * 2.f - 1.f, screen.y * 2.f - 1.f, 1.f, 1.f);
			glm::vec4 q = inv_proj_view * view;
			glm::vec3 p = glm::vec3(q * (1.f / q.w));

			ray_hit r(cam_pos, glm::normalize(p - cam_pos));
		
			if(intersect(r)) color += r.ray.tfar;
		}

		color /= samples;
	}

	return frame;
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

#endif // GPROSHAN_EMBREE

