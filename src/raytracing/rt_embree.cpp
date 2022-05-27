#include <gproshan/raytracing/rt_embree.h>


#include <gproshan/util.h>

#include <random>
#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


embree::ray_hit::ray_hit(const glm::vec3 & p_org, const glm::vec3 & v_dir, float near, float far)
{
	ray.org_x = p_org.x;
	ray.org_y = p_org.y;
	ray.org_z = p_org.z;
	ray.tnear = near;

	ray.dir_x = v_dir.x;
	ray.dir_y = v_dir.y;
	ray.dir_z = v_dir.z;

	ray.time = 0.0f;

	ray.tfar = far;
	ray.mask = 0;
	ray.flags = 0;

	//hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
}

glm::vec3 embree::ray_hit::org() const
{
	return {ray.org_x, ray.org_y, ray.org_z};
}

glm::vec3 embree::ray_hit::dir() const
{
	return {ray.dir_x, ray.dir_y, ray.dir_z};
}

glm::vec3 embree::ray_hit::color(const rt_mesh & mesh) const
{
	const vertex & c = mesh.pointcloud ? mesh->color(hit.primID) :
										 mesh->shading_color(hit.primID, 1.0 - hit.u - hit.v, hit.u, hit.v);
	return glm_vec3(c);
}

glm::vec3 embree::ray_hit::normal(const rt_mesh & mesh, const bool & flat) const
{
	if(flat || mesh.pointcloud)
		return glm::normalize(glm::vec3(hit.Ng_x, hit.Ng_y, hit.Ng_z));

	const vertex & n = mesh->shading_normal(hit.primID, 1.0 - hit.u - hit.v, hit.u, hit.v);
	return glm::normalize(glm_vec3(n));
}

index_t embree::ray_hit::closest_vertex(const rt_mesh & mesh) const
{
	if(mesh.pointcloud) return hit.primID;

	index_t he = che::mtrig * hit.primID;
	float w = 1 - hit.u - hit.v;

	if(w < hit.u)
	{
		he = che::mtrig * hit.primID + 1;
		w = hit.u;
	}

	if(w < hit.v)
	{
		he = che::mtrig * hit.primID + 2;
		w = hit.v;
	}

	return mesh->vt(he);
}

glm::vec3 embree::ray_hit::position() const
{
	return org() + ray.tfar * dir();
}


void embree_error(void *, RTCError, const char * str)
{
	fprintf(stderr, "EMBREE ERROR: %s\n", str);
}


embree::embree()
{
	device = rtcNewDevice(NULL);
	scene = rtcNewScene(device);

	rtcSetSceneFlags(scene, RTC_SCENE_FLAG_COMPACT);

	rtcInitIntersectContext(&intersect_context);
	rtcSetDeviceErrorFunction(device, embree_error, NULL);
}

embree::embree(const std::vector<che *> & meshes, const std::vector<glm::mat4> & model_mats, const bool & pointcloud, const float & pcr): embree()
{
	pc_radius = pcr;
	build_bvh(meshes, model_mats, pointcloud);
}

embree::~embree()
{
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);
}

index_t embree::cast_ray(const glm::vec3 & org, const glm::vec3 & dir)
{
	ray_hit r(org, dir);
	return intersect(r) ? r.closest_vertex(geomID_mesh[r.hit.geomID]) : NIL;
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

void embree::build_bvh(const std::vector<che *> & meshes, const std::vector<glm::mat4> & model_mats, const bool & pointcloud)
{
	for(index_t i = 0; i < meshes.size(); ++i)
	{
		che * mesh = meshes[i];
		const glm::mat4 & model_mat = model_mats[i];

		if(mesh->is_pointcloud() || pointcloud)
			geomID_mesh[add_pointcloud(mesh)] = {mesh, true};
		else
			geomID_mesh[add_mesh(mesh, model_mat)] = {mesh, false};
	}

	rtcCommitScene(scene);
}

index_t embree::add_sphere(const glm::vec4 & xyzr)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

	glm::vec4 * pxyzr = (glm::vec4 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_VERTEX, 0,
																RTC_FORMAT_FLOAT4, 4 * sizeof(float), 1);
	*pxyzr = xyzr;

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

index_t embree::add_mesh(const che * mesh, const glm::mat4 & model_mat)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

	glm::vec3 * vertices = (glm::vec3 *) rtcSetNewGeometryBuffer(	geom,
																	RTC_BUFFER_TYPE_VERTEX, 0,
																	RTC_FORMAT_FLOAT3, 3 * sizeof(float),
																	mesh->n_vertices
																	);

	index_t * tri_idxs = (index_t *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_INDEX, 0,
																RTC_FORMAT_UINT3, 3 * sizeof(index_t),
																mesh->n_faces
																);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
		vertices[i] = glm::vec3(model_mat * glm::vec4(glm_vec3(mesh->gt(i)), 1));

	memcpy(tri_idxs, &mesh->vt(0), mesh->n_half_edges * sizeof(index_t));

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

index_t embree::add_pointcloud(const che * mesh)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

	glm::vec4 * pxyzr = (glm::vec4 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_VERTEX, 0,
																RTC_FORMAT_FLOAT4,
																4 * sizeof(float),
																mesh->n_vertices
																);

	glm::vec3 * normal = (glm::vec3 *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_NORMAL, 0,
																RTC_FORMAT_FLOAT3,
																3 * sizeof(float),
																mesh->n_vertices
																);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		pxyzr[i] = glm::vec4(glm_vec3(mesh->gt(i)), pc_radius);
		normal[i] = glm_vec3(mesh->normal(i));
	}

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

float embree::pointcloud_hit(glm::vec3 & position, glm::vec3 & normal, glm::vec3 & color, ray_hit r)
{
	float w, sum_w = 0;
	position = glm::vec3(0);
	normal = glm::vec3(0);
	color = glm::vec3(0);

	do
	{
		glm::vec4 * xyzr = (glm::vec4 *) rtcGetGeometryBufferData(rtcGetGeometry(scene, r.hit.geomID), RTC_BUFFER_TYPE_VERTEX, 0);

		sum_w += w = pc_radius - glm::length(r.position() - glm::vec3(xyzr[r.hit.primID]));
		position += w * r.position();
		normal += w * r.normal(geomID_mesh[r.hit.geomID]);
		color += w * r.color(geomID_mesh[r.hit.geomID]);

		r = ray_hit(r.position(), r.dir());
	}
	while(intersect(r) && r.ray.tfar < pc_radius);

	position /= sum_w;
	normal /= sum_w;
	color /= sum_w;

	return sum_w;
}

glm::vec4 embree::li(const glm::vec3 & light, const glm::vec3 & position, const glm::vec3 & normal, const glm::vec3 & color, const float & near)
{
	const glm::vec3 wi = normalize(light - position);
	const float dot_wi_normal = glm::dot(wi, normal);
	const glm::vec3 L = (dot_wi_normal < 0 ? -dot_wi_normal : dot_wi_normal) * color;

	ray_hit r(position, wi, near);
	return glm::vec4((occluded(r) ? 0.4f : 1.f) * L, 1);
}

glm::vec4 embree::li(ray_hit r, const glm::vec3 & light, const bool & flat)
{
	float total_tfar = 0;

	float near;
	glm::vec3 position, normal, color;

	glm::vec4 L(0);
//	while(total_tfar < 0.1)
	{
		total_tfar += r.ray.tfar;

		position = r.position();
		normal = r.normal(geomID_mesh[r.hit.geomID], flat);
		color = r.color(geomID_mesh[r.hit.geomID]);

		near = 1e-5f;
		if(geomID_mesh[r.hit.geomID].pointcloud)
			near += pointcloud_hit(position, normal, color, r);

		L += r.ray.tfar * li(light, position, normal, color, near);

		r = ray_hit(r.position(), r.dir());
//		if(!intersect(r))
//			break;
	}

	return L / total_tfar;
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

