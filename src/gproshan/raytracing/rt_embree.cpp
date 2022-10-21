#include <gproshan/raytracing/rt_embree.h>


#include <gproshan/util.h>

#include <random>
#include <cstring>


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


embree::ray_hit::ray_hit(const vertex & p_org, const vertex & v_dir, float near, float far)
{
	ray.org_x = p_org.x();
	ray.org_y = p_org.y();
	ray.org_z = p_org.z();
	ray.tnear = near;

	ray.dir_x = v_dir.x();
	ray.dir_y = v_dir.y();
	ray.dir_z = v_dir.z();

	ray.time = 0.0f;

	ray.tfar = far;
	ray.mask = 0;
	ray.flags = 0;

	//hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
}

vertex embree::ray_hit::org() const
{
	return {ray.org_x, ray.org_y, ray.org_z};
}

vertex embree::ray_hit::dir() const
{
	return {ray.dir_x, ray.dir_y, ray.dir_z};
}

vertex embree::ray_hit::normal() const
{
	return normalize(vec3{hit.Ng_x, hit.Ng_y, hit.Ng_z});
}

vertex embree::ray_hit::position() const
{
	return org() + ray.tfar * dir();
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

	return mesh->VT[he];
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

embree::embree(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats, const bool & pointcloud, const float & pcr): embree()
{
	pc_radius = pcr;
	build_bvh(meshes, model_mats, pointcloud);
}

embree::~embree()
{
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);
}

index_t embree::closest_vertex(const vertex & org, const vertex & dir)
{
	ray_hit r(org, dir);
	return intersect(r) ? r.closest_vertex(geomID_mesh[r.hit.geomID]) : NIL;
}

hit embree::intersect(const vertex & org, const vertex & dir)
{
	ray_hit r(org, dir);
	if(intersect(r))
	{
		/*
		const rt_mesh & mesh = geomID_mesh[r.hit.geomID];
		const vertex & color = r.color(mesh);
		const vertex & normal = r.normal(mesh);

		return	{	r.closest_vertex(mesh),
					r.ray.tfar,
					{color.x(), color.y(), color.z()},
					{normal.x(), normal.y(), normal.z()}
					};*/
	}

	return hit();
}

void embree::build_bvh(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats, const bool & pointcloud)
{
	for(index_t i = 0; i < meshes.size(); ++i)
	{
		CHE * mesh = new CHE(meshes[i]);
		const mat4 & model_mat = model_mats[i];

		if(!mesh->n_faces || pointcloud)
			geomID_mesh[add_pointcloud(meshes[i], model_mat)] = {mesh, true};
		else
			geomID_mesh[add_mesh(meshes[i], model_mat)] = {mesh, false};
	}

	rtcCommitScene(scene);
}

index_t embree::add_sphere(const vec4 & xyzr)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

	vec4 * pxyzr = (vec4 *) rtcSetNewGeometryBuffer(	geom,
														RTC_BUFFER_TYPE_VERTEX, 0,
														RTC_FORMAT_FLOAT4, 4 * sizeof(float), 1
														);
	*pxyzr = xyzr;

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

index_t embree::add_mesh(const che * mesh, const mat4 & model_mat)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

	vertex * vertices = (vertex *) rtcSetNewGeometryBuffer(	geom,
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
		vertices[i] = model_mat * vec4(mesh->point(i), 1);

	memcpy(tri_idxs, &mesh->halfedge(0), mesh->n_half_edges * sizeof(index_t));

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

index_t embree::add_pointcloud(const che * mesh, const mat4 & model_mat)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

	vec4 * pxyzr = (vec4 *) rtcSetNewGeometryBuffer(	geom,
														RTC_BUFFER_TYPE_VERTEX, 0,
														RTC_FORMAT_FLOAT4,
														4 * sizeof(float),
														mesh->n_vertices
														);

	vertex * normal = (vertex *) rtcSetNewGeometryBuffer(	geom,
															RTC_BUFFER_TYPE_NORMAL, 0,
															RTC_FORMAT_FLOAT3,
															3 * sizeof(float),
															mesh->n_vertices
															);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		pxyzr[i] = model_mat * vec4(mesh->point(i), 1);
		pxyzr[i][3] = pc_radius;
		normal[i] = mesh->normal(i);
	}

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

vec3 embree::closesthit_radiance(const vertex & org, const vertex & dir, const vertex * lights, const int & n_lights, const bool & flat)
{
	ray_hit r(org, dir);
	if(!intersect(r)) return {};

	auto shadow = [&](const vec3 & position, const vec3 & wi, const float & light_dist) -> bool
					{
						ray_hit ro(position, wi, 1e-3f, light_dist - 1e-3f);
						return occluded(ro);
					};

	eval_hit<float, decltype(shadow)> hit(*geomID_mesh[r.hit.geomID].mesh, r.hit.primID, r.hit.u, r.hit.v);
	hit.position = r.position();
	hit.normal = flat ? r.normal() : hit.normal;

	return hit.eval_li(lights,  n_lights, shadow);
}

float embree::intersect_depth(const vertex & org, const vertex & dir)
{
	ray_hit r(org, dir);
	return intersect(r) ? r.ray.tfar : 0.f;
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


} // namespace gproshan

