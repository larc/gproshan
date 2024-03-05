#include <gproshan/raytracing/embree.h>


#include <gproshan/util.h>
#include <gproshan/pointcloud/knn.h>

#include <random>
#include <cstring>


// geometry processing and shape analysis framework
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
	ray.mask = -1;
	ray.flags = 0;

	hit.geomID = RTC_INVALID_GEOMETRY_ID;
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
	return normalize(vertex{hit.Ng_x, hit.Ng_y, hit.Ng_z});
}

vertex embree::ray_hit::pos() const
{
	return org() + ray.tfar * dir();
}

void embree_error(void *, RTCError, const char * str)
{
	fprintf(stderr, "EMBREE ERROR: %s\n", str);
}


embree::embree()
{
	rtc_device = rtcNewDevice(NULL);
	rtc_scene = rtcNewScene(rtc_device);

	rtcSetSceneFlags(rtc_scene, RTC_SCENE_FLAG_COMPACT);

	rtcSetDeviceErrorFunction(rtc_device, embree_error, NULL);
}

embree::embree(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats, const pc_opts & pc): embree()
{
	build_bvh(meshes, model_mats, pc);
}

embree::~embree()
{
	for(CHE * m: g_meshes)
		delete m;

	rtcReleaseScene(rtc_scene);
	rtcReleaseDevice(rtc_device);
}

index_t embree::closest_vertex(const vertex & org, const vertex & dir) const
{
	ray_hit r(org, dir);
	if(!intersect(r)) return NIL;

	return closest_hit_vertex(*g_meshes[r.hit.geomID], r.hit);
}

eval_hit embree::intersect(const vertex & org, const vertex & dir, const bool flat) const
{
	ray_hit r(org, dir);
	if(!intersect(r)) return {};

	eval_hit hit(*g_meshes[r.hit.geomID], r.hit.primID, r.hit.u, r.hit.v, sc);
	hit.dist = r.ray.tfar;
	hit.position = r.pos();
	hit.normal = flat ? r.normal() : hit.normal;

	return hit;
}

void embree::build_bvh(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats, const pc_opts & pc)
{
	g_meshes.resize(size(meshes));
	for(index_t i = 0; i < size(meshes); ++i)
	{
		g_meshes[i] = new CHE(meshes[i]);

		if(!meshes[i]->n_trigs || pc.enable)
			g_meshes[i]->n_trigs = 0;

		[[maybe_unused]]
		const index_t geomID = g_meshes[i]->n_trigs || meshes[i]->is_scene() ?
											add_mesh(meshes[i], model_mats[i]) :
											add_pointcloud(meshes[i], model_mats[i], pc);

		gproshan_debug_var(i == geomID);
	}

	rtcCommitScene(rtc_scene);
}

index_t embree::add_sphere(const vec4 & xyzr)
{
	RTCGeometry geom = rtcNewGeometry(rtc_device, RTC_GEOMETRY_TYPE_SPHERE_POINT);

	vec4 * pxyzr = (vec4 *) rtcSetNewGeometryBuffer(	geom,
														RTC_BUFFER_TYPE_VERTEX, 0,
														RTC_FORMAT_FLOAT4, 4 * sizeof(float), 1
														);
	*pxyzr = xyzr;

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(rtc_scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

index_t embree::add_mesh(const che * mesh, const mat4 & model_mat)
{
	RTCGeometry geom = rtcNewGeometry(rtc_device, RTC_GEOMETRY_TYPE_TRIANGLE);

	vertex * vertices = (vertex *) rtcSetNewGeometryBuffer(	geom,
															RTC_BUFFER_TYPE_VERTEX, 0,
															RTC_FORMAT_FLOAT3, 3 * sizeof(float),
															mesh->n_vertices
															);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
		vertices[i] = model_mat * (mesh->point(i), 1);

	index_t * tri_idxs = (index_t *) rtcSetNewGeometryBuffer(	geom,
																RTC_BUFFER_TYPE_INDEX, 0,
																RTC_FORMAT_UINT3, 3 * sizeof(index_t),
																mesh->is_scene() ? mesh->n_vertices / 3 : mesh->n_trigs
																);


	if(mesh->is_scene())
	{
		#pragma omp parallel for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
			tri_idxs[i] = i;
	}
	else
	{
		memcpy(tri_idxs, mesh->trigs_ptr(), mesh->n_half_edges * sizeof(index_t));
	}

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(rtc_scene, geom);
	rtcReleaseGeometry(geom);

	if(mesh->is_scene())
	{
		g_meshes[geom_id]->VT = tri_idxs;
		g_meshes[geom_id]->n_trigs = mesh->n_vertices / 3;
		g_meshes[geom_id]->n_half_edges = mesh->n_vertices;

		scene * psc = (scene *) mesh;
		sc.materials = psc->materials.data();
		sc.textures = psc->textures.data();
		sc.trig_mat = psc->trig_mat;
		sc.texcoords = psc->texcoords;
	}

	return geom_id;
}

index_t embree::add_pointcloud(const che * mesh, const mat4 & model_mat, const pc_opts & pc)
{
	RTCGeometry geom = rtcNewGeometry(rtc_device, pc.normals ? RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT 
															: RTC_GEOMETRY_TYPE_DISC_POINT);

	vec4 * pxyzr = (vec4 *) rtcSetNewGeometryBuffer(	geom,
														RTC_BUFFER_TYPE_VERTEX, 0,
														RTC_FORMAT_FLOAT4,
														4 * sizeof(float),
														mesh->n_vertices
														);

	vertex * normals = !pc.normals	? nullptr
									: (vertex *) rtcSetNewGeometryBuffer(	geom,
																			RTC_BUFFER_TYPE_NORMAL, 0,
																			RTC_FORMAT_FLOAT3,
																			3 * sizeof(float),
																			mesh->n_vertices
																			);


	knn::k3tree * nn = pc.radius > 0 ? nullptr : new knn::k3tree(&mesh->point(0), mesh->n_vertices, pc.knn + 1);


	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
		pxyzr[i] = model_mat * (mesh->point(i), 1);

	#pragma omp parallel for
	for(index_t i = 0; i < mesh->n_vertices; ++i)
	{
		if(!nn)
		{
			pxyzr[i][3] = pc.radius;
			continue;
		}

		const real_t r = length(model_mat * (mesh->point(i) - mesh->point((*nn)(i, pc.knn)), 0));

		pxyzr[i][3] = pc.knn_area ? sqrt(r * r / pc.knn) : r;
	}


	delete nn;


	if(normals)
	{
		#pragma omp parallel for
		for(index_t i = 0; i < mesh->n_vertices; ++i)
			normals[i] = mesh->normal(i);
	}

	rtcCommitGeometry(geom);

	index_t geom_id = rtcAttachGeometry(rtc_scene, geom);
	rtcReleaseGeometry(geom);

	return geom_id;
}

bool embree::closesthit_radiance(	vertex & color,
									vertex & attenuation,
									vertex & position,
									vertex & ray_dir,
									float &,
									random<real_t> & rnd,
									const render_params & params,
									const bool flat
									) const
{
	ray_hit r(position, ray_dir);
	if(!intersect(r)) return false;

	const CHE & mesh = *g_meshes[r.hit.geomID];

	eval_hit hit(mesh, r.hit.primID, r.hit.u, r.hit.v, sc);
	hit.position = r.pos();
	hit.normal = flat ? r.normal() : hit.normal;

	color = eval_li(hit, params.ambient, params.lights, params.n_lights, params.cam_pos,
					[&](const vec3 & position, const vec3 & wi, const float light_dist) -> bool
					{
						ray_hit ro(position, wi, 1e-3f, light_dist - 1e-3f);
						return occluded(ro);
					});

	color *= attenuation;
	position = hit.position;

	if(!hit.scatter_mat(ray_dir, rnd))
		attenuation = 0;

	attenuation /= 2;

	return true;
}

float embree::intersect_depth(const vertex & org, const vertex & dir) const
{
	ray_hit r(org, dir);
	return intersect(r) ? r.ray.tfar : 0.f;
}

bool embree::intersect(ray_hit & r) const
{
	rtcIntersect1(rtc_scene, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}

bool embree::occluded(ray_hit & r) const
{
	rtcIntersect1(rtc_scene, &r);
	return r.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}


} // namespace gproshan

