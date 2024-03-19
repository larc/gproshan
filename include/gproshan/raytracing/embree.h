#ifndef RT_EMBREE_H
#define RT_EMBREE_H

#include <gproshan/mesh/che.h>
#include <gproshan/scenes/scene.h>
#include <gproshan/raytracing/raytracing.h>

#include <cfloat>

#include <embree4/rtcore.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class embree : public raytracing
{
	public:
		enum knn_opt {	NONE,
						MAX,
						MEAN,
						MEDIAN,
						AREA,
						MEDIAN_PAIRS
						};
		struct pc_opts
		{
			bool enable		= false;
			bool normals	= false;
			float radius	= 0.01;
			knn_opt opt		= NONE;
			float scale		= 1;
			int knn			= 8;


			pc_opts() {};
		};

	protected:
		struct ray_hit : public RTCRayHit
		{
			ray_hit(const vertex & p_org = {0, 0, 0},
					const vertex & v_dir = {0, 0, 0},
					float near = 1e-5f,
					float far = 1e20f
					);

			vertex org() const;
			vertex dir() const;
			vertex pos() const;
			vertex normal() const;
		};


		RTCDevice rtc_device;
		RTCScene rtc_scene;

		std::vector<CHE *> g_meshes;
		scene_data sc;

	public:
		embree();
		embree(	const std::vector<che *> & meshes,
				const std::vector<mat4> & model_mats,
				const pc_opts & pc = pc_opts()
				);
		virtual ~embree();

		virtual index_t closest_vertex(const vertex & org, const vertex & dir) const;
		virtual eval_hit intersect(const vertex & org, const vertex & dir, const bool flat = true) const;

		vec4 * pc_data(const index_t geomID = 0);


	protected:
		void build_bvh(	const std::vector<che *> & meshes,
						const std::vector<mat4> & model_mats,
						const pc_opts & pc = pc_opts()
						);

		index_t add_sphere(const vec4 & xyzr);
		index_t add_mesh(const che * mesh, const mat4 & model_mat);

		virtual index_t add_pointcloud(const che * mesh, const mat4 & model_mat, const pc_opts & pc);

		virtual bool closesthit_radiance(	vertex & color,
											vertex & attenuation,
											vertex & position,
											vertex & ray_dir,
											real_t & dist,
											random<real_t> & rnd,
											const render_params & params,
											const bool flat
											) const;

		float intersect_depth(const vertex & org, const vertex & dir) const;

		bool intersect(ray_hit & r) const;
		bool occluded(ray_hit & r) const;
};


} // namespace gproshan

#endif // RT_EMBREE_H

