#ifndef SPLAT_H
#define SPLAT_H

#include <gproshan/raytracing/splat_utils.h>

// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat
{
	public:
		std::vector<che *> pointclouds;
		static int k;
		static real_t n_threshold;
		static real_t delta;

	protected:
		std::vector<splats_data *> splats_pcs;

	public:
		splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats);
		virtual ~splat();

	private:
		void add_splats(che * pc, const mat4 & model_mat);

		std::vector<index_t> planar_segmentation(che * pc, std::vector<index_t> & vertices);

		std::vector<index_t> voronoi_subdivision(	std::vector<index_t> & voronoi_set,
													const vertex * points,
													const std::vector<index_t> & vertices,
													const index_t & seg_begin,
													const index_t & seg_end
													);

		void init_splats(const che * mesh, const mat4 & model_mat, std::vector<index_t> & vertices, const std::vector<index_t> & idx_splats);

};


} // namespace gproshan

#endif // SPLAT_H

