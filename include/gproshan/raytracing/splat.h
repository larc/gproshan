#ifndef SPLAT_H
#define SPLAT_H

#include <gproshan/raytracing/splat_utils.h>
#include <gproshan/pointcloud/knn.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat
{
	public:
		std::vector<che *> pointclouds;
		static size_t k_nn;
		static real_t t_normal;
		static real_t d_overlap;

	private:
		double time = 0;
		double time_knn = 0;
		double time_segmentation = 0;
		double time_subdivision = 0;
		double time_initsplats = 0;

	protected:
		std::vector<splats_data> splats_pcs;

	public:
		splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats);
		virtual ~splat();

	private:
		void add_splats(che * pc, const mat4 & model_mat);

		std::vector<index_t> planar_segmentation(	const che * pc,
													std::vector<index_t> & vertices,
													const knn::k3tree & k3tree
													);

		std::vector<index_t> voronoi_subdivision(	std::vector<index_t> & voronoi_set,
													const vertex * points,
													const std::vector<index_t> & vertices,
													const index_t & seg_begin,
													const index_t & seg_end
													);

		che * init_splats(	const che * mesh,
							const mat4 & model_mat,
							std::vector<index_t> & vertices,
							const std::vector<index_t> & idx_splats
							);

		void display_sets(che * mesh, const std::vector<index_t> & sets, const index_t * mapid = nullptr);

		void save_stats(const std::string & file) const;
};


} // namespace gproshan

#endif // SPLAT_H

