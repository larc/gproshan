#ifndef SPLAT_H
#define SPLAT_H

#include <gproshan/raytracing/raytracing.h>
#include <gproshan/raytracing/splat_utils.h>

// geometry processing and shape analysis framework
namespace gproshan::rt {


class splat
{
	public:
		std::vector<che *> pointclouds;

	protected:
		std::vector<splats_data *> splats_pcs;

	public:
		splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats);
		virtual ~splat();

	protected:
		void add_splats_mesh(che * mesh, const mat4 & model_mat);
		void build_splats_ch(splats_data * s, const che * mesh, const std::vector<index_t> & vertices);
};


} // namespace gproshan

#endif // SPLAT_H

