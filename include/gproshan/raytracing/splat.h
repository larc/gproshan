#ifndef SPLAT_H
#define SPLAT_H

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

	private:
		void add_splats_mesh(che * mesh, const mat4 & model_mat);
};


} // namespace gproshan

#endif // SPLAT_H

