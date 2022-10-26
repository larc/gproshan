#include <gproshan/raytracing/splat.h>

#include <gproshan/raytracing/splat_utils.h>


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat::splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats)
{
	for(che * m: meshes)
		add_splats_mesh(m);

	
}

splat::~splat()
{
	for(che * m: pointclouds)
		delete m;

	for(splats_data * spc: splats_pcs)
	{
		delete [] spc->pc;
		delete [] spc->morton_codes;
		delete [] spc->idx_splats;
	}
}

void splat::add_splats_mesh(const che * mesh)
{
	std::vector<index_t> vertices;
	
	che * pc = nullptr;
	splats_data * spc = new splats_data;

	

	pointclouds.push_back(pc);
	splats_pcs.push_back(spc);
}

void splat::build_splats_ch(splats_data * s, const che * pc, const std::vector<index_t> & vertices)
{

}


} // namespace gproshan

