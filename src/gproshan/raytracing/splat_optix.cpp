#include <gproshan/raytracing/splat_optix.h>


#ifdef GPROSHAN_OPTIX


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat_optix::splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats), optix("/src/splat_optix.ptx")
{
	optix_params.traversable = build_as(pointclouds, {mat4::identity()});
	build_sbt();

	gproshan_error_var(splats_pcs.size());

	d_splats_pcs.resize(splats_pcs.size());

	for(index_t i = 0; i < splats_pcs.size(); ++i)
	{
		const che & p = *pointclouds[i];
		const splats_data & h = splats_pcs[i];
		splats_data & d = d_splats_pcs[i];

	gproshan_error_var(p.n_vertices);	

		d.n_splats = h.n_splats;
	gproshan_error(SO);	
		cudaMalloc(&d.morton_codes, sizeof(unsigned int) * p.n_vertices);
	gproshan_error(SO);	
		cudaMalloc(&d.primID_splat, sizeof(index_t) * p.n_trigs);
	gproshan_error(SO);	
		cudaMalloc(&d.splats, sizeof(splat_t<real_t>) * d.n_splats);
	gproshan_error(SO);	

		cudaMemcpy(d.morton_codes, h.morton_codes, sizeof(unsigned int) * p.n_vertices, cudaMemcpyHostToDevice);
	gproshan_error(SO);	
		cudaMemcpy(d.primID_splat, h.primID_splat, sizeof(index_t) * p.n_trigs, cudaMemcpyHostToDevice);
	gproshan_error(SO);	
		cudaMemcpy(d.splats, h.splats, sizeof(splat_t<real_t>) * d.n_splats, cudaMemcpyHostToDevice);

		gproshan_error_var(d.n_splats);
		gproshan_error_var(h.n_splats);
	}

	cudaMalloc(&dd_splats_pcs, sizeof(splats_data) * splats_pcs.size());
	cudaMemcpy(dd_splats_pcs, d_splats_pcs.data(), sizeof(splats_data) * splats_pcs.size(), cudaMemcpyHostToDevice);

	optix_params.other = dd_splats_pcs;
}


splat_optix::~splat_optix()
{
	for(splats_data & sd: d_splats_pcs)
	{
		cudaFree(sd.morton_codes);
		cudaFree(sd.primID_splat);
		cudaFree(sd.splats);

		sd.morton_codes = nullptr;
		sd.primID_splat = nullptr;
		sd.splats = nullptr;
	}

	cudaFree(dd_splats_pcs);
}


} // namespace gproshan

#endif // GPROSHAN_OPTIX

