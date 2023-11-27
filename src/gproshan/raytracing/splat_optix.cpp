#include <gproshan/raytracing/splat_optix.h>


#ifdef GPROSHAN_OPTIX


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat_optix::splat_optix(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats): splat(meshes, model_mats), optix("/src/splat_optix.ptx")
{
	optix_params.traversable = build_as(pointclouds, {mat4::identity()});
	build_sbt();

	d_splats_pcs.resize(size(splats_pcs));

	for(index_t i = 0; i < size(splats_pcs); ++i)
	{
		const che & p = *pointclouds[i];
		const splats_data & h = splats_pcs[i];
		splats_data & d = d_splats_pcs[i];

		d.n_splats = h.n_splats;
		cudaMalloc(&d.morton_codes, sizeof(unsigned int) * p.n_vertices);
		cudaMalloc(&d.primID_splat, sizeof(index_t) * p.n_trigs);
		cudaMalloc(&d.splats, sizeof(splat_t<real_t>) * d.n_splats);

		cudaMemcpy(d.morton_codes, h.morton_codes, sizeof(unsigned int) * p.n_vertices, cudaMemcpyHostToDevice);
		cudaMemcpy(d.primID_splat, h.primID_splat, sizeof(index_t) * p.n_trigs, cudaMemcpyHostToDevice);
		cudaMemcpy(d.splats, h.splats, sizeof(splat_t<real_t>) * d.n_splats, cudaMemcpyHostToDevice);

		gproshan_error_var(d.n_splats);
		gproshan_error_var(h.n_splats);
	}

	cudaMalloc(&dd_splats_pcs, sizeof(splats_data) * size(splats_pcs));
	cudaMemcpy(dd_splats_pcs, d_splats_pcs.data(), sizeof(splats_data) * size(splats_pcs), cudaMemcpyHostToDevice);

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

