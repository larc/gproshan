#include <gproshan/scenes/scanner.h>

#include <cmath>
#include <thread>

#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-truncation"
	#include <CImg.h>
#pragma GCC diagnostic pop
#else
	#include <CImg.h>
#endif // __clang__

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


che * scanner_ptx(const rt::raytracing * rt, const size_t n_rows, const size_t n_cols, const vertex & cam_pos, const bool dist_error)
{
	che * mesh_ptx = new che(n_cols * n_rows);

	const real_t delta_phi = (2 * M_PI) / n_rows;
	const real_t delta_theta = M_PI / n_cols;

	#pragma omp parallel for
	for(index_t i = 0; i < n_rows; ++i)
	for(index_t j = 0; j < n_cols; ++j)
	{
		const index_t v = i * n_cols + j;

		const real_t phi = i * delta_phi;
		const real_t theta = j * delta_theta;
		const vertex & dir = {std::sin(theta) * std::sin(phi), std::cos(theta), std::sin(theta) * std::cos(phi)};

		const rt::eval_hit & h = rt->intersect(cam_pos, dir);

		mesh_ptx->point(v) = h.position;
		mesh_ptx->normal(v) = h.normal;
		mesh_ptx->heatmap(v) = dist_error ? h.dist : h.heatmap;
		mesh_ptx->rgb(v) = h.Kd;
	}

	real_t max_dist = 0;

	#pragma omp parallel for reduction(max: max_dist)
	for(index_t v = 0; v < mesh_ptx->n_vertices; ++v)
		max_dist = std::max(max_dist, mesh_ptx->heatmap(v));

	#pragma omp parallel for
	for(index_t v = 0; v < mesh_ptx->n_vertices; ++v)
		mesh_ptx->heatmap(v) /= max_dist;

	return mesh_ptx;
}

che * scanner_ptx_jpg(const rt::raytracing * rt, const size_t n_rows, const size_t n_cols, const vertex & cam_pos, const std::string & file_jpg)
{
	che * mesh_ptx = scanner_ptx(rt, n_rows, n_cols, cam_pos);

	CImg<unsigned char> img((unsigned char *) &mesh_ptx->rgb(0), 3, n_cols, n_rows);
	img.permute_axes("zycx");

	img.save((file_jpg + "_scan.jpg").c_str());

	std::thread([](const CImg<unsigned char> & img) { img.display(); }, img).detach();

	return mesh_ptx;
}


} // namespace gproshan

