#include <gproshan/scenes/scanner.h>

#include <thread>
#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


che * scanner_ptx(const raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam_pos)
{
	std::vector<vertex> ptx;

	return new che(ptx.data(), ptx.size(), nullptr, 0);
}

che * scanner_ptx(const che * mesh, raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam_pos, const std::string & file_jpg)
{
	che * mesh_ptx = new che(n_cols * n_rows);

	const real_t delta_phi = (2 * M_PI) / n_rows;
	const real_t delta_theta = M_PI / n_cols;

	#pragma omp parallel for
	for(index_t i = 0; i < n_rows; ++i)
	for(index_t j = 0; j < n_cols; ++j)
	{
		const index_t & v = i * n_cols + j;

		const real_t & phi = i * delta_phi;
		const real_t & theta = j * delta_theta;
		const vertex & dir = {	sin(theta) * cos(phi),
									sin(theta) * sin(phi),
									cos(theta)
									};

		const hit & h = rt->intersect(cam_pos, dir);

		if(h.idx != NIL)
		{
			const vertex & pos = cam_pos + dir * h.dist;
			mesh_ptx->get_vertex(v) = {pos.x, pos.y, pos.z};
			mesh_ptx->normal(v) = h.normal;
			mesh_ptx->heatmap(v) = h.dist / M_SQRT2;
			mesh_ptx->rgb(v) = {	(unsigned char) (h.color.x * 255),
									(unsigned char) (h.color.y * 255),
									(unsigned char) (h.color.z * 255)
									};
		}
		else
		{
			mesh_ptx->rgb(v) = {0, 0, 0};
		}
	}

	CImg<unsigned char> img((unsigned char *) &mesh_ptx->rgb(0), 3, n_cols, n_rows);
	img.permute_axes("zycx");

	std::string img_filename = file_jpg + mesh->name() + ".jpg";
	img.save(img_filename.c_str());

	std::thread([](CImg<real_t> img) { img.display(); }, img).detach();

	return mesh_ptx;
}


} // namespace gproshan

