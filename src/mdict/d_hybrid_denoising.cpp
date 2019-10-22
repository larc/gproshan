#include "d_hybrid_denoising.h"

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_hybrid_denoising(const string & file)
{
	size_t N = 128;

	CImg<real_t> image(file.c_str());
	image.resize(N, N);
	image.save("../tmp/image_128.jpg");
	image = image.get_normalize(0, 1);

	size_t p = 8;							// square side of each patche
	size_t rows = image.width();
	size_t cols = image.height();
	size_t n = p * p;						// size of each patche
	size_t n_basis = 4;
	size_t m = n_basis * n_basis;							// number of atoms
	size_t M = rows * cols;					// number of patches
	size_t L = 10;							// sparsity OMP norm L_0
	size_t K = 10;							// KSVD iterations


	a_mat X(n, M);
	che * mesh = new che_img("../tmp/image_128.jpg");
	che_off::write_file(mesh,"../tmp/image_128");
	std::vector<patch> patches(M);				///< vector of patches.
	std::vector<vpatches_t> patches_map(M);		///< invert index vertex to patches.

	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t i = x + y * rows;

		for(index_t b = y; b < cols && b < y + p; b++)
		for(index_t a = x; a < rows && a < x + p; a++)
			patches[i].vertices.push_back(a + b * rows);
	}
	
	a_mat A;
	a_mat alpha;
	basis * phi_basis = new basis_dct(n_basis);
	A.eye(m, m);
	alpha.zeros(m, M);
	
	for(index_t s = 0; s < M; s++)
		patches[s].reset_xyz(mesh, patches_map, s, nullptr);

	//#pragma omp parallel for
	for(index_t s = 0; s < M; s++)
	{
		patch & p = patches[s];
		//p.T.eye(3, 3);
		//p.transform();
		p.phi.set_size(p.xyz.n_cols, n_basis*n_basis);
	//	gproshan_debug_var(p.xyz);
		phi_basis->discrete(p.phi, p.xyz);
	}

	OMP_all_patches_ksvt(alpha, A, patches, M, L);
	gproshan_debug_var(size(alpha));

	//Mesh reconstruction
	for(index_t p = 0; p < M; p++)
	{
		patch & rp = patches[p];

		if(rp.phi.n_rows)
		{
			a_vec x = rp.phi * A * alpha.col(p);
		//	gproshan_debug_var(x);          
			rp.xyz.row(2) = x.t();
		//	rp.itransform();
		}
	}
	a_mat V(3, mesh->n_vertices(), arma::fill::zeros);
	#pragma omp parallel for
	for(index_t v = 0; v < mesh->n_vertices(); v++)
	{
		if(patches_map[v].size())
			V.col(v) = mdict::simple_means_vertex(v, patches, patches_map);
	}
	vertex * new_vertices = (vertex *) V.memptr();
	mesh->set_vertices(new_vertices, mesh->n_vertices(), 0);


	CImg<double> image_out = image;
	image_out.fill(0);

	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t i = x + y * rows;

		image_out(x, y) = mesh->gt(i).z;
		gproshan_debug_var(mesh->gt(i).z);
	}
/*
	rows = image.width();
	cols = image.height();
	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t dx = p, dy = p;
		if(x < p && x < dx) dx = x + 1;
		if(y < p && y < dy) dy = y + 1;
		if((rows - x) < p && (rows - x) < dx) dx = rows - x;
		if((cols - y) < p && (cols - y) < dy) dy = cols - y;

		image_out(x, y) /= dx * dy;
	}
*/
//	CImg<double> diff = abs(image - image_out);
//	(image, image_out, diff).display();
//	(image_out).display();
//	image_out = image_out.get_normalize(0, 255);	
//	(image, image_out).display();
}


} // namespace gproshan::mdict

