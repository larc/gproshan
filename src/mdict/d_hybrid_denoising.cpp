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
	size_t rows = image.width() - p + 1;
	size_t cols = image.height() - p + 1;
	size_t col = image.height();	
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
	std::vector<vpatches_t> patches_map;		///< invert index vertex to patches.
	patches_map.resize(M);
	gproshan_debug(patches);

	index_t s = 0;

	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t i = x + y * rows;

		for(index_t b = y; b < y + p; b++)
		for(index_t a = x; a < x + p; a++)
		{
			patches[s].vertices.push_back(a * col + b);

			//cout<< a*col+b<<" ";
		}
		
		s++;
	}
	
	a_mat A;
	a_mat alpha;
	basis * phi_basis = new basis_dct(n_basis);
	A.eye(m, m);
	alpha.zeros(m, M);
	

	for(index_t s = 0; s < M; s++)
	{
		patches[s].reset_xyz(mesh, patches_map, s, nullptr);
	}
	

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
	/*
	a_mat D(n, m, arma::fill::randu);
	D = normalise(D);
	
	CImg<real_t> imdict;
	for(index_t i = 0; i < 16; i++)
	{
		CImg<real_t> imrow;
		for(index_t j = 0; j < 16; j++)
			imrow.append(CImg<real_t>(D.colptr(i * 16 + j), p, p, 1, 1, true), 'x');

		imdict.append(imrow, 'y');
	}
	imdict.display();

	gproshan_log(KSVD);

	double time;

	TIC(time)
	KSVD(D, X, L, K);
	TOC(time)
	
	gproshan_log_var(time);
	
	CImg<real_t> imdictlearned;
	for(index_t i = 0; i < 16; i++)
	{
		CImg<real_t> imrow;
		for(index_t j = 0; j < 16; j++)
			imrow.append(CImg<real_t>(D.colptr(i * 16 + j), p, p, 1, 1, true), 'x');

		imdictlearned.append(imrow, 'y');
	}
	(imdict, imdictlearned).display();

	a_mat alpha(m, M);

	gproshan_log(OMP);

	TIC(time)
	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
		alpha.col(i) = OMP(X.col(i), D, L);
	TOC(time)
	
	gproshan_log_var(time);

	a_mat Y = D * alpha;

	CImg<double> image_out = image;
	image_out.fill(0);

	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t i = x + y * rows;
		index_t k = 0;

		for(index_t b = y; b < y + p; b++)
		for(index_t a = x; a < x + p; a++)
		{
			image_out(a, b) += Y(k, i);
			k++;
		}
	}

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

	CImg<double> diff = abs(image - image_out);
	(image, image_out, diff).display();*/
	(image).display();
}


} // namespace gproshan::mdict

