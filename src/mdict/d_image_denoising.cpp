#include "d_image_denoising.h"

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_image_denoising(const string & file)
{
	CImg<real_t> image(file.c_str());
	image.resize(128, 128);
	image = image.get_normalize(0, 1);

	size_t p = 8;							// square side of each patche
	size_t rows = image.width() - p + 1;
	size_t cols = image.height() - p + 1;	
	size_t n = p * p;						// size of each patche
	size_t m = 256;							// number of atoms
	size_t M = rows * cols;					// number of patches
	size_t L = 10;							// sparsity OMP norm L_0
	size_t K = 10;							// KSVD iterations

	a_mat X(n, M);

	for(index_t x = 0; x < rows; x++)
	for(index_t y = 0; y < cols; y++)
	{
		index_t i = x + y * rows;
		index_t k = 0;

		for(index_t b = y; b < y + p; b++)
		for(index_t a = x; a < x + p; a++)
		{
			X(k, i) = image(a, b);
			k++;
		}
	}
	
	a_mat D(n, m, arma::fill::randu);
	D = normalise(D);
	
	CImgList<real_t> imlist;
	for(index_t i = 0; i < p; i++)
		imlist.push_back(CImg<real_t>(D.colptr(i), p, p, 1, 1, true));
	imlist.display();

	gproshan_log(KSVD);

	double time;

	TIC(time)
	KSVD(D, X, L, K);
	TOC(time)
	
	gproshan_log_var(time);
	
	imlist.clear();
	for(index_t i = 0; i < p; i++)
		imlist.push_back(CImg<real_t>(D.colptr(i), p, p, 1, 1, true));
	imlist.display();

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
	(image, image_out, diff).display();
}


} // namespace gproshan::mdict

