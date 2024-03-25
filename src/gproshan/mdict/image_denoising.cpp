#include <gproshan/mdict/image_denoising.h>

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
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_image_denoising(const std::string & file)
{
	CImg<float> image(file.c_str());
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

	arma::fmat X(n, M);

	for(index_t x = 0; x < rows; ++x)
	for(index_t y = 0; y < cols; ++y)
	{
		index_t i = x + y * rows;
		index_t k = 0;

		for(index_t b = y; b < y + p; ++b)
		for(index_t a = x; a < x + p; ++a)
			X(k++, i) = image(a, b);
	}

	arma::fmat D(n, m, arma::fill::randu);
	D = normalise(D);

	arma::fmat spD = D;

	CImg<float> imdict;
	for(index_t i = 0; i < 16; ++i)
	{
		CImg<float> imrow;
		for(index_t j = 0; j < 16; ++j)
			imrow.append(CImg<float>(D.colptr(i * 16 + j), p, p, 1, 1, true), 'x');

		imdict.append(imrow, 'y');
	}
	imdict.display();

	gproshan_log(KSVD);

	double time;

	TIC(time)
	KSVD(D, X, L, K);
	TOC(time)

	gproshan_log_var(time);

	TIC(time)
	sp_KSVD(spD, X, L, K);
	TOC(time)

	gproshan_log_var(time);

	gproshan_log_var(norm(D - spD));

	CImg<float> imdictlearned;
	for(index_t i = 0; i < 16; ++i)
	{
		CImg<float> imrow;
		for(index_t j = 0; j < 16; ++j)
			imrow.append(CImg<float>(D.colptr(i * 16 + j), p, p, 1, 1, true), 'x');

		imdictlearned.append(imrow, 'y');
	}
	(imdict, imdictlearned).display();


	gproshan_log(OMP);

	TIC(time)
	arma::fmat Y = D * OMP_all(X, D, L);
	TOC(time)

	gproshan_log_var(time);

	TIC(time)
	std::vector<locval_t> locval;
	arma::fmat spY = D * OMP_all(locval, X, D, L);
	TOC(time)

	gproshan_log_var(time);

	gproshan_log_var(norm(Y - spY));

	CImg<double> image_out = image;
	image_out.fill(0);

	for(index_t x = 0; x < rows; ++x)
	for(index_t y = 0; y < cols; ++y)
	{
		index_t i = x + y * rows;
		index_t k = 0;

		for(index_t b = y; b < y + p; ++b)
		for(index_t a = x; a < x + p; ++a)
			image_out(a, b) += Y(k++, i);
	}

	rows = image.width();
	cols = image.height();
	for(index_t x = 0; x < rows; ++x)
	for(index_t y = 0; y < cols; ++y)
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

