#include "d_image_denoising.h"

#include <CImg.h>

using namespace cimg_library;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

void test_image_denoising(string file)
{
	CImg<double> image(file.c_str());
	image.resize(128, 128);

	size_t p = 8;
	size_t rows = image.width() - p + 1;
	size_t cols = image.height() - p + 1;
	size_t n = p * p;
	size_t m = 256;
	size_t M = rows * cols;
	size_t L = 10;

	mat X(n, M);

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

	mat D(n, m);
	D.randu();

	double time = omp_get_wtime();
	
	KSVD(D, X, L);
	
	time = omp_get_wtime() - time;
	cout << "time KSVD: " << time << endl;
	
	mat alpha(m, M);

	#pragma omp parallel for
	for(index_t i = 0; i < M; i++)
	{
		vec a;
		vec x = X.col(i);
		OMP(a, x, D, L);
		alpha.col(i) = a;
	}

	mat Y = D * alpha;
	
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
		int dx = p, dy = p;
		if(x < p && x < dx) dx = x + 1; 
		if(y < p && y < dy) dy = y + 1; 
		if((rows - x) < p && (rows - x) < dx) dx = rows - x; 
		if((cols - y) < p && (cols - y) < dy) dy = cols - y;
		
		image_out(x, y) /= dx * dy;
	}

	CImg<double> diff = abs(image - image_out);
	(image, image_out, diff).display();
}

} // mdict

