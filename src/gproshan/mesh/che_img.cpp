#include <gproshan/mesh/che_img.h>

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <thread>

#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


che_img::che_img(const std::string & file)
{
	init(file);
}

void che_img::read_file(const std::string & file)
{
	CImg<unsigned char> img(file.c_str());

	alloc(img.height() * img.width(), 2 * (img.height() - 1) * (img.width() - 1));

	index_t v = 0, he = 0;
	for(int i = 0; i < img.width(); ++i)
	for(int j = 0; j < img.height(); ++j)
	{
		if(i && j)
		{
			VT[he++] = i * img.height() + j - 1;
			VT[he++] = (i - 1) * img.height() + j;
			VT[he++] = v;

			VT[he++] = (i - 1) * img.height() + j - 1;
			VT[he++] = (i - 1) * img.height() + j;
			VT[he++] = i * img.height() + j - 1;
		}

		const int c = img.width() - i - 1;
		const int r = img.height() - j - 1;

		GT[v] = {float(i), float(j), float(img.spectrum() == 1 ? - img(c, r) : 0)};

		if(img.spectrum() == 3)
			VC[v] = {img(c, r, 0), img(c, r, 1), img(c, r, 2)};

		++v;
	}

//	std::thread([](const CImg<unsigned char> & img) { img.display(); }, img).detach();
}


} // namespace gproshan

