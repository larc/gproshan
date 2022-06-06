#include <gproshan/mesh/che_img.h>

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <thread>

#include <CImg.h>


using namespace std;
using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


che_img::che_img(const string & file)
{
	init(file);
}

void che_img::read_file(const string & file)
{
	CImg<real_t> img(file.c_str());

	alloc(img.height() * img.width(), 2 * (img.height() - 1) * (img.width() - 1));

	index_t v = 0, he = 0;
	for(int i = 0; i < img.width(); ++i)
	for(int j = 0; j < img.height(); ++j)
	{
		if(i && j)
		{
			VT[he++] = v;
			VT[he++] = (i - 1) * img.height() + j;
			VT[he++] = i * img.height() + j - 1;

			VT[he++] = i * img.height() + j - 1;
			VT[he++] = (i - 1) * img.height() + j;
			VT[he++] = (i - 1) * img.height() + j - 1;
		}

		GT[v++] = {real_t(i), real_t(j), img(i, j)};
	}

	thread([](CImg<real_t> img) { img.display(); }, img).detach();
}


} // namespace gproshan

