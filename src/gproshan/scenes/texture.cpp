#include <gproshan/scenes/texture.h>


#include <CImg.h>

using namespace cimg_library;


// geometry processing and shape analysis framework
namespace gproshan {


texture::texture(const std::string & file)
{
	try
	{
		CImg<unsigned char> img(file.c_str());
		img.mirror('y');

		width = img.width();
		height = img.height();
		spectrum = img.spectrum();
		data = new unsigned char[width * height * spectrum];

		img.permute_axes("cxyz");
		memcpy(data, img.data(), width * height * spectrum);
	}
	catch(CImgException & e)
	{
		delete [] data;
		data = nullptr;
//		gproshan_error_var(e.what());
	}
}


} // namespace gproshan

