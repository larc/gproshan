#include "che_img.h"

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>
#include <thread>

#include <CImg.h>

using namespace cimg_library;

che_img::che_img(const size_t & n_v, const size_t & n_f)
{
	init(n_v, n_f);
}

che_img::che_img(const string & file)
{
	debug(file)
	init(file);
}

che_img::~che_img()
{

}

void che_img::read_file(const string & file)
{
	CImg<real_t> img(file.c_str());
	
	init(img.height() * img.width(), 2 * (img.height() - 1) * (img.width() - 1));
	
	index_t v = 0, he = 0;
	for(index_t i = 0; i < img.width(); i++)
	for(index_t j = 0; j < img.height(); j++)
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
		
		GT[v++] = vertex(i, img.height() - j - 1, -img(i, j));
	}

	thread([](CImg<real_t> img) { img.display(); }, img).detach();
}

void che_img::write_file(const string & file) const
{
	ofstream os(file);

	os << "OFF" << endl;
	os << n_vertices_ << " " << n_faces_ << " 0" << endl;

	for(size_t v = 0; v < n_vertices_; v++)
		os << GT[v] << endl;

	for(index_t he = 0; he < n_half_edges_; he++)
	{
		if(!(he % che::P)) os << che::P;
		os << " " << VT[he];
		if(he % che::P == che::P - 1) os << endl;
	}

	os.close();
}
