#include "d_mesh_denoising.h"



// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_mesh_denoising(string file)
{
	che_off mesh(file.c_str());

	size_t n = 4; 
	size_t m = 16;
	size_t M = 0;
	distance_t f = 1;
	bool learn = 0;

//	gproshan_input(n m M f learn);
//	cin >> n >> m >> M >> f >> learn;
    basis * phi = new basis_dct(n);
	denoising dict(mesh, phi, m, M, f, learn);
	dict.execute();
	
	delete phi;
}


} // namespace gproshan::mdict

