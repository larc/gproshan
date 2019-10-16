#include "mdict/d_mesh_denoising.h"



// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_mesh_denoising(string file)
{
	che * mesh = new  che_off(file.c_str());

	size_t n = 4; 
	size_t m = 16;
	size_t M = 0;
	distance_t f = 1;
	bool learn = false;
	distance_t error;
//	gproshan_input(n m M f learn);
//	cin >> n >> m >> M >> f >> learn;
    basis * phi = new basis_dct(n);

	ofstream os("../tmp/test_mesh.txt");
	for(;f<2; f+=0.2)
	{
		denoising dict(mesh, phi, m, M, f, learn,0); 	
		error = dict.execute();
		os<< f << "\t"<<error<<endl;
	}
	os.close();
	delete phi;
}


} // namespace gproshan::mdict

