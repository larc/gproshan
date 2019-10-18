#include "mdict/d_mesh_denoising.h"

#include "mdict/basis_dct.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


void test_mesh_denoising(const string & file)
{
	che * mesh = new che_off(file.c_str());

	size_t n = 4; 
	size_t m = 16;
	size_t M = 0;
	distance_t f = 1.2;
	bool learn = false;
	distance_t error;
	//dictionary::L = 20;
	basis * phi = new basis_dct(n);

	//denoising dict(mesh, phi, m, M, f, learn,0);
	//dict.execute(); 

	ofstream os("../tmp/test.txt");

	for(; f<1.4; f+=0.1)
	{
		os<< f ;
		for(size_t i = 10; i<26; i+=5)
		{
			dictionary::L = i;
			mesh = new che_off(file.c_str());
			denoising dict(mesh, phi, m, M, f, learn,0); 	
			error = dict.execute();
			os<< "\t"<<error;
			//os<< "\t"<<0.001;
		}
		os<<endl;
	}

	os.close();
	
	system("gnuplot -persist ../tmp/test_mesh.gp &");
	
	delete phi;
}


} // namespace gproshan::mdict

