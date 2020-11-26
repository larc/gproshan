#include "mesh/che_ptx.h"

#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


che_ptx::che_ptx(const string & file)
{
	init(file);
}

che_ptx::che_ptx(const che_ptx & mesh): che(mesh)
{
}

che_ptx::~che_ptx()
{
}

void che_ptx::read_file(const string & file)
{
	size_t n_rows, n_cols;
	vertex T[4], R[3], tr;
	real_t s;

	ifstream is(file);

	assert(is.good());

	is >> n_rows >> n_cols;
	init(n_rows * n_cols, 0);

	is >> T[0] >> T[1] >> T[2] >> T[3];
	is >> R[0] >> s;
	is >> R[1] >> s;
	is >> R[2] >> s;
	is >> tr >> s;

	VN = new vertex[n_vertices_];
	VC = new vertex[n_vertices_];
	
	char line[256];
	is.getline(line, sizeof(line));

	real_t alpha;	// color
	for(index_t i = 0; i < n_vertices_; i++)
	{
		is.getline(line, sizeof(line));
		stringstream ss(line);
		
		ss >> GT[i] >> alpha >> VC[i];

		VN[i] = tr - GT[i];
		VN[i].unit();
	}
	
	#pragma omp parallel for
	for(index_t i = 0; i < n_vertices_; i++)
		VC[i] /= 255;
	
	is.close();
}

void che_ptx::write_file(const che * mesh, const string & file)
{
	ofstream os(file + ".off");


	os.close();
}


} // namespace gproshan

