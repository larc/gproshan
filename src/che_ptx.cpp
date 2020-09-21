#include "che_ptx.h"

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

	char line[256];
	is.getline(line, sizeof(line));

	for(index_t i = 0; i < n_vertices_; i++)
	{
		is.getline(line, sizeof(line));
		stringstream ss(line);
		
		ss >> GT[i];
	}
	
	is.close();
}

void che_ptx::write_file(const che * mesh, const string & file)
{
	ofstream os(file + ".off");


	os.close();
}


} // namespace gproshan

