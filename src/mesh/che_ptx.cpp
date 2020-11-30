#include "mesh/che_ptx.h"

#include <fstream>
#include <sstream>
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
	ifstream is(file);

	assert(is.good());

	size_t n_rows, n_cols;
	vertex T[4], R[3], tr;
	real_t s;

	is >> n_rows >> n_cols;

	is >> T[0] >> T[1] >> T[2] >> T[3];
	is >> R[0] >> s;
	is >> R[1] >> s;
	is >> R[2] >> s;
	is >> tr >> s;

	real_t * points = new real_t[7 * n_rows * n_cols];

	index_t p = 0;
	while(is >> points[p++]);

	// p is offset for [x y z a] or [x y z a r g b] point data
	p = n_rows * n_cols * 4 == p ? 4 : 7;

	init(n_rows * n_cols, 2 * (n_rows - 1) * (n_cols - 1));

	#pragma omp parallel for
	for(index_t v = 0; v < n_vertices_; v++)
	{
		const index_t & i = v * p;

		GT[v] = { points[i], points[i + 1], points[i + 2] };

		if(p == 4)
			VC[v] = { points[i + 3], points[i + 3], points[i + 3] };
		else
			VC[v] = { points[i + 4], points[i + 5], points[i + 6] };
	}

	delete [] points;

	index_t he = 0;
	auto add_trig = [&](const index_t & i, const index_t & j, const index_t & k)
	{
		if(GT[i].is_zero() || GT[j].is_zero() || GT[k].is_zero())
			return;
		
		VT[he++] = i;
		VT[he++] = j;
		VT[he++] = k;
		
		if(pdetriq(trig(he - 1)) < 0.1)
			he -= 3;
	};

	for(index_t r = 0; r < n_rows - 1; r++)
	for(index_t c = 0; c < n_cols - 1; c++)
	{
		add_trig((c    ) + (r    ) * n_cols,
				 (c    ) + (r + 1) * n_cols,
				 (c + 1) + (r    ) * n_cols);
		
		add_trig((c + 1) + (r + 1) * n_cols,
				 (c + 1) + (r    ) * n_cols,
				 (c    ) + (r + 1) * n_cols);
	}

	n_half_edges_ = he;
	n_faces_ = he / che::mtrig;

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

