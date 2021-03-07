#include "geodesics/sampling.h"

#include "geodesics/geodesics_ptp.h"
#include "mesh/che_off.h"

#include <fstream>

using namespace std;


// geometry processing and shape analysis framework
namespace gproshan {


index_t ** sampling_shape(vector<index_t> & points, size_t *& sizes, vertex *& normals, che * mesh, size_t n_points, real_t radio)
{
	normals = new vertex[n_points];
	sizes = new size_t[n_points];
	index_t ** indexes = new index_t * [n_points];
	
	geodesics::params params;
	params.radio = radio;
	
	#pragma omp parallel for
	for(index_t i = 0; i < n_points; ++i)
	{
		const index_t & v = points[i];
		normals[i] = mesh->normal(v);

		geodesics fm(mesh, { v }, params);

		indexes[i] = new index_t[fm.n_sorted_index()];

		fm.copy_sorted_index(indexes[i], fm.n_sorted_index());
		sizes[i] = fm.n_sorted_index();
	}

	return indexes;
}

bool load_sampling(vector<index_t> & points, real_t & radio, che * mesh, size_t n)
{
	const string & filename = mesh->filename;

	string file = filename.substr(filename.find_last_of('/'), filename.size() - filename.find_last_of('/')) + "." + to_string(n);

	ifstream is(tmp_file_path(file));
	gproshan_log_var(tmp_file_path(file));

	if(is.good())
	{
		is >> radio;

		size_t n, p;
		is >> n;

		while(n--)
		{
			is >> p;
			points.push_back(p);
		}
	}
	else
	{
		if(!points.size())
			points.push_back(0);

#ifdef GPROSHAN_CUDA
		double time_fps;
		radio = farthest_point_sampling_ptp_gpu(mesh, points, time_fps, n);
		gproshan_debug_var(time_fps);
#else
		radio = 0; // IMPLEMENT: farthest_point_sampling_ptp_cpu(mesh, points, time_fps, n);
#endif // GPROSHAN_CUDA

		ofstream os(tmp_file_path(file));
		os << radio << endl;
		os << points.size() << endl;
		for(const index_t & i: points)
			os << i << endl;

		os.close();
	}

	is.close();

	return true;
}


} // namespace gproshan

