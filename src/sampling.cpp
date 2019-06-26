#include "sampling.h"
#include "geodesics_ptp.h"
#include "che_off.h"

#include <fstream>

index_t ** sampling_shape(vector<index_t> & points, size_t *& sizes, vertex *& normals, che * shape, size_t n_points, distance_t radio)
{
	normals = new vertex[n_points];
	sizes = new size_t[n_points];
	index_t ** indexes = new index_t * [n_points];

	index_t v;

	#pragma omp parallel for private(v)
	for(index_t i = 0; i < n_points; i++)
	{
		v = points[i];
		normals[i] = shape->normal(v);

		geodesics fm(shape, {v}, geodesics::FM, NULL, NIL, radio);

		indexes[i] = new index_t[fm.n_sorted_index()];

		fm.copy_sorted_index(indexes[i], fm.n_sorted_index());
		sizes[i] = fm.n_sorted_index();
	}

	return indexes;
}

bool load_sampling(vector<index_t> & points, distance_t & radio, che * mesh, size_t n)
{
	const string filename = mesh->filename();

	string file = "tmp";
	file += filename.substr(filename.find_last_of('/'), filename.size() - filename.find_last_of('/'));
	file += "." + to_string(n);

	ifstream is(file);
	debug(file)
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

		double time_fps;
		radio = farthest_point_sampling_ptp_gpu(mesh, points, time_fps, n);
		debug(time_fps)

		ofstream os(file);
		os << radio << endl;
		os << points.size() << endl;
		for(const index_t & i: points)
			os << i << endl;

		os.close();
	}

	is.close();

	return true;
}

