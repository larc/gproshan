#include "sampling.h"
#include "test.h"
#include "che_off.h"

#include <fstream>

inline index_t farthest(distance_t * d, size_t n)
{
	index_t f = 0;
	
	#pragma omp parallel for
	for(index_t v = 0; v < n; v++)
		#pragma omp critical
		if(d[v] < INFINITY && d[f] < d[v])
			f = v;
	
	return f;
}

distance_t parallel_farthest_point_sampling(vector<index_t> & points, che * shape, size_t n, distance_t radio)
{
	float time;
	index_t f;
	distance_t max_dis = INFINITY;

	ofstream os(PATH_TEST + "fastmarching/" + shape->name() + ".fps");
	
	index_t * rings = new index_t[shape->n_vertices()];
	index_t * sorted = new index_t[shape->n_vertices()];

	if(n >= shape->n_vertices())
		n = shape->n_vertices() / 2;

	n -= points.size() - 1;
	
	while(n-- && max_dis > radio)
	{
		vector<index_t> limites;
		shape->sort_by_rings(rings, sorted, limites, points);
		
		distance_t * distances = parallel_fastmarching(shape, points.data(), points.size(), time, limites, sorted);
		
		f = farthest(distances, shape->n_vertices());

		os << points.size() << " " << time << endl;
		if(n)
		{
			points.push_back(f);
			max_dis = distances[f];
		}

		delete [] distances;
	}

	os.close();

	delete [] rings;
	delete [] sorted;

	return max_dis;
}

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

		geodesics fm(shape, { v }, NIL, radio);
		
		indexes[i] = new index_t[fm.n_sorted_index()];
		
		fm.copy_sorted_index(indexes[i], fm.n_sorted_index());
		sizes[i] = fm.n_sorted_index();

	}

	return indexes;
}



void sampling_shape(int nargs, char ** args)
{
	if(nargs < 2) return;

	cout<<__FUNCTION__<<endl;
	
	string file = args[1];

	che * shape = new che_off(PATH_MDATA + file);

	size_t n_points;
	cin >> n_points;
	

	vector<index_t> points;
	points.push_back(0);
	cin >> points[0];

	size_t * sizes;
	vertex * normals;

	double time = omp_get_wtime();
	
	distance_t max_dis = parallel_farthest_point_sampling(points, shape, n_points);
	
	time = omp_get_wtime() - time;
	cout<<"time fps: "<<time<<endl;
	cout<<"max_dis: "<<max_dis<<endl;
	
	index_t ** indexes = sampling_shape(points, sizes, normals, shape, n_points, 1.5 * max_dis);
	
	cout<<n_points<<endl;

	ofstream os(PATH_TEST + file + "-indexes");
	ofstream pos(PATH_TEST + file + "-normals");

	for(index_t v, i = 0; i < n_points; i++)
	{
		pos<<normals[i]<<endl;

		v = points[i];

		for(index_t k = 0; k < sizes[i]; k++)
			os<<" "<<indexes[i][k];

		os<<endl;
	}

	os.close();
	pos.close();

	delete shape;
	delete [] normals;
	for(index_t i = 0; i < n_points; i++)
		delete [] indexes[i];
	delete [] indexes;
	delete [] sizes;
}

void sampling_shape(const char * name)
{
	cout<<__FUNCTION__<<endl;
	
	string file = name;

	che * shape = new che_off(PATH_MDATA + file);

	size_t n_points, K;
	cin >> n_points >> K;
	

	vector<index_t> points;
	points.push_back(0);
	size_t * sizes;
	vertex * normals;

	double time = omp_get_wtime();
	
	distance_t max_dis = parallel_farthest_point_sampling(points, shape, n_points);
	
	time = omp_get_wtime() - time;
	cout<<"time fps: "<<time<<endl;
	cout<<"max_dis: "<<max_dis<<endl;

	index_t ** indexes = sampling_shape(points, sizes, normals, shape, n_points, 1.5 * max_dis);
	
	cout<<n_points<<endl;

	ofstream os(PATH_TEST + file + "-indexes");
	ofstream pos(PATH_TEST + file + "-normals");

	size_t max_K = 0;
	for(index_t i = 0; i < n_points; i++)
		max_K = max_K < sizes[i] ? sizes[i] : max_K;

	cout<<"K: "<<max_K<<endl;
	max_K = K;
	size_t n_, n, m;
	index_t k;
	for(index_t v, i = 0; i < n_points; i++)
	{
		pos<<normals[i]<<endl;

		v = points[i];
		n = max_K / sizes[i];
		m = max_K % sizes[i];

		int c = 0;
		for(k = 0; k < sizes[i]; k++)
		{
			n_ = n;
			while(n_--) os<<" "<<indexes[i][k];
			if(m)
			{
				os<<" "<<indexes[i][k];
				m--;
			}
		}

		os<<endl;
	}

	os.close();
	pos.close();

	delete shape;
	delete [] normals;
	for(index_t i = 0; i < n_points; i++)
		delete [] indexes[i];
	delete [] indexes;
	delete [] sizes;
}

bool load_sampling(vector<index_t> & points, distance_t & radio, che * mesh, size_t M)
{
	const string filename = mesh->filename();
	string file = PATH_TEST;
	file += "sampling";
	file += filename.substr(filename.find_last_of('/'), filename.size() - filename.find_last_of('/'));
	file += "." + to_string(M);

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
		is.close();

		if(!points.size())
			points.push_back(0);

		float time;
//		radio = parallel_farthest_point_sampling(points, mesh, M);
		TIC(time)
		radio = farthest_point_sampling_gpu(points, time, mesh, M, 0);
		TOC(time)
		debug(time)

		ofstream os(file);
		os << radio << endl;
		os << points.size() << endl;
		for(index_t i: points)
			os << i << endl;

		os.close();
	}

	return true;
}

