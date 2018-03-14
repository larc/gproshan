#include "test_geodesics_ptp.h"

#include "che_off.h"
#include "geodesics_ptp.h"

#include <cassert>

void main_test_geodesics_ptp(int nargs, const char ** args)
{
	if(nargs < 4)
	{
		printf("./gproshan [data_path] [test_path] [exact_dist_path]\n");
		return;
	}
	
	const char * data_path = args[1];
	const char * test_path = args[2];
	const char * exact_dist_path = args[3];

	int n_test = 10;

	string filename;
	while(cin >> filename)
	{
		debug(filename)
		
		che * mesh = new che_off(data_path + filename + ".off");
		size_t n_vertices = mesh->n_vertices();
		printf("%20s & %12lu &", ("\\verb|" + filename + '|').c_str(), n_vertices);
		
		index_t * toplesets = new index_t[n_vertices];
		index_t * sorted_index = new index_t[n_vertices];
		vector<index_t> limits;
		mesh->compute_toplesets(toplesets, sorted_index, limits, {0});
		
		vector<index_t> source = { 0 };

		// PERFORMANCE
		float time_fm, time_ptp_cpu, time_ptp_gpu, time;
		time_fm = time_ptp_cpu = time_ptp_gpu = 0;
		
		// fast marching
		for(int t = 0; t < n_test; t++)
		{
			TIC(time) geodesics fm(mesh, source, geodesics::FM); TOC(time)
			time_fm += time;
		}
		time_fm /= n_test;
		
		distance_t * ptp_cpu, * ptp_gpu;
		ptp_cpu = ptp_gpu = NULL;

		// ptp cpu
		for(int t = 0; t < n_test; t++)
		{
			if(ptp_cpu) delete [] ptp_cpu;

			TIC(time)
			ptp_cpu = parallel_toplesets_propagation_cpu(mesh, source, limits, sorted_index);
			TOC(time)
			time_ptp_cpu += time;
		}
		time_ptp_cpu /= n_test;
		
		// ptp gpu
		for(int t = 0; t < n_test; t++)
		{
			if(ptp_gpu) delete [] ptp_gpu;

			ptp_gpu = parallel_toplesets_propagation_gpu(mesh, source, limits, sorted_index, time);
			time_ptp_gpu += time;
		}
		time_ptp_gpu /= n_test;
		
		printf("%18.3f &%18.3f &%18.3f &%18.3f &%18.3f &", time_fm, time_ptp_cpu, time_fm / time_ptp_cpu, time_ptp_gpu, time_fm / time_ptp_gpu);
		
		// ACCURACY
		distance_t * exact = load_exact_geodesics(exact_dist_path + filename, n_vertices);
		geodesics fm(mesh, source, geodesics::FM);

		distance_t error_fm, error_ptp_cpu, error_ptp_gpu;

		error_fm = error_ptp_cpu = error_ptp_gpu = 0;
		for(index_t v = 0; v < n_vertices; v++)
		{
			if(exact[v] > 0)
			{
				error_fm += abs(fm[v] - exact[v]) / exact[v];
				error_ptp_cpu += abs(ptp_cpu[v] - exact[v]) / exact[v];
				error_ptp_gpu += abs(ptp_gpu[v] - exact[v]) / exact[v];
			}
		}

		error_fm /= n_vertices;
		error_ptp_cpu /= n_vertices;
		error_ptp_gpu /= n_vertices;

		printf("%18.3e &%18.3e &%18.3e", error_fm * 100, error_ptp_cpu * 100, error_ptp_gpu * 100);
		printf("\\\\\\hline\n");

		delete mesh;
		delete [] toplesets;
		delete [] sorted_index;
		delete [] ptp_cpu;
		delete [] ptp_gpu;
		delete [] exact;
	}
}

distance_t * load_exact_geodesics(const string & file, const size_t & n)
{
	ifstream is(file);
	assert(is.good());

	distance_t * exact = new distance_t[n];

	for(index_t i = 0; i < n; i++)
		is >> exact[i];
	is.close();
	
	return exact;
}

