#include "test_geodesics_ptp.h"

#include "che_off.h"
#include "geodesics_ptp.h"

#include <cassert>

void main_test_geodesics_ptp(const int & nargs, const char ** args)
{
	if(nargs < 4)
	{
		printf("./test_geodesics [data_path] [test_path] [exact_dist_path] [n_test = 10]\n");
		return;
	}
	
	const char * data_path = args[1];
	const char * test_path = args[2];
	const char * exact_dist_path = args[3];

	int n_test = nargs == 5 ? atoi(args[4]) : 10;
	bool cpu = 0;

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

		
		// PERFORMANCE _____________________________________________________________________________
		
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
		if(cpu)
		{
			for(int t = 0; t < n_test; t++)
			{
				if(ptp_cpu) delete [] ptp_cpu;

				TIC(time)
				ptp_cpu = parallel_toplesets_propagation_cpu(mesh, source, limits, sorted_index);
				TOC(time)
				time_ptp_cpu += time;
			}
			time_ptp_cpu /= n_test;
		}
		
		// ptp gpu
		for(int t = 0; t < n_test; t++)
		{
			if(ptp_gpu) delete [] ptp_gpu;

			ptp_gpu = parallel_toplesets_propagation_coalescence_gpu(mesh, source, limits, sorted_index, time);
			time_ptp_gpu += time;
		}
		time_ptp_gpu /= n_test;
		

		printf("%18.3f &", time_fm);
		if(cpu) printf("%18.3f &%18.3f &", time_ptp_cpu, time_fm / time_ptp_cpu);
		printf("%18.3f &%18.3f &", time_ptp_gpu, time_fm / time_ptp_gpu);
		
		
		// ACCURACY ________________________________________________________________________________

		distance_t * exact = load_exact_geodesics(exact_dist_path + filename, n_vertices);
		geodesics fm(mesh, source, geodesics::FM);

		distance_t error_fm, error_ptp_cpu, error_ptp_gpu;

		error_fm = error_ptp_cpu = error_ptp_gpu = 0;
		for(index_t v = 0; v < n_vertices; v++)
		{
			if(exact[v] > 0)
			{
				error_fm += abs(fm[v] - exact[v]) / exact[v];
				if(cpu) error_ptp_cpu += abs(ptp_cpu[v] - exact[v]) / exact[v];
				error_ptp_gpu += abs(ptp_gpu[v] - exact[v]) / exact[v];
			}
		}

		error_fm /= n_vertices;
		if(cpu) error_ptp_cpu /= n_vertices;
		error_ptp_gpu /= n_vertices;

		printf("%18.3e &", error_fm * 100);
		if(cpu) printf("%18.3e &", error_ptp_cpu * 100);
		printf("%18.3e ", error_ptp_gpu * 100);
		printf("\\\\\\hline\n");
		
		
		// DEGREE HISTOGRAM ________________________________________________________________________

		index_t dv;
		map<index_t, index_t> deg;
		for(index_t v = 0; v < n_vertices; v++)
		{
			dv = mesh->ot_evt(v) == NIL ? 1 : 0;
			for_star(he, mesh, v) dv++;
			deg[dv]++;
		}

		ofstream os(test_path + filename + ".deg");
		for(auto & ii: deg)
			os << ii.first << " " << ii.second << endl;
		os.close();
			

		// TOPLESETS DISTRIBUTION __________________________________________________________________

		index_t * toplesets_dist = new index_t[limits.size() - 1];

		os.open(test_path + filename + "_toplesets.dist");
		for(index_t i = 1; i < limits.size(); i++)
		{
			toplesets_dist[i - 1] = limits[i] - limits[i - 1];
			os << i - 1 << " " << toplesets_dist[i - 1] << endl;
		}
		os.close();

		sort(toplesets_dist, toplesets_dist + limits.size() - 1);

		os.open(test_path + filename + "_toplesets_sorted.dist");
		for(index_t i = 0; i < limits.size() - 1; i++)
			os << i << " " << toplesets_dist[i] << endl;
		os.close();

		
		// PTP ITERATION ERROR _____________________________________________________________________
		
		distance_t * iter_error = iter_error_parallel_toplesets_propagation_gpu(mesh, source, limits, sorted_index, exact, time);

		os.open(test_path + filename + "_error.iter");
		index_t n_iter = iterations(limits);
		for(index_t j = 0, i = limits.size(); i < n_iter; i++, j++)
			os << i << " " << iter_error[j] << endl;
		os.close();

		// FARTHEST POINT SAMPLING _________________________________________________________________
		
		size_t i_samples = source.size();
		size_t n_samples = 1001;
		float * times_fps = times_farthest_point_sampling_ptp_gpu(mesh, source, n_samples);
		
		os.open(test_path + filename + ".fps");
		for(index_t i = i_samples; i < n_samples; i++)
			os << i << " " << times_fps[i] << endl;
		os.close();

		// FREE MEMORY

		delete mesh;
		delete [] toplesets;
		delete [] sorted_index;
		if(cpu) delete [] ptp_cpu;
		delete [] ptp_gpu;
		delete [] exact;
		delete [] toplesets_dist;
		delete [] iter_error;
		delete [] times_fps;
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

