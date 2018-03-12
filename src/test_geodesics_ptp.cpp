#include "test_geodesics_ptp.h"

#include "che_off.h"
#include "geodesics_ptp.h"

void main_test_geodesics_ptp(int nargs, const char ** args)
{
	if(nargs < 3)
	{
		printf("./gproshan [path_data] [test_path]\n");
		return;
	}
	
	int n_test = 10;

	string filename;
	while(cin >> filename)
	{
		printf("%20s &", ("\\verb|" + filename + '|').c_str());

		filename = args[1] + filename;
		debug(filename)
		
		che * mesh = new che_off(filename + ".off");
		size_t n_vertices = mesh->n_vertices();
		
		index_t * toplesets = new index_t[n_vertices];
		index_t * sorted_index = new index_t[n_vertices];
		vector<index_t> limits;
		mesh->compute_toplesets(toplesets, sorted_index, limits, {0});
		
		// performance test
		float time_fm, time_ptp_cpu, time_ptp_gpu, time;
		time_fm = time_ptp_cpu = time_ptp_gpu = 0;
		
		// fast marching
		for(int t = 0; t < n_test; t++)
		{
			TIC(time) geodesics fm(mesh, {0}, geodesics::FM); TOC(time)
			time_fm += time;
		}
		time_fm /= n_test;
		
		// ptp cpu
		for(int t = 0; t < n_test; t++)
		{
			TIC(time)
			distance_t * ptp_cpu = parallel_toplesets_propagation_cpu(mesh, {0}, limits, sorted_index);
			TOC(time)
			time_ptp_cpu += time;
			delete [] ptp_cpu;
		}
		time_ptp_cpu /= n_test;
		
		// ptp gpu
		for(int t = 0; t < n_test; t++)
		{
			distance_t * ptp_gpu = parallel_toplesets_propagation_gpu(mesh, {0}, limits, sorted_index, time);
			time_ptp_gpu += time;
			delete [] ptp_gpu;
		}
		time_ptp_gpu /= n_test;
		
		printf("%18.3f &%18.3f &%18.3f &%18.3f &%18.3f ", time_fm, time_ptp_cpu, time_fm / time_ptp_cpu, time_ptp_gpu, time_fm / time_ptp_gpu);
		
		printf("\n");

		delete [] toplesets;
		delete [] sorted_index;
	}
}

