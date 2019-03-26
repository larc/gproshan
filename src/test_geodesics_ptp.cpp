#include "test_geodesics_ptp.h"

#include "che_off.h"
#include "geodesics_ptp.h"
#include "heat_flow.h"

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

	FILE * ftable;
	#ifdef SINGLE_P
		ftable = fopen("ptp_results.tex", "w");
	#else
		ftable = fopen("ptp_results_double.tex", "w");
	#endif
	
	const char * str[2] = {"", "\\bf"};
	const char * ptime = "& %6.3lfs ";
	const char * pspeedup = "& \\bf (%.1lfx) ";
	const char * pbtime = "& %6s %6.3lfs ";
	const char * pberror = "& %6s %6.3lf\\%% ";

	string filename;
	while(cin >> filename)
	{
		debug(filename)
		
		vector<index_t> source = { 0 };
		
		che * mesh = new che_off(data_path + filename + ".off");
		size_t n_vertices = mesh->n_vertices();
		
		index_t * toplesets = new index_t[n_vertices];
		index_t * sorted_index = new index_t[n_vertices];
		vector<index_t> limits;
		mesh->compute_toplesets(toplesets, sorted_index, limits, source);
		
		
		// PERFORMANCE & ACCURACY ___________________________________________________________________

		double Time[7], time;	// FM, PTP GPU, HEAT cholmod, HEAT cusparse
		distance_t Error[5];	// FM, PTP GPU, HEAT cholmod, HEAT cusparse

		distance_t * exact = load_exact_geodesics(exact_dist_path + filename, n_vertices);

		Time[0] = test_fast_marching(Error[0], exact, mesh, source, n_test);
		Time[1] = test_ptp_cpu(Error[1], exact, mesh, source, limits, sorted_index, n_test);
		Time[2] = test_ptp_gpu(Error[2], exact, mesh, source, limits, sorted_index, n_test);
		
		#ifdef SINGLE_P
			Time[6] = Time[5] = Time[4] = Time[3] = INFINITY;
			Error[4] = Error[3] = INFINITY;
		#else
			Time[4] = test_heat_method_cholmod(Error[3], Time[3], exact, mesh, source, n_test);
			if(n_vertices < 100000)	
				Time[6] = test_heat_method_gpu(Error[4], Time[5], exact, mesh, source, n_test);
			else
			{
				Time[6] = Time[5] = INFINITY;
				Error[4] = INFINITY;
			}
		#endif
		
		int t_min = 0;
		for(int i = 1; i < sizeof(Time) / sizeof(double); i++)
			if(Time[t_min] > Time[i]) t_min = i;
		
		int e_min = 0;
		for(int i = 1; i < sizeof(Error) / sizeof(distance_t); i++)
			if(Error[e_min] > Error[i]) e_min = i;
		
		fprintf(ftable, "%20s ", ("\\verb|" + filename + '|').c_str());
		fprintf(ftable, "& %12lu ", n_vertices);

		// FM
		fprintf(ftable, pbtime, str[0 == t_min], Time[0]);
		fprintf(ftable, pberror, str[0 == e_min], Error[0]);

		// PTP CPU
		fprintf(ftable, pbtime, str[1 == t_min], Time[1]);
		fprintf(ftable, pspeedup, Time[0] / Time[1]);
		fprintf(ftable, pberror, str[1 == e_min], Error[1]);
		
		#ifndef SINGLE_P
			fprintf(ftable, "& OpenMP ");
		#endif

		#ifdef SINGLE_P
			// PTP GPU
			fprintf(ftable, pbtime, str[2 == t_min], Time[2]);
			fprintf(ftable, pspeedup, Time[0] / Time[2]);
			fprintf(ftable, pberror, str[2 == e_min], Error[2]);
		#endif

		#ifndef SINGLE_P
			// HEAT FLOW cholmod
			fprintf(ftable, ptime, Time[4]);
			fprintf(ftable, pbtime, str[3 == t_min], Time[3]);
			fprintf(ftable, pspeedup, Time[0] / Time[3]);
			fprintf(ftable, pberror, str[3 == e_min], Error[3]);
			fprintf(ftable, "& Cholmod ");
		#endif
		fprintf(ftable, "\\\\\n");

		#ifndef SINGLE_P
			// PTP GPU
			fprintf(ftable, "&&& ");
			fprintf(ftable, pbtime, str[2 == t_min], Time[2]);
			fprintf(ftable, pspeedup, Time[0] / Time[2]);
			fprintf(ftable, pberror, str[2 == e_min], Error[2]);
			fprintf(ftable, "& Cuda ");
			
			// HEAT FLOW cusparse
			fprintf(ftable, ptime, Time[6]);
			fprintf(ftable, pbtime, str[5 == t_min], Time[5]);
			fprintf(ftable, pspeedup, Time[0] / Time[5]);
			fprintf(ftable, pberror, str[4 == e_min], Error[4]);
			fprintf(ftable, "& cusolverSp ");

			fprintf(ftable, "\\\\\\hline\n");
		#endif

		
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
		
		vector<pair<index_t, distance_t> > iter_error = iter_error_parallel_toplesets_propagation_gpu(mesh, source, limits, sorted_index, exact, time);

		os.open(test_path + filename + "_error.iter");
		for(auto & p: iter_error)
			os << p.first << " " << p.second << endl;
		os.close();


		// FARTHEST POINT SAMPLING _________________________________________________________________
		
		size_t i_samples = source.size();
		size_t n_samples = 1001;
		double * times_fps = times_farthest_point_sampling_ptp_gpu(mesh, source, n_samples);
		
		os.open(test_path + filename + ".fps");
		for(index_t i = i_samples; i < n_samples; i++)
			os << i << " " << times_fps[i] << endl;
		os.close();

		// FREE MEMORY

		delete mesh;
		delete [] toplesets;
		delete [] sorted_index;
		delete [] exact;
		delete [] toplesets_dist;
		delete [] times_fps;
	}
	
	fclose(ftable);
}

double test_fast_marching(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test)
{
	double t, seconds = INFINITY;

	for(int i = 0; i < n_test; i++)
	{
		TIC(t) geodesics fm(mesh, source, geodesics::FM); TOC(t);
		seconds = min(seconds, t);
	}

	geodesics fm(mesh, source, geodesics::FM);

	error = compute_error(&fm[0], exact, mesh->n_vertices(), source.size());

	return seconds;
}

double test_ptp_gpu(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const vector<index_t> & limits, const index_t * sorted_index, const int & n_test)
{
	double t, seconds = INFINITY;
	
	distance_t * dist = new distance_t[mesh->n_vertices()];
	for(int i = 0; i < n_test; i++)
	{
		t = parallel_toplesets_propagation_coalescence_gpu(dist, mesh, source, limits, sorted_index);
		seconds = min(seconds, t);
	}

	error = compute_error(dist, exact, mesh->n_vertices(), source.size());

	delete [] dist;
	
	return seconds;
}

double test_ptp_cpu(distance_t & error, const distance_t * exact, che * mesh, const vector<index_t> & source, const vector<index_t> & limits, const index_t * sorted_index, const int & n_test)
{
	double t, seconds = INFINITY;
	
	distance_t * dist = new distance_t[mesh->n_vertices()];
	for(int i = 0; i < n_test; i++)
	{
		TIC(t) parallel_toplesets_propagation_cpu(dist, mesh, source, limits, sorted_index); TOC(t)
		seconds = min(seconds, t);
	}

	error = compute_error(dist, exact, mesh->n_vertices(), source.size());
	
	delete [] dist;

	return seconds;
}

double test_heat_method_cholmod(distance_t & error, double & stime, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test)
{
	double t, st, ptime;
	ptime = stime = INFINITY;
	
	distance_t * dist = NULL;
	for(int i = 0; i < n_test; i++)
	{
		if(dist) delete [] dist;
		
		TIC(t) dist = heat_flow(mesh, source, st); TOC(t)
		ptime = min(t - st, ptime);
		stime = min(st, stime);
	}

	error = compute_error(dist, exact, mesh->n_vertices(), source.size());

	delete [] dist;

	return ptime;
}

double test_heat_method_gpu(distance_t & error, double & stime, const distance_t * exact, che * mesh, const vector<index_t> & source, const int & n_test)
{
	double t, st, ptime;
	ptime = stime = INFINITY;
	
	distance_t * dist = NULL;
	for(int i = 0; i < n_test; i++)
	{
		if(dist) delete [] dist;
		
		TIC(t) dist = heat_flow_gpu(mesh, source, st); TOC(t)

		ptime = min(t - st, ptime);
		stime = min(st, stime);
	}

	error = compute_error(dist, exact, mesh->n_vertices(), source.size());

	delete [] dist;

	return ptime;
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

distance_t compute_error(const distance_t * dist, const distance_t * exact, const size_t & n, const size_t & s)
{
	distance_t error = 0;

	#pragma omp parallel for reduction(+: error)
	for(index_t v = 0; v < n; v++)
		if(exact[v] > 0) error += abs(dist[v] - exact[v]) / exact[v];

	return error * 100 / (n - s);
}

