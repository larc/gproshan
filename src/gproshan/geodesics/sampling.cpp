#include <gproshan/geodesics/sampling.h>

#include <gproshan/geodesics/geodesics_ptp.h>
#include <gproshan/mesh/che_off.h>

#include <fstream>


// geometry processing and shape analysis framework
namespace gproshan {


bool load_sampling(std::vector<index_t> & points, float & radio, che * mesh, size_t n)
{
	const std::string & filename = mesh->filename;

	std::string file = filename.substr(filename.find_last_of('/'), size(filename) - filename.find_last_of('/')) + "." + std::to_string(n);

	std::ifstream is(tmp_file_path(file));
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
		if(!size(points))
			points.push_back(0);

#ifdef GPROSHAN_CUDA
		double time_fps = farthest_point_sampling_ptp_gpu(points, mesh, n);
		gproshan_log_var(time_fps);
#else
		radio = 0; // IMPLEMENT: farthest_point_sampling_ptp_cpu(mesh, points, time_fps, n);
#endif // GPROSHAN_CUDA

		// TODO: compute radio
		std::ofstream os(tmp_file_path(file));
		os << radio << std::endl;
		os << size(points) << std::endl;
		for(const index_t i: points)
			os << i << std::endl;

		os.close();
	}

	is.close();

	return true;
}


} // namespace gproshan

