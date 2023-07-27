#include <gproshan/pointcloud/knn.h>

#include <unordered_map>


// geometry processing and shape analysis framework
namespace gproshan {


knn::knn(const point * pc, const size_t & n_points, const mat4 & transform): points(n_points)
{
	double build_time = 0;

	TIC(build_time);
	
	for(index_t i = 0; i < n_points; ++i)
	{
		point & p = points[i];
		p = transform * vec4(pc[i], 1);
		
		grid[hash(p)].push_back(i); 
	}

	TOC(build_time);
	
	gproshan_log_var(build_time);
	gproshan_log_var(grid.size());
	gproshan_log_var(double(n_points) / grid.size());
}


} // namespace gproshan

