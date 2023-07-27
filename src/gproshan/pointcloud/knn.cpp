#include <gproshan/pointcloud/knn.h>

#include <unordered_map>
#include <queue>


// geometry processing and shape analysis framework
namespace gproshan {


grid_knn::grid_knn(const point * pc, const size_t & n_points, const mat4 & transform): points(n_points)
{
	double build_time = 0;

	TIC(build_time);

	res = sqrt(n_points / 10);

	for(index_t i = 0; i < n_points; ++i)
	{
		point & p = points[i];
		p = transform * vec4(pc[i], 1);

		grid[hash(p, res)].push_back(i); 
	}

	TOC(build_time);

	gproshan_log_var(sizeof(size_t));
	gproshan_log_var(build_time);
	gproshan_log_var(res);
	gproshan_log_var(grid.size());
	gproshan_log_var(double(n_points) / grid.size());
}

std::vector<index_t> grid_knn::operator () (const point & p, int knn)
{
	const uvec3 key = hash(p, res);

	std::priority_queue<std::pair<real_t, index_t> > q;

	for(int i = -1; i < 2; ++i)
	for(int j = -1; j < 2; ++j)
	for(int k = -1; k < 2; ++k)
	{
		const uvec3 cell = {key.x() + i, key.y() + j, key.z() + k};

		if(cell.x() == NIL || cell.y() == NIL || cell.z() == NIL)
			continue;

		if(grid.find(cell) == grid.end())
			continue;

		for(const index_t & v: grid[cell])
			q.push({-length(p - points[v]), v});
	}

	std::vector<index_t> nn;
	nn.reserve(knn);

	while(!q.empty() && knn)
	{
		nn.push_back(q.top().second);
		q.pop();

		--knn;
	}

	return nn;
}


} // namespace gproshan

