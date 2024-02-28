#include <gproshan/pointcloud/knn.h>

#include <unordered_map>
#include <queue>


// geometry processing and shape analysis framework
namespace gproshan::knn {


grid::grid(const point * pc, const size_t n_points, const mat4 & transform): points(n_points)
{
	double build_time = 0;

	TIC(build_time);

	res = sqrt(n_points / 10);

	for(index_t i = 0; i < n_points; ++i)
	{
		point & p = points[i];
		p = transform * (pc[i], 1);

		voxels[hash(p, res)].push_back(i);
	}

	TOC(build_time);

	gproshan_log_var(sizeof(size_t));
	gproshan_log_var(build_time);
	gproshan_log_var(res);
	gproshan_log_var(size(voxels));
	gproshan_log_var(double(n_points) / size(voxels));
}

std::vector<index_t> grid::operator () (const point & p, int k) const
{
	const uvec3 key = hash(p, res);

	std::priority_queue<std::pair<real_t, index_t> > q;

	for(int i = -1; i < 2; ++i)
	for(int j = -1; j < 2; ++j)
	for(int k = -1; k < 2; ++k)
	{
		const uvec3 pos = {key.x() + i, key.y() + j, key.z() + k};

		if(pos.x() == NIL || pos.y() == NIL || pos.z() == NIL)
			continue;

		const auto iter = voxels.find(pos);
		if(iter == end(voxels))
			continue;

		for(const index_t v: iter->second)
			q.push({-length(p - points[v]), v});
	}

	std::vector<index_t> nn;
	nn.reserve(k);

	while(!q.empty() && k)
	{
		nn.push_back(q.top().second);
		q.pop();

		--k;
	}

	return nn;
}


k3tree::k3tree(const point * pc, const size_t n_points, const size_t k, const std::vector<point> & query)
		: k3tree(pc, n_points, size(query) ? query.data() : nullptr, size(query), k) {}

///< Implementation using flann, by default compute all knn
k3tree::k3tree(const point * pc, const size_t n_points, const point * query, const size_t n_query, const size_t k)
{
	double time_build, time_query;

	TIC(time_build);
		flann::Matrix<real_t> mpc((real_t *) pc, n_points, 3);
		flann::Index<flann::L2<real_t> > index(mpc, flann::KDTreeSingleIndexParams());
		index.buildIndex();
	TOC(time_build);
	gproshan_log_var(time_build);

	TIC(time_query);
		const point * q = query && n_query ? query : pc;
		const size_t n_results = query && n_query ? n_query : n_points;

		flann::Matrix<real_t> mq((real_t *) q, n_results, 3);

		indices = flann::Matrix(new int[n_results * k], n_results, k);
		flann::Matrix dists(new real_t[n_results * k], n_results, k);

		flann::SearchParams params;
		params.cores = 16;
		index.knnSearch(mq, indices, dists, k, params);
	TOC(time_query);
	gproshan_log_var(time_query);

	gproshan_log_var(time_build + time_query);

	delete [] dists.ptr();
}

k3tree::~k3tree()
{
	delete [] indices.ptr();
}

const int * k3tree::operator [] (const index_t i) const
{
	return indices[i];
}

const int * k3tree::operator () (const index_t i) const
{
	return indices[i];
}

int k3tree::operator () (const index_t i, const index_t j) const
{
	return indices[i][j];
}


real_t pc_median_pairwise_distant(const point * pc, const size_t n_points, const mat4 & model_mat)
{
	k3tree p2nn(pc, n_points, 2);

	std::vector<real_t> dist(n_points);

	#pragma omp parallel for
	for(index_t i = 0; i < n_points; ++i)
		dist[i] = length(model_mat * (pc[i] - pc[p2nn(i, 1)], 0));

	std::ranges::sort(dist);

	return dist[size(dist) >> 1];
}


} // namespace gproshan

