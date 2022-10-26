#include <gproshan/raytracing/splat.h>

#include <gproshan/raytracing/splat_utils.h>

#include <queue>
#include <numeric>
#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat::splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats)
{
	for(che * m: meshes)
		add_splats_mesh(m);
}

splat::~splat()
{
	for(che * m: pointclouds)
		delete m;

	for(splats_data * spc: splats_pcs)
	{
		delete [] spc->pc;
		delete [] spc->morton_codes;
		delete [] spc->idx_splats;
	}
}

void splat::add_splats_mesh(che * mesh)
{
	const real_t n_threshold = 0.9;
	const real_t max_neigs = 1000;

	std::vector<index_t> vertices;
	std::vector<index_t> idx_splats({0});
	
	std::vector<bool> visited;
	visited.assign(mesh->n_vertices, 0);
	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(visited[v]) continue;

		const vertex & vnormal = mesh->normal(v);
		
		std::queue<index_t> q; q.push(v);
		while(!q.empty())
		{
			if(dot(vnormal, mesh->normal(q.front())) < n_threshold)
				break;

			if(vertices.size() - idx_splats.back() >= max_neigs)
				break;

			vertices.push_back(q.front());
			visited[q.front()] = true;

			for(const index_t & he: mesh->star(q.front()))
			{
				const index_t & u = mesh->halfedge(prev(he));
				if(!visited[u]) q.push(u);
			}

			q.pop();
		}
		
		// splat verification
		if(vertices.size() - idx_splats.back() < 3)
		{
			for(index_t i = idx_splats.back(); i < vertices.size(); ++i)
				visited[vertices[i]] = false;
			vertices.resize(idx_splats.back());
			continue;
		}
		
		// new splat limit
		idx_splats.push_back(vertices.size());
	}
	
	std::vector<int> color(idx_splats.size() - 1);
	std::iota(color.begin(), color.end(), 0);
	std::random_shuffle(color.begin(), color.end());
	for(index_t i = 1; i < idx_splats.size(); ++i)
	for(index_t j = idx_splats[i - 1]; j < idx_splats[i]; ++j)
		mesh->heatmap(vertices[j]) = real_t(color[i - 1]) / (color.size() - 1);
	
//	pointclouds.push_back(pc);
//	splats_pcs.push_back(spc);
}

void splat::build_splats_ch(splats_data * s, const che * pc, const std::vector<index_t> & vertices)
{

}


} // namespace gproshan

