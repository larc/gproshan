#include <gproshan/raytracing/splat.h>

#include <gproshan/raytracing/splat_utils.h>
#include <gproshan/geometry/convex_hull.h>

#include <queue>
#include <numeric>
#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan::rt {


splat::splat(const std::vector<che *> & meshes, const std::vector<mat4> & model_mats)
{
	for(index_t i = 0; i < meshes.size(); ++i)
		add_splats_mesh(meshes[i], model_mats[i]);
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
		delete [] spc->tbns;
		delete [] spc->centers;
	}
}

void splat::add_splats_mesh(che * mesh, const mat4 & model_mat)
{
	const real_t n_threshold = 0.81;
	const real_t max_neigs = 1000;

	std::vector<index_t> vertices;
	std::vector<index_t> idx_splats({0});

	std::vector<unsigned int> visited;
	visited.assign(mesh->n_vertices, -1);

	for(index_t v = 0; v < mesh->n_vertices; ++v)
	{
		if(visited[v] != NIL) continue;

		const vertex & vnormal = mesh->normal(v);
		const vertex & vpoint = mesh->point(v);

		std::queue<index_t> q; q.push(v);

		while(!q.empty())
		{
			index_t front = q.front();
			q.pop();

			if(visited[front] == idx_splats.size())
				continue;

			if(visited[front] != NIL)
				break;

			if(vertices.size() - idx_splats.back() >= max_neigs)
				break;
			
			real_t dist = length(vec3(model_mat * vec4(vpoint, 1) - model_mat * vec4(mesh->point(front), 1))) / (2 * M_SQRT2);
			if(dot(vnormal, mesh->normal(front)) < dist) // vs n_threshold
				break;

			vertices.push_back(front);
			visited[front] = idx_splats.size();

			for(const index_t & he: mesh->star(front))
			{
				const index_t & u = mesh->halfedge(prev(he));
				if(visited[u] == NIL) q.push(u);
			}
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

	gproshan_error_var(vertices.size());
	gproshan_error_var(idx_splats.size());

	std::vector<vertex> points(vertices.size());
	std::vector<index_t> faces;
	
	#pragma omp parallel for
	for(index_t i = 0; i < vertices.size(); ++i)
		points[i] = model_mat * vec4(mesh->point(vertices[i]), 1);
	
	splats_data * spc = new splats_data;
	spc->morton_codes = new unsigned int[vertices.size()];
	spc->n_splats = idx_splats.size() - 1;
	spc->idx_splats = new unsigned int[spc->n_splats];
	spc->tbns = new mat3[spc->n_splats];
	spc->centers = new vertex[spc->n_splats];
	
	std::vector<convex_hull *> splat_chs(spc->n_splats);
	
	#pragma omp parallel for
	for(index_t i = 0; i < spc->n_splats; ++i)
	{
		const unsigned int & begin = idx_splats[i];
		const unsigned int & end = idx_splats[i + 1];
		vertex & center = spc->centers[i];
		mat3 & tbn = spc->tbns[i];

		center = points[begin];
		tbn[2] = mesh->normal(vertices[begin]);
		tbn[0] = points[end - 1] - center;
		tbn[0] = normalize(tbn[0] - dot(tbn[0], tbn[2]) * tbn[2]);
		tbn[1] = normalize(tbn[2] * tbn[0]);
		
		for(index_t j = begin; j < end; ++j)
		{
			vertex & p = points[j];
			p = tbn * (p - center);
			spc->morton_codes[j] = morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2);
		}

		std::sort(vertices.begin() + begin, vertices.begin() + end,
					[&](const index_t & a, const index_t & b)
					{
						return spc->morton_codes[a] < spc->morton_codes[b];
					});

		for(index_t j = begin; j < end; ++j)
		{
			vertex & p = points[j];
			p = model_mat * vec4(mesh->point(vertices[j]), 1);
			p = tbn * (p - center);
			spc->morton_codes[j] = morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2);
		}
		
		splat_chs[i] = new convex_hull(points.data() + begin, end - begin);
	}
	
	for(convex_hull * ch: splat_chs)
		delete ch;
//	pointclouds.push_back(pc);
//	splats_pcs.push_back(spc);
}

void splat::build_splats_ch(splats_data * s, const che * pc, const std::vector<index_t> & vertices)
{

}


} // namespace gproshan

