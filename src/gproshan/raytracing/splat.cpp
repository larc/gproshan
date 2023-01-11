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
	const real_t n_threshold = 0.9;

	std::vector<index_t> vertices;
	std::vector<index_t> segmentation({0});
	std::vector<vertex> centers;

	std::vector<index_t> visited;
	visited.assign(mesh->n_vertices, -1);

	std::vector<index_t> shuffle(mesh->n_vertices);
	std::iota(begin(shuffle), end(shuffle), 0);
	std::random_shuffle(begin(shuffle), end(shuffle));

	const index_t & idx = segmentation.size();

	std::queue<index_t> q;
	for(const index_t & v: shuffle)
	{
		if(visited[v] != NIL) continue;

		const vertex & vnormal = mesh->normal(v);

		q.push(v);
		visited[v] = 0;

		while(!q.empty())
		{
			index_t front = q.front();
			q.pop();

			vertices.push_back(front);
			visited[front] = idx;

			for(const index_t & he: mesh->star(front))
			{
				const index_t & u = mesh->halfedge(he_prev(he));
				if(visited[u] == NIL &&
					dot(vnormal, mesh->normal(front)) > n_threshold)
				{
					q.push(u);
					visited[u] = 0;
				}

			}
		}

		while(!q.empty())
		{
			visited[q.front()] = NIL;
			q.pop();
		}

		if(vertices.size() - segmentation.back() < 3)
		{
			for(index_t i = segmentation.back(); i < vertices.size(); ++i)
				visited[vertices[i]] = NIL;

			vertices.resize(segmentation.back());

			continue;
		}

		segmentation.push_back(vertices.size());
	}

	gproshan_error_var(vertices.size());
	gproshan_error_var(segmentation.size());

	std::vector<index_t> idx_splats({0});
	visited.assign(mesh->n_vertices, -1);

	std::vector<index_t> seeds;
	std::vector<std::vector<index_t> > voronoi;
	std::vector<real_t> dist;
	dist.assign(mesh->n_vertices, INFINITY);

	for(index_t i = 1; i < segmentation.size(); ++i)
	{
		const index_t & begin = segmentation[i - 1];
		const index_t & end = segmentation[i];

		seeds.clear();
		seeds.push_back(vertices[begin]);

		real_t radio = INFINITY;
		real_t radio_threshold = 0;
		index_t next_seed;
		while(radio > radio_threshold)
		{
			radio = 0;

			const index_t & s = seeds.back();
			for(index_t j = begin; j < end; ++j)
			{
				const index_t & u = vertices[j];
				const real_t & d = length(mesh->point(u) - mesh->point(s));

				if(d < dist[u])
				{
					visited[u] = seeds.size() - 1;
					dist[u] = d;
				}

				if(radio < dist[u])
				{
					radio = dist[u];
					next_seed = u;
				}
			}

			if(seeds.size() == 1)
				radio_threshold = std::max(0.2, radio * 0.1);

			seeds.push_back(next_seed);
		}

		voronoi.assign(seeds.size(), {});
		for(index_t j = begin; j < end; ++j)
		{
			const index_t & v = vertices[j];
			voronoi[visited[v]].push_back(v);
		}

		for(auto & region: voronoi)
		{
			if(region.size() < 3) continue;

			for(index_t i = 0; i < region.size(); ++i)
				vertices[i + idx_splats.back()] = region[i];
			idx_splats.push_back(idx_splats.back() + region.size());
		}
	}


	auto display = [&mesh, &vertices](const std::vector<index_t> & sets)
	{
		std::vector<int> color(sets.size() - 1);
		std::iota(color.begin(), color.end(), 0);
		std::random_shuffle(color.begin(), color.end());

		for(index_t i = 1; i < sets.size(); ++i)
		for(index_t j = sets[i - 1]; j < sets[i]; ++j)
			mesh->heatmap(vertices[j]) = real_t(color[i - 1]) / (color.size() - 1);
	};

	display(idx_splats);

	gproshan_error_var(idx_splats.size());

	init_splats(mesh, model_mat, vertices, idx_splats);
}

void splat::init_splats(const che * mesh, const mat4 & model_mat, std::vector<index_t> & vertices, const std::vector<index_t> & idx_splats)
{
	std::vector<vertex> points(vertices.size());
	std::vector<index_t> trigs;

	splats_data * spc = new splats_data;
	spc->morton_codes = new unsigned int[mesh->n_vertices];
	spc->n_splats = idx_splats.size() - 1;
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
		vec3 & normal = tbn[2];

		center = {0, 0, 0};
		normal = {0, 0, 0};
		for(index_t j = begin; j < end; ++j)
		{
			const index_t & v = vertices[j];
			center += mesh->point(v);
			normal += mesh->normal(v);
		}
		center /= end - begin;
		normal /= length(normal);

		tbn[0] = points[end - 1] - center;
		tbn[0] = normalize(tbn[0] - dot(tbn[0], tbn[2]) * tbn[2]);
		tbn[1] = normalize(cross(tbn[2], tbn[0]));

		center = model_mat * vec4(center, 1);
		for(index_t j = begin; j < end; ++j)
		{
			const index_t & v = vertices[j];
			vertex p = model_mat * vec4(mesh->point(v), 1);
			p = tbn * (p - center);
			spc->morton_codes[v] = morton_2d((p.x() + 1) / 2, (p.y() + 1) / 2);
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
		}

		splat_chs[i] = new convex_hull(points.data() + begin, end - begin);

		for(index_t j = begin; j < end; ++j)
		{
			vertex & p = points[j];
			p = mat3::transpose(tbn) * p + center;
		}
	}

	std::vector<index_t> primID_splat;
	for(index_t i = 0; i < spc->n_splats; ++i)
	{
		const index_t & begin = idx_splats[i];
		const vertex & center = spc->centers[i];
		const mat3 & tbn = spc->tbns[i];

		std::vector<index_t> sch = *splat_chs[i];
		for(index_t & v: sch)
		{
			vertex p = points[v + begin] - center;
			p = p - dot(p, tbn[2]) * tbn[2];
			p = p + center;

			v = points.size();
			points.push_back(p);
		}

		index_t f = -1;
		for(const index_t & v: che::trig_convex_polygon(sch.data(), sch.size()))
		{
			trigs.push_back(v);
			if(!(++f % 3))
				primID_splat.push_back(i);
		}
	}

	for(convex_hull * ch: splat_chs)
		delete ch;

	gproshan_error_var(points.size());
	gproshan_error_var(trigs.size());

	che * pc = new che(points.data(), points.size(), trigs.data(), trigs.size() / 3);

	#pragma omp parallel for
	for(index_t i = 0; i < vertices.size(); ++i)
	{
		pc->heatmap(i) = mesh->heatmap(vertices[i]);
		pc->normal(i) = mesh->normal(vertices[i]);
		pc->rgb(i) = mesh->rgb(vertices[i]);
	}

	spc->primID_splat = new unsigned int[primID_splat.size()];
	memcpy(spc->primID_splat, primID_splat.data(), sizeof(unsigned int) * primID_splat.size());

	spc->idx_splats = new unsigned int[idx_splats.size()];
	memcpy(spc->idx_splats, idx_splats.data(), sizeof(unsigned int) * idx_splats.size());

	pointclouds.push_back(pc);
	splats_pcs.push_back(spc);
}


} // namespace gproshan

