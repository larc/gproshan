#include <gproshan/raytracing/splat.h>

#include <gproshan/raytracing/splat_utils.h>
#include <gproshan/geometry/convex_hull.h>
#include <gproshan/pointcloud/knn.h>

#include <queue>
#include <numeric>
#include <algorithm>

#include <flann/flann.hpp>


// geometry processing and shape analysis framework
namespace gproshan::rt {


int splat::k = 1;
real_t splat::n_threshold = 0.9;
real_t splat::delta = 0.01;

splat::splat(const std::vector<che *> & pcs, const std::vector<mat4> & model_mats)
{
	for(index_t i = 0; i < pcs.size(); ++i)
		add_splats(pcs[i], model_mats[i]);
}

splat::~splat()
{
	for(che * m: pointclouds)
		delete m;
}

void splat::add_splats(che * pc, const mat4 & model_mat)
{
	std::vector<index_t> vertices;
	vertices.reserve(pc->n_vertices);

	std::vector<index_t> segs = planar_segmentation(pc, vertices, model_mat);

	gproshan_error_var(segs.size() - 1);
	gproshan_error_var(vertices.size());

	display_sets(pc, segs, vertices.data());


	std::vector<index_t> voronois[segs.size() - 1];
	std::vector<index_t> voronoi_sets[segs.size() - 1];

	#pragma omp parallel for
	for(index_t i = 1; i < segs.size(); ++i)
		voronois[i - 1] = voronoi_subdivision(voronoi_sets[i - 1], &pc->point(0), vertices, segs[i - 1], segs[i]);

	size_t n_points = 0;
	for(const auto & vs: voronoi_sets)
		n_points += vs.size();

	std::vector<index_t> splats({0});
	for(const auto & voronoi: voronois)
	for(const auto & size: voronoi)
		splats.push_back(splats.back() + size);

	gproshan_error_var(n_points);
	gproshan_error_var(splats.back());
	gproshan_error_var(splats.size() - 1);

	vertices.resize(n_points);

	n_points = 0;
	for(const auto & vs: voronoi_sets)
	{
		memcpy(vertices.data() + n_points, vs.data(), sizeof(index_t) * vs.size());
		n_points += vs.size();
	}

	gproshan_log_var(n_points);
	gproshan_log_var(vertices.size());


	init_splats(pc, model_mat, vertices, splats);
}

std::vector<index_t> splat::planar_segmentation(che * pc, std::vector<index_t> & vertices, const mat4 & model_mat)
{
	vertices.clear();
	vertices.reserve(pc->n_vertices);

	std::vector<index_t> shuffle(pc->n_vertices);
	std::iota(begin(shuffle), end(shuffle), 0);
	std::random_shuffle(begin(shuffle), end(shuffle));

	std::vector<index_t> segs({0});
	const index_t & idx = segs.size();

	std::vector<index_t> visited;
	visited.assign(pc->n_vertices, -1);


	double nn_time = 0;
	const size_t nn = 8;

/*
	TIC(nn_time);
		flann::Matrix<real_t> kpc((real_t *) &pc->point(0), pc->n_vertices, 3);

		flann::Matrix<int> indices(new int[pc->n_vertices * nn], pc->n_vertices, nn);
		flann::Matrix<real_t> dists(new real_t[pc->n_vertices * nn], pc->n_vertices, nn);

		// construct an randomized kd-tree index using 4 kd-trees
		flann::Index<flann::L2<real_t> > index(kpc, flann::KDTreeIndexParams(4));
		index.buildIndex();
	TOC(nn_time);
	gproshan_log_var(nn_time);

	TIC(nn_time);
		// do a knn search, using 128 checks
		flann::SearchParams sparams(128);
		sparams.cores = 16;
		index.knnSearch(kpc, indices, dists, nn, sparams);

		//delete [] indices.ptr();
		delete [] dists.ptr();

	TOC(nn_time);
	gproshan_log_var(nn_time);
*/

	grid_knn knn(&pc->point(0), pc->n_vertices, model_mat);
	std::vector<std::vector<index_t> > kpc(pc->n_vertices);

	TIC(nn_time);
	#pragma omp parallel for
	for(index_t v = 0; v < pc->n_vertices; ++v)
		kpc[v] = knn(vec3(model_mat * vec4(pc->point(v), 1)), nn);

	TOC(nn_time);
	gproshan_log_var(nn_time);


	vertex vnormal;
	vertex vcenter;

	std::queue<index_t> q;
	for(const index_t & v: shuffle)
	{
		if(visited[v] != NIL) continue;

		vnormal = 0;

		q.push(v);
		visited[v] = 0;

		while(!q.empty())
		{
			index_t front = q.front();
			q.pop();

			vertices.push_back(front);
			visited[front] = idx;

			const size_t & n = vertices.size() - segs.back();
			vnormal = (vnormal * (n - 1) + pc->normal(front)) / n;
			vcenter = (vcenter * (n - 1) + pc->point(front)) / n;
/*
			for(const index_t & he: pc->star(front))
			{
				const index_t & u = pc->halfedge(he_prev(he));
*/
/*
			for(index_t i = 0; i < nn; ++i)
			{
				const int & u = indices[front][i];
*/

			for(const index_t & u: kpc[front])
			{
				const vertex & p = model_mat * vec4(pc->point(u), 1);	// for adapt noisy
				if(visited[u] == NIL &&
					dot(vnormal, pc->normal(u)) > n_threshold)
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

		if(vertices.size() - segs.back() < 16)
		{
			for(index_t i = segs.back(); i < vertices.size(); ++i)
				visited[vertices[i]] = NIL;

			vertices.resize(segs.back());

			continue;
		}

		segs.push_back(vertices.size());
	}

	return segs;
}

std::vector<index_t> splat::voronoi_subdivision(std::vector<index_t> & voronoi_set,
												const vertex * points,
												const std::vector<index_t> & vertices,
												const index_t & seg_begin,
												const index_t & seg_end
												)
{
	std::vector<real_t> dist;
	dist.assign(seg_end - seg_begin, INFINITY);

	std::vector<index_t> seeds;
	seeds.push_back(vertices[seg_begin]);

	real_t radio = INFINITY;
	real_t radio_threshold = 0;
	index_t new_seed;

	while(radio > radio_threshold)
	{
		radio = 0;

		const index_t & s = seeds.back();
		for(index_t i = seg_begin; i < seg_end; ++i)
		{
			const index_t & v = vertices[i];

			real_t & vdist = dist[i - seg_begin];
			vdist = std::min(vdist, length(points[v] - points[s]));

			if(radio < vdist)
			{
				radio = vdist;
				new_seed = v;
			}
		}

		if(seeds.size() == 1)
			radio_threshold = std::max(0.2, radio * 0.1);

		seeds.push_back(new_seed);
	}


	std::vector<std::vector<index_t> > regions(seeds.size());

	for(index_t i = seg_begin; i < seg_end; ++i)
	for(index_t j = 0; j < seeds.size(); ++j)
	{
		const index_t & s = seeds[j];
		const index_t & v = vertices[i];
		const real_t & d = length(points[v] - points[s]);

		if(d < dist[i - seg_begin] + delta * radio)
			regions[j].push_back(v);
	}


	std::vector<index_t> voronoi;
	voronoi_set.clear();

	for(const auto & r: regions)
	{
		if(r.size() < 16) continue;

		for(const index_t & v: r)
			voronoi_set.push_back(v);

		voronoi.push_back(r.size());
	}

	return voronoi;
}

void splat::init_splats(const che * mesh, const mat4 & model_mat, std::vector<index_t> & vertices, const std::vector<index_t> & idx_splats)
{
	std::vector<vertex> points(vertices.size());
	std::vector<index_t> trigs;

	splats_data spc(points.size(), idx_splats.size() - 1);

	std::vector<convex_hull *> splat_chs(spc.n_splats);

	#pragma omp parallel for
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		splat_t<real_t> & s = spc.splats[i];

		s.begin = idx_splats[i];
		s.end = idx_splats[i + 1];
		vertex & center = s.center;
		mat3 & tbn = s.tbn;
		vec3 & normal = s.tbn[2];

		center = {0, 0, 0};
		normal = {0, 0, 0};
		for(index_t j = s.begin; j < s.end; ++j)
		{
			const index_t & v = vertices[j];
			vertex & p = points[j];
			p = model_mat * vec4(mesh->point(v), 1);
			center += p;
			normal += mesh->normal(v);
		}
		center /= s.end - s.begin;
		normal /= length(normal);

		tbn[0] = points[s.end - 1] - center;
		tbn[0] = normalize(tbn[0] - dot(tbn[0], tbn[2]) * tbn[2]);
		tbn[1] = normalize(cross(tbn[2], tbn[0]));

		s.radius = 0;
		for(index_t j = s.begin; j < s.end; ++j)
			s.radius = std::max(s.radius, length(points[j] - center));

		std::sort(vertices.begin() + s.begin, vertices.begin() + s.end,
					[&](const index_t & a, const index_t & b)
					{
						const vertex & p = model_mat * vec4(mesh->point(a), 1);
						const vertex & q = model_mat * vec4(mesh->point(b), 1);
						return s.morton2d(p) < s.morton2d(q);
					});

		for(index_t j = s.begin; j < s.end; ++j)
		{
			vertex & p = points[j];
			p = model_mat * vec4(mesh->point(vertices[j]), 1);
			spc.morton_codes[j] = s.morton2d(p);
			p = tbn * (p - center);
		}

		splat_chs[i] = new convex_hull(points.data() + s.begin, s.end - s.begin);

		real_t h = INFINITY;
		for(index_t j = s.begin; j < s.end; ++j)
		{
			vertex & p = points[j];
			p = mat3::transpose(tbn) * p + center;
			h = std::min(h, dot(p - center, normal));
		}

		center += h * normal;

		{
			for(index_t j = s.begin + 1; j < s.end; ++j)
			{
				if(spc.morton_codes[j - 1] > spc.morton_codes[j])
				{
					gproshan_error(FATAL ERROR);
					gproshan_log_var(spc.morton_codes[j - 1]);
					gproshan_log_var(spc.morton_codes[j]);
					break;
				}
				if(spc.morton_codes[j] >= (1 << 20))
				{
					gproshan_error(FATAL ERROR);
					gproshan_error(spc.morton_codes[j]);
					exit(0);
				}
			}
		}
	}

	std::vector<index_t> primID_splat;
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		const splat_t<real_t> & s = spc.splats[i];

		std::vector<index_t> sch = *splat_chs[i];
		for(index_t & v: sch)
		{
			vertex p = points[v + s.begin] - s.center;
			p = p - dot(p, s.tbn[2]) * s.tbn[2];
			p = p + s.center;

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
		const index_t & v = vertices[i];

		pc->heatmap(i) = mesh->heatmap(v);
		pc->normal(i) = mesh->normal(v);
		pc->rgb(i) = mesh->rgb(v);
	}

	spc.primID_splat = new unsigned int[primID_splat.size()];
	memcpy(spc.primID_splat, primID_splat.data(), sizeof(unsigned int) * primID_splat.size());

	pointclouds.push_back(pc);
	splats_pcs.emplace_back(std::move(spc));

	display_sets(pc, idx_splats);
}

void splat::display_sets(che * pc, const std::vector<index_t> & sets, const index_t * mapid)
{
	std::vector<int> color(sets.size() - 1);
	std::iota(color.begin(), color.end(), 0);
	std::random_shuffle(color.begin(), color.end());

	for(index_t i = 1; i < sets.size(); ++i)
	for(index_t j = sets[i - 1]; j < sets[i]; ++j)
		pc->heatmap(mapid ? mapid[j] : j) = real_t(color[i - 1]) / (color.size() - 1);
}


} // namespace gproshan

