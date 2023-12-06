#include <gproshan/raytracing/splat.h>

#include <gproshan/raytracing/splat_utils.h>
#include <gproshan/geometry/convex_hull.h>

#include <queue>
#include <numeric>
#include <algorithm>


// geometry processing and shape analysis framework
namespace gproshan::rt {


size_t splat::k_nn = 9;
real_t splat::t_normal = 0.9;
real_t splat::d_overlap = 0.1;

splat::splat(const std::vector<che *> & pcs, const std::vector<mat4> & model_mats)
{
	for(index_t i = 0; i < size(pcs); ++i)
		add_splats(pcs[i], model_mats[i]);
}

splat::~splat()
{
	for(che * m: pointclouds)
		delete m;
}

void splat::add_splats(che * pc, const mat4 & model_mat)
{
	TIC(time);
	knn::k3tree k3tree(&pc->point(0), pc->n_vertices, splat::k_nn);
	TOC(time);
	time_knn += time;


	TIC(time);
	std::vector<index_t> vertices;
	vertices.reserve(pc->n_vertices);

	const std::vector<index_t> segs = planar_segmentation(pc, vertices, k3tree);
	TOC(time);
	time_segmentation += time;


	display_sets(pc, segs, vertices.data());
	gproshan_error_var(size(segs) - 1);
	gproshan_error_var(size(vertices));


	TIC(time);
	std::vector<index_t> voronois[size(segs) - 1];
	std::vector<index_t> voronoi_sets[size(segs) - 1];

	#pragma omp parallel for
	for(index_t i = 1; i < size(segs); ++i)
		voronois[i - 1] = voronoi_subdivision(voronoi_sets[i - 1], vertices, &pc->point(0), k3tree, segs[i - 1], segs[i]);
	TOC(time);
	time_subdivision += time;

	gproshan_error_var(time);

	TIC(time);
	size_t n_points = 0;
	for(const auto & vs: voronoi_sets)
		n_points += size(vs);

	std::vector<index_t> splats({0});
	for(const auto & voronoi: voronois)
	for(const auto & size: voronoi)
		splats.push_back(splats.back() + size);

	vertices.resize(n_points);

	n_points = 0;
	for(const auto & vs: voronoi_sets)
	{
		memcpy(vertices.data() + n_points, vs.data(), sizeof(index_t) * size(vs));
		n_points += size(vs);
	}

	che * new_pc = init_splats(pc, model_mat, vertices, splats);
	pointclouds.push_back(new_pc);
	TOC(time);
	time_initsplats += time;


	display_sets(new_pc, splats);

	gproshan_error_var(size(vertices));
	gproshan_error_var(splats.back());
	gproshan_error_var(size(splats) - 1);

	gproshan_error_var(new_pc->n_vertices);
	gproshan_error_var(new_pc->n_trigs);

	time = time_knn + time_segmentation + time_subdivision + time_initsplats;

	gproshan_error_var(time);
}

std::vector<index_t> splat::planar_segmentation(const che * pc, std::vector<index_t> & vertices, const knn::k3tree & k3tree)
{
	vertices.clear();
	vertices.reserve(pc->n_vertices);

	std::vector<index_t> shuffle(pc->n_vertices);
	std::iota(begin(shuffle), end(shuffle), 0);

	std::random_device rd;
	std::mt19937 gen{rd()};
	std::shuffle(begin(shuffle), end(shuffle), gen);

	std::vector<index_t> segs({0});
	const index_t & idx = size(segs);

	std::vector<index_t> visited;
	visited.assign(pc->n_vertices, -1);

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

			const size_t & n = size(vertices) - segs.back();
			vnormal = (vnormal * (n - 1) + pc->normal(front)) / n;
			vcenter = (vcenter * (n - 1) + pc->point(front)) / n;

			const int * nn = k3tree(front);
			for(index_t i = 0; i < splat::k_nn; ++i)
			{
				const int & u = nn[i];

//				const vertex & p = model_mat * (pc->point(u), 1);	// for adapt noisy
				if(visited[u] == NIL && dot(vnormal, pc->normal(u)) > splat::t_normal)
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

		if(size(vertices) - segs.back() < splat::k_nn)
		{
			for(index_t i = segs.back(); i < size(vertices); ++i)
				visited[vertices[i]] = NIL;

			vertices.resize(segs.back());

			continue;
		}

		segs.push_back(size(vertices));
	}

	return segs;
}

std::vector<index_t> splat::voronoi_subdivision(std::vector<index_t> & voronoi_set,
												const std::vector<index_t> & vertices,
												const vertex * points,
												const knn::k3tree & k3tree,
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
	index_t new_seed = NIL;

	const size_t max_seeds = 3 * (log10(seg_end - seg_begin) + 1); 

	while(radio > radio_threshold || seeds.size() < max_seeds)
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

		if(size(seeds) == 1)
			radio_threshold = std::max(0.2, radio * 0.1);

		if(new_seed != NIL)
		{
			seeds.push_back(new_seed);
			new_seed = NIL;
		}
	}


	std::vector<std::vector<index_t> > regions(size(seeds));

	bool in = false;
	for(index_t i = seg_begin; i < seg_end; ++i)
	for(index_t j = 0; j < size(seeds); ++j)
	{
		const index_t & s = seeds[j];
		const index_t & v = vertices[i];

		in = false;

		const int * nn = k3tree(v);
		for(index_t k = 0; k < splat::k_nn; ++k)
		{
			const int & u = nn[k];
			const real_t & d = length(points[u] - points[s]);
			in |= d < (dist[i - seg_begin] + 1e-5);
		}

		if(in) regions[j].push_back(v);
	}


	std::vector<index_t> voronoi;
	voronoi_set.clear();

	for(const auto & r: regions)
	{
		if(size(r) < splat::k_nn) continue;

		for(const index_t & v: r)
			voronoi_set.push_back(v);

		voronoi.push_back(size(r));
	}

	return voronoi;
}

che * splat::init_splats(const che * mesh, const mat4 & model_mat, std::vector<index_t> & vertices, const std::vector<index_t> & idx_splats)
{
	std::vector<vertex> points(size(vertices));
	std::vector<index_t> trigs;

	splats_data spc(size(idx_splats) - 1);

	std::vector<convex_hull *> splat_chs(spc.n_splats);

	#pragma omp parallel for
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		auto & s = spc.splats[i];

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
			p = model_mat * (mesh->point(v), 1);
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

		std::sort(begin(vertices) + s.begin, begin(vertices) + s.end,
					[&](const index_t & a, const index_t & b)
					{
						const vertex & p = model_mat * (mesh->point(a), 1);
						const vertex & q = model_mat * (mesh->point(b), 1);
						return s.morton2d(p) < s.morton2d(q);
					});

		for(index_t j = s.begin; j < s.end; ++j)
		{
			vertex & p = points[j];
			p = model_mat * (mesh->point(vertices[j]), 1);
			p = tbn * (p - center);
		}

		splat_chs[i] = new convex_hull(points.data() + s.begin, s.end - s.begin);

//		real_t h = INFINITY;
		for(index_t j = s.begin; j < s.end; ++j)
		{
			vertex & p = points[j];
			p = mat3::transpose(tbn) * p + center;
//			h = std::min(h, dot(p - center, normal));
		}

//		center += h * normal;
	}

	std::vector<index_t> primID_splat;
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		const auto & s = spc.splats[i];

		std::vector<index_t> sch = *splat_chs[i];
		for(index_t & v: sch)
		{
			vertex p = points[v + s.begin] - s.center;
			p = p - dot(p, s.tbn[2]) * s.tbn[2];
			p = p + s.center;

			v = size(points);
			points.push_back(p);
		}

		index_t f = -1;
		for(const index_t & v: che::trig_convex_polygon(sch.data(), size(sch)))
		{
			trigs.push_back(v);
			if(!(++f % 3))
				primID_splat.push_back(i);
		}
	}

	for(convex_hull * ch: splat_chs)
		delete ch;

	che * pc = new che(points.data(), size(points), trigs.data(), size(trigs) / 3);

	#pragma omp parallel for
	for(index_t i = 0; i < size(vertices); ++i)
	{
		const index_t & v = vertices[i];

		pc->heatmap(i) = mesh->heatmap(v);
		pc->normal(i) = mesh->normal(v);
		pc->rgb(i) = mesh->rgb(v);
	}

	spc.primID_splat = new unsigned int[size(primID_splat)];
	memcpy(spc.primID_splat, primID_splat.data(), sizeof(unsigned int) * size(primID_splat));

	spc.morton_codes = new unsigned int[size(points)];

	#pragma omp parallel for
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		auto & s = spc.splats[i];

		s.begin = idx_splats[i];
		s.end = idx_splats[i + 1];
		for(index_t j = s.begin; j < s.end; ++j)
			spc.morton_codes[j] = s.morton2d(points[j]);
	}

	splats_pcs.emplace_back(std::move(spc));

	return pc;
}

void splat::display_sets(che * pc, const std::vector<index_t> & sets, const index_t * mapid)
{
	std::vector<int> color(size(sets));
	std::iota(begin(color), end(color), 0);

	std::random_device rd;
	std::mt19937 gen{rd()};
	std::shuffle(begin(color), end(color), gen);

	for(index_t i = 1; i < size(sets); ++i)
	for(index_t j = sets[i - 1]; j < sets[i]; ++j)
		pc->heatmap(mapid ? mapid[j] : j) = real_t(color[i]) / size(color);
}

void splat::save_histogram(const std::string & file) const
{
	gproshan_error_var(file);

	FILE * fp = fopen(file.c_str(), "a");

	const splats_data & spc = splats_pcs.back();
	for(index_t i = 0; i < spc.n_splats; ++i)
	{
		auto & s = spc.splats[i];
		fprintf(fp, "%u %u\n", i, s.end - s.begin);
	}

	fclose(fp);
}


} // namespace gproshan

