#include <gproshan/pointcloud/occlusion.h>

#include <gproshan/scenes/scanner.h>
#include <gproshan/pointcloud/knn.h>

#include <armadillo>
#include <map>


const std::vector<std::string> occlusion::rt_opts_str = {	"mesh"
															, "area_npoints"
															, "mean_median_knn_distant_8"
															, "knn_area_per_point_1.5"
															, "knn_median_pairs_per_point_8"
															, "voronoi_8"
															, "voronoi_8_anisotropy"
															};

const size_t occlusion::n_tracers = size(rt_opts_str);

const std::vector<std::string> occlusion::patterns_str = {	"center",
															"min",
															"max",
															"centerMin",
															"centerMax",
															"centerMinMax"
															};

const size_t occlusion::n_patterns = size(patterns_str);


bool occlusion::volume_dynamic_rays = false;
bool occlusion::use_inside_ray = true;


occlusion::data::data(const gp::che * m): mesh(m)
{
	if(!mesh) return;

	get_min_max_vertex(mesh, min_vertex, max_vertex);

	gp::vertex center = (min_vertex + max_vertex) / 2;
	gp::vertex min_mid = (min_vertex + center) / 2;
	gp::vertex max_mid = (max_vertex + center) / 2;

	patterns[0] = {center};
	patterns[1] = {min_mid};
	patterns[2] = {max_mid};
	patterns[3] = {min_mid, center};
	patterns[4] = {max_mid, center};
	patterns[5] = {min_mid, max_mid, center};

	gp::vertex box = max_vertex - min_vertex;
	box_min = std::min({box.x(), box.y(), box.z()});
	area = 2 * (box.x() * box.y() + box.y() * box.z() + box.x() * box.z());
	volume = box.x() * box.y() * box.z();
	radius = sqrt(area / m->n_vertices);
}

gp::rt::embree::pc_opts occlusion::rt_opts(const data & scene, const unsigned id)
{
	gp::rt::embree::pc_opts pc_opts;
	pc_opts.enable = true;

	switch(id)
	{
		case 0: pc_opts.enable = false;
				break;

		case 1: pc_opts.opt = gp::rt::embree::NONE;
				pc_opts.radius = scene.radius;
				break;

		case 2: pc_opts.opt = gp::rt::embree::NONE;
				pc_opts.radius = gp::knn::mean_median_knn_distant(&scene.mesh->point(0), scene.mesh->n_vertices, 8);
				break;

		case 3: pc_opts.opt = gp::rt::embree::AREA;
				pc_opts.scale = 1.5;
				pc_opts.knn = 64;
				break;

		case 4: pc_opts.opt = gp::rt::embree::MEDIAN_PAIRS;
				pc_opts.knn = 8;
				break;

		case 5: pc_opts.opt = gp::rt::embree::VORONOI;
				pc_opts.knn = 8;
				break;

		case 6: pc_opts.opt = gp::rt::embree::VORONOI;
				pc_opts.anisotropy = true;
				pc_opts.knn = 8;
				break;
	};

	return pc_opts;
}

void occlusion::get_min_max_vertex(const gp::che * mesh, gp::vertex & min_vertex, gp::vertex & max_vertex)
{
	min_vertex = INFINITY;
	max_vertex = -INFINITY;
	for(unsigned v = 0; v < mesh->n_vertices; ++v)
	{
		const gp::vertex & p = mesh->point(v);

		min_vertex.x() = std::min(min_vertex.x(), p.x());
		min_vertex.y() = std::min(min_vertex.y(), p.y());
		min_vertex.z() = std::min(min_vertex.z(), p.z());

		max_vertex.x() = std::max(max_vertex.x(), p.x());
		max_vertex.y() = std::max(max_vertex.y(), p.y());
		max_vertex.z() = std::max(max_vertex.z(), p.z());
	}
}

float occlusion::halton(int index, const int base)
{
	float result = 0.0;
	float f = 1.0 / base;

	while(index > 0)
	{
		result += f * (index % base);
		index /= base;
		f /= base;
	}

	return result;
}

std::vector<occlusion::pn_sample> occlusion::halton_sample_trigs( gp::partitions & trig_samples
																, const gp::che * mesh
																, const int total_samples
																)
{
	std::vector<pn_sample> samples;

	unsigned n_trigs = mesh->is_scene() ? mesh->n_vertices / 3 : mesh->n_trigs;

	auto vtrig = [&](const unsigned i)
	{
		return mesh->is_scene() ? i : mesh->halfedge(i);
	};

	float total_area = 0;

	#pragma omp parallel for reduction(+: total_area)
	for(unsigned t = 0; t < n_trigs; ++t)
	{
		const size_t & he = t * 3;
		const gp::vertex & a = mesh->point(vtrig(he));
		const gp::vertex & b = mesh->point(vtrig(he + 1));
		const gp::vertex & c = mesh->point(vtrig(he + 2));

		total_area += length(cross(a - c, b - c)) / 2;
	}

	const float samples_per_unit_area = total_samples / total_area;

	size_t idx = 0;
	for(unsigned t = 0; t < n_trigs; ++t)
	{
		const size_t & he = t * 3;
		const gp::vertex & a = mesh->point(vtrig(he));
		const gp::vertex & b = mesh->point(vtrig(he + 1));
		const gp::vertex & c = mesh->point(vtrig(he + 2));

		const gp::vertex & an = mesh->normal(vtrig(he));
		const gp::vertex & bn = mesh->normal(vtrig(he + 1));
		const gp::vertex & cn = mesh->normal(vtrig(he + 2));

		const float area = length(cross(a - c, b - c)) / 2;
		size_t num_samples = area * samples_per_unit_area + 0.5f;

		trig_samples.add(num_samples);

		while(num_samples--)
		{
			float r1 = halton(idx + 1, 2); // use 2 as base
			float r2 = halton(idx + 1, 3);

			float sqrtR1 = std::sqrt(r1);
			float alpha = 1 - sqrtR1;
			float beta = r2 * sqrtR1;
			float gamma = 1 - alpha - beta;

			gp::vertex p = alpha * a + beta * b + gamma * c;
			gp::vertex n = alpha * an + beta * bn + gamma * cn;
			samples.emplace_back(p, n);
			++idx;
		}
	}

	return samples;
}

std::vector<int> occlusion::raycast_vpoint(	const gp::rt::raytracing * rt
											, const std::vector<occlusion::pn_sample> & samples
											, const std::vector<gp::vertex> & pviews
											, const float min_dist
											)
{
	std::vector<int> vocc(size(samples));

	#pragma omp parallel for
	for(unsigned i = 0; i < size(samples); ++i)
	{
		const auto & s = samples[i];

		unsigned visible = 0;
		unsigned occluded = 0;
		for(const gp::vertex & vpoint: pviews)
		{
			const float dist = length(vpoint - s.p);
			const gp::vec3 & dir = (vpoint - s.p) / dist;

			if(dot(dir, s.n) < 0) continue;

			const auto & hit = rt->intersect(s.p + 1e-04 * dir, dir);

			if(hit.dist < min_dist) continue;

			++visible;
			if(hit.dist < dist)
				++occluded;
		}

		vocc[i] = (visible > 0) + (occluded == size(pviews));
	}

	return vocc;
}

gp::che * occlusion::scan(	const gp::rt::raytracing * rt
							, const std::vector<gp::vertex> & vo
							, const size_t rows
							, const size_t cols
							)
{
	gp::che * out = gp::scanner_ptx(rt, rows, cols, vo[0], true);
	for(size_t i = 1; i < size(vo); ++i)
	{
		gp::che * p = gp::scanner_ptx(rt, rows, cols, vo[i], true);
		gp::che * q = out->merge(p);
		delete out;
		delete p;

		out = q;
	}

	return out;
}

float occlusion::inside_ray(const gp::rt::raytracing * rt, const gp::vertex & org, const int n_inside_ray)
{
	if(!occlusion::use_inside_ray)
		return 1;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0, 1);

	int hits = 0;
	for(int i = 0; i < n_inside_ray; ++i)
	{
		const float theta = dis(gen) * 2.f * M_PI;
		const float phi = acosf(2.f * dis(gen) - 1.f);
		const gp::vertex dir = {sinf(phi) * cosf(theta), sinf(phi) * sinf(theta), cosf(phi)};
		if(rt->intersect(org, dir).primID != NIL)
			++hits;
	}

	return float(hits) / n_inside_ray;
}

int occlusion::generate_random_rays(std::vector<gp::vertex> & origins
									, std::vector<gp::vertex> & directions
									, const gp::vertex & min_vertex
									, const gp::vertex & max_vertex
									, int num_rays
									)
{
	if(	min_vertex.x() > max_vertex.x()
		|| min_vertex.y() > max_vertex.y()
		|| min_vertex.z() > max_vertex.z())
	{
		std::cerr << "Invalid min_pt and max_pt values." << std::endl;
		return 0;
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis_x(min_vertex.x(), max_vertex.x());
	std::uniform_real_distribution<float> dis_y(min_vertex.y(), max_vertex.y());
	std::uniform_real_distribution<float> dis_z(min_vertex.z(), max_vertex.z());

	if(occlusion::volume_dynamic_rays)
	{
		const gp::vertex c = max_vertex - min_vertex;
		num_rays *= c.x() * c.y() * c.z();
	}

	origins.clear();
	origins.reserve(num_rays);

	directions.clear();
	directions.reserve(num_rays);

	for(int i = 0; i < num_rays; ++i)
	{
		gp::vertex org = {dis_x(gen), dis_y(gen), dis_z(gen)};
		gp::vertex look_at = {dis_x(gen), dis_y(gen), dis_z(gen)};

		origins.push_back(org);
		directions.push_back(normalize(look_at - org));
	}

	return num_rays;
}

gp::vec2 occlusion::raycast_random(	const gp::rt::raytracing * rt
									, const gp::vertex & min_vertex
									, const gp::vertex & max_vertex
									, int num_rays
									, gp::che ** out
									)
{
	std::vector<gp::vertex> origins;
	std::vector<gp::vertex> directions;

	num_rays = generate_random_rays(origins, directions, min_vertex, max_vertex, num_rays);

	const int bin_res = 2;
	std::map<gp::ivec3, float> bins;

	float nohits = 0;
	float rays = 0;

	if(out)
	{
		delete out[0];
		out[0] = new gp::che(num_rays);
	}

	arma::frowvec p(num_rays);
	#pragma omp parallel for
	for(int i = 0; i < num_rays; ++i)
		p[i] = inside_ray(rt, origins[i]);

//	gproshan_log_var(hist(p));

	#pragma omp parallel for
	for(int i = 0; i < num_rays; ++i)
	{
		const auto & hit = rt->intersect(origins[i], directions[i]);
		const float w = p[i];

		#pragma omp atomic
		rays += w;

		if(hit.primID == NIL)
		{
			#pragma omp atomic
			nohits += w;
		}

		if(out && hit.primID == NIL)
		{
			out[0]->point(i) = origins[i];
			out[0]->heatmap(i) = w;
			out[0]->rgb(i) = gp::vertex{0,w,0};
		}
		else if(out)
		{
			out[0]->heatmap(i) = -1;
		}
	}

	if(out)
	{
		const gp::che & pc = *out[0];
		for(int i = 0; i < num_rays; ++i)
			if(pc.heatmap(i) >= 0)
			{
				auto v = origins[i] * bin_res;
				gp::ivec3 b = {v.x(), v.y(), v.z()};
				bins[b] += pc.heatmap(i);
			}
	}

	gproshan_error_var(min_vertex);
	gproshan_error_var(max_vertex);
	gproshan_log_var(bins.size());

	if(out)
	{
		delete out[1];
		out[1] = new gp::che(bins.size());
	}

	float max_w = 0;
	for(const auto & p: bins)
		max_w = std::max(max_w, p.second);

	int i = 0;
	for(const auto & p: bins)
	{
		auto b = p.first;
		gp::vertex v = {b.x(), b.y(), b.z()};
		out[1]->point(i) = (v + 0.5f) / bin_res;
		out[1]->heatmap(i) = p.second / max_w;
		++i;
	}

	return {nohits, rays};
}

float occlusion::occlusion_random(	const gp::rt::raytracing * rt
									, const gp::vertex & min_vertex
									, const gp::vertex & max_vertex
									, const int num_rays
									, gp::che ** out
									)
{
	auto occ = raycast_random(rt, min_vertex, max_vertex, num_rays, out);
	return occlusion_random(occ.x() / occ.y());
}

float occlusion::occlusion_random(const float ratio)
{
	return std::pow(ratio, 2.0 / 3.0);
}


int occlusion::main_test_bbr(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occlusion::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};
	const std::vector<float> vradius = {};//0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04;
//	bool scan_patterns = false;
	const size_t n_tests = 1;
//	const size_t n_rows = 1000;
//	const size_t n_cols = 1000;


	gp::rt::embree::pc_opts pc_opts;

	std::vector<std::ofstream> results(size(vnum_rays));
	for(unsigned i = 0; i < size(vnum_rays); ++i)
	{
		std::ofstream & os = results[i];
		os.open("bbr_" + std::to_string(vnum_rays[i]) + "_" + input);
	}

	std::vector<float> stats;	// radius | min max mean median stddev
	std::ifstream is(input);

	std::string file;
	while(is >> file)
	{
		std::cerr << "processing: " << file << "\n";

		gp::che * pc = gp::che::load_mesh(file);
		const occlusion::data scene(pc);

		for(auto & os: results)
			os << pc->name();

		for(unsigned i = 0; i < n_tracers; ++i)
		{
			if(!pc->n_trigs && !i) continue;

			gp::rt::embree rt({pc}, {gp::mat4::identity()}, rt_opts(scene, i));

/*
			stats.clear();
			if(i == 2 || i == 3) // radius stats
			{
				const arma::fmat xyzr((float *) rt.pc_data(), 4, pc->n_vertices, false, true);
				const auto & vradii = xyzr.row(3);

				stats.push_back(min(vradii));
				stats.push_back(max(vradii));
				stats.push_back(mean(vradii));
				stats.push_back(median(vradii));
				stats.push_back(stddev(vradii));
			}
			else if(i) stats.push_back(pc_opts.radius);
*/

			for(unsigned i = 0; i < size(vnum_rays); ++i)
			{
				const int num_rays = vnum_rays[i];

				float sum_occlusion = 0;
				for(unsigned t = 0; t < n_tests; ++t)
					sum_occlusion += occlusion_random(&rt, scene.min_vertex, scene.max_vertex, num_rays);
				sum_occlusion /= n_tests;

				std::ofstream & os = results[i];

				for(float v: stats)
					os << " " << v;

				os << " " << sum_occlusion;
			}

/*
			if(!scan_patterns) continue;

			std::string filename = file + "_scan_" + std::to_string(i);
			for(int i = 0; i < 1; ++i)
			{
				gp::che *out = scan(&rt, scene.patterns[i], n_rows, n_cols);
				out->heatmap_scale(pc->heatmap_scale());

				gp::che_xyz::write_file(out, filename + "_" + std::to_string(i), true);

				delete out;
			}

			scan_patterns = false;
*/
		}

		for(auto & os: results)
			os << "\n";

		delete pc;

//		scan_patterns = true;
	}

	for(auto & os: results)
		os.close();

	return 0;
}

int occlusion::main_test_inside(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occlusion::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};


	gp::rt::embree::pc_opts pc_opts;

	std::ifstream is(input);
	std::ofstream results[size(vnum_rays)][n_tracers];

	std::vector<gp::vertex> orgs;
	std::vector<gp::vertex> dirs;
	std::vector<gp::vec<float,10>> inside;

	std::string file;
	while(is >> file)
	{
		std::cerr << "processing: " << file << "\n";

		for(unsigned i = 0; i < size(vnum_rays); ++i)
		for(unsigned t = 0; t < n_tracers; ++t)
		{
			std::ofstream & os = results[i][t];
			os.open(file + "_" + std::to_string(vnum_rays[i]) + "_inside_" + std::to_string(t + 1));
		}

		const gp::che * pc = gp::che::load_mesh(file);
		const occlusion::data scene(pc);

		for(unsigned t = 0; t < n_tracers; ++t)
		{
			if(!pc->n_trigs && !t) continue;

			gp::rt::embree rt({pc}, {gp::mat4::identity()}, rt_opts(scene, t));

			for(unsigned i = 0; i < size(vnum_rays); ++i)
			{
				const int num_rays = generate_random_rays(orgs, dirs, scene.min_vertex, scene.max_vertex, vnum_rays[i]);
				std::ofstream & os = results[i][t];

				inside.resize(num_rays);

				#pragma omp parallel for
				for(int j = 0; j < num_rays; ++j)
				for(int k = 0; k < 10; ++k)
				{
					inside[j][k] = inside_ray(&rt, orgs[j],  (k + 1) * 10);
				}

				for(const auto & p: inside)
				{
					os << p[0];
					for(int k = 1; k < 10; ++k)
						os << " " << p[k];
					os << "\n";
				}
			}
		}

		delete pc;

		for(unsigned i = 0; i < size(vnum_rays); ++i)
		for(unsigned t = 0; t < n_tracers; ++t)
			results[i][t].close();
	}

	return 0;
}

int occlusion::main_test_intersection(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occlusion::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};


	gp::rt::embree::pc_opts pc_opts;

	std::vector<std::ofstream> results(size(vnum_rays));
	for(unsigned i = 0; i < size(vnum_rays); ++i)
	{
		std::ofstream & os = results[i];
		os.open("hit_" + std::to_string(vnum_rays[i]) + "_" + input);
	}

	std::ifstream is(input);

	std::vector<gp::vertex> orgs;
	std::vector<gp::vertex> dirs;

	int miss_hit, extra_hit;


	std::string file;
	while(is >> file)
	{
		std::cerr << "processing: " << file << "\n";

		const gp::che * pc = gp::che::load_mesh(file);
		const gp::rt::embree mrt({pc}, {gp::mat4::identity()});
		const occlusion::data scene(pc);

		for(auto & os: results)
			os << pc->n_vertices << " " << pc->n_trigs;

		for(unsigned t = 1; t < n_tracers; ++t)
		{
			const gp::rt::embree prt({pc}, {gp::mat4::identity()}, rt_opts(scene, t));

			for(unsigned i = 0; i < size(vnum_rays); ++i)
			{
				const int num_rays = generate_random_rays(orgs, dirs, scene.min_vertex, scene.max_vertex, vnum_rays[i]);

				miss_hit = extra_hit = 0;

				#pragma omp parallel for
				for(int j = 0; j < num_rays; ++j)
				{
					const auto & mhit = mrt.intersect(orgs[j], dirs[j]);
					const auto & phit = prt.intersect(orgs[j], dirs[j]);

					if(mhit.primID != NIL && phit.primID == NIL)
					{
						#pragma omp atomic
						++miss_hit;
					}
					if(mhit.primID == NIL && phit.primID != NIL)
					{
						#pragma omp atomic
						++extra_hit;
					}
				}

				std::ofstream & os = results[i];

				if(t == 1) os << " " << num_rays;

				os << " " << 100000 * float(miss_hit) / num_rays;
				os << " " << 100000 * float(extra_hit) / num_rays;
			}
		}

		for(auto & os: results)
			os << "\n";

		delete pc;
	}

	for(auto & os: results)
		os.close();

	return 0;
}

int occlusion::main_test_reconstruction(const std::string & input)
{
	const std::string home = std::getenv("HOME");
	const std::string path = home + "/scannet++/";

	const int num_rays = 1000;
	occlusion::volume_dynamic_rays = true;

	std::string file;
	std::string filename;

	std::ifstream is(input);
	std::ofstream os("fscore_" + input);

	while(is >> file)
	{
		filename = path + file + "/scans/mesh_aligned_0.05.ply";
		gproshan_error_var(filename);

		gp::che * orig = gp::che::load_mesh(filename);
		gp::che * scan = gp::che::load_mesh(filename + "_scan_-2.000000_0.xyz");
		gp::che * reco = gp::che::load_mesh(filename + "_scan_-2.000000_0.xyz.off");	// reconstructed

		const occlusion::data scene(orig);
		const float d = 0.01 * scene.box_min;

		os << file << " " << f_score(&orig->point(0), orig->n_vertices, &reco->point(0), reco->n_vertices, d);

		for(unsigned i = 1; i < occlusion::n_tracers; ++i)
		{
			gp::rt::embree rt({scan}, {gp::mat4::identity()}, rt_opts(scene, i));

			os << " " << occlusion_random(&rt, scene.min_vertex, scene.max_vertex, num_rays);
		}

		os << "\n";

		delete orig;
		delete scan;
		delete reco;
	}

	is.close();
	os.close();

	return 0;
}


// https://dl.acm.org/doi/10.1145/3072959.3073599
float f_score(const point * G, const size_t nG, const point * R, const size_t nR, const float d)
{
	float p = f_score_percent(R, nR, G, nG, d);	// precision
	float r = f_score_percent(G, nG, R, nR, d);	// recall

	return 2 * p * r / (p + r);
}

float f_score_percent(const point * Q, const size_t nQ, const point * P, const size_t nP, const float d)
{
	gp::knn::k3tree nn(P, nP, Q, nQ, 1);

	int sum = 0;

	#pragma omp parallel for
	for(unsigned i = 0; i < nQ; ++i)
		if(gp::length(Q[i] - P[nn(i, 0)]) < d)
		{
			#pragma omp atomic
			++sum;
		}

	return 100.f * sum / nQ;
}

