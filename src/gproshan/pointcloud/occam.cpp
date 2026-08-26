#include <gproshan/pointcloud/occam.h>

#include <gproshan/scenes/scanner.h>
#include <gproshan/pointcloud/knn.h>

#include <armadillo>
#include <map>


// geometry processing and shape analysis framework
namespace gproshan {



const std::vector<std::string> occam::rt_opts_str = { "mesh"
													, "area_npoints"
													, "mean_median_knn_distant_8"
													, "knn_area_per_point_1.5"
													, "knn_median_pairs_per_point_8"
													, "voronoi_8"
													, "voronoi_8_anisotropy"
													};

const size_t occam::n_tracers = size(rt_opts_str);

const std::vector<std::string> occam::patterns_str = {	"center"
														, "min"
														, "max"
														, "centerMin"
														, "centerMax"
														, "centerMinMax"
														};

const size_t occam::n_patterns = size(patterns_str);


bool occam::volume_dynamic_rays = false;
bool occam::use_inside_ray = true;


occam::data::data(const che * m): mesh(m)
{
	if(!mesh) return;

	get_min_max_vertex(mesh, min_vertex, max_vertex);

	vertex center = (min_vertex + max_vertex) / 2;
	vertex min_mid = (min_vertex + center) / 2;
	vertex max_mid = (max_vertex + center) / 2;

	patterns[0] = {center};
	patterns[1] = {min_mid};
	patterns[2] = {max_mid};
	patterns[3] = {min_mid, center};
	patterns[4] = {max_mid, center};
	patterns[5] = {min_mid, max_mid, center};

	vertex box = max_vertex - min_vertex;
	box_min = std::min({box.x(), box.y(), box.z()});
	area = 2 * (box.x() * box.y() + box.y() * box.z() + box.x() * box.z());
	volume = box.x() * box.y() * box.z();
	radius = sqrt(area / m->n_vertices);
}

rt::embree::pc_opts occam::rt_opts(const data & scene, const unsigned id)
{
	rt::embree::pc_opts pc_opts;
	pc_opts.enable = true;

	switch(id)
	{
		case 0: pc_opts.enable = false;
				break;

		case 1: pc_opts.opt = rt::embree::NONE;
				pc_opts.radius = scene.radius;
				break;

		case 2: pc_opts.opt = rt::embree::NONE;
				pc_opts.radius = knn::mean_median_knn_distant(&scene.mesh->point(0), scene.mesh->n_vertices, 8);
				break;

		case 3: pc_opts.opt = rt::embree::AREA;
				pc_opts.scale = 1.5;
				pc_opts.knn = 64;
				break;

		case 4: pc_opts.opt = rt::embree::MEDIAN_PAIRS;
				pc_opts.knn = 8;
				break;

		case 5: pc_opts.opt = rt::embree::VORONOI;
				pc_opts.knn = 8;
				break;

		case 6: pc_opts.opt = rt::embree::VORONOI;
				pc_opts.anisotropy = true;
				pc_opts.knn = 8;
				break;
	};

	return pc_opts;
}

void occam::get_min_max_vertex(const che * mesh, vertex & min_vertex, vertex & max_vertex)
{
	min_vertex = INFINITY;
	max_vertex = -INFINITY;
	for(unsigned v = 0; v < mesh->n_vertices; ++v)
	{
		const vertex & p = mesh->point(v);

		min_vertex.x() = std::min(min_vertex.x(), p.x());
		min_vertex.y() = std::min(min_vertex.y(), p.y());
		min_vertex.z() = std::min(min_vertex.z(), p.z());

		max_vertex.x() = std::max(max_vertex.x(), p.x());
		max_vertex.y() = std::max(max_vertex.y(), p.y());
		max_vertex.z() = std::max(max_vertex.z(), p.z());
	}
}

float occam::halton(int index, const int base)
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

std::vector<occam::pn_sample> occam::halton_sample_trigs( partitions & trig_samples
																, const che * mesh
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
		const vertex & a = mesh->point(vtrig(he));
		const vertex & b = mesh->point(vtrig(he + 1));
		const vertex & c = mesh->point(vtrig(he + 2));

		total_area += length(cross(a - c, b - c)) / 2;
	}

	const float samples_per_unit_area = total_samples / total_area;

	size_t idx = 0;
	for(unsigned t = 0; t < n_trigs; ++t)
	{
		const size_t & he = t * 3;
		const vertex & a = mesh->point(vtrig(he));
		const vertex & b = mesh->point(vtrig(he + 1));
		const vertex & c = mesh->point(vtrig(he + 2));

		const vertex & an = mesh->normal(vtrig(he));
		const vertex & bn = mesh->normal(vtrig(he + 1));
		const vertex & cn = mesh->normal(vtrig(he + 2));

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

			vertex p = alpha * a + beta * b + gamma * c;
			vertex n = alpha * an + beta * bn + gamma * cn;
			samples.emplace_back(p, n);
			++idx;
		}
	}

	return samples;
}

std::vector<int> occam::raycast_vpoint(	const rt::raytracing * rt
											, const std::vector<occam::pn_sample> & samples
											, const std::vector<vertex> & pviews
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
		for(const vertex & vpoint: pviews)
		{
			const float dist = length(vpoint - s.p);
			const vec3 & dir = (vpoint - s.p) / dist;

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

che * occam::scan(	const rt::raytracing * rt
							, const std::vector<vertex> & vo
							, const size_t rows
							, const size_t cols
							)
{
	che * out = scanner_ptx(rt, rows, cols, vo[0], true);
	for(size_t i = 1; i < size(vo); ++i)
	{
		che * p = scanner_ptx(rt, rows, cols, vo[i], true);
		che * q = out->merge(p);
		delete out;
		delete p;

		out = q;
	}

	return out;
}

float occam::inside_ray(const rt::raytracing * rt, const vertex & org, const int n_inside_ray)
{
	if(!occam::use_inside_ray)
		return 1;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0, 1);

	int hits = 0;
	for(int i = 0; i < n_inside_ray; ++i)
	{
		const float theta = dis(gen) * 2.f * M_PI;
		const float phi = acosf(2.f * dis(gen) - 1.f);
		const vertex dir = {sinf(phi) * cosf(theta), sinf(phi) * sinf(theta), cosf(phi)};
		if(rt->intersect(org, dir).primID != NIL)
			++hits;
	}

	return float(hits) / n_inside_ray;
}

int occam::generate_random_rays(std::vector<vertex> & origins
									, std::vector<vertex> & directions
									, const vertex & min_vertex
									, const vertex & max_vertex
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

	if(occam::volume_dynamic_rays)
	{
		const vertex c = max_vertex - min_vertex;
		num_rays *= c.x() * c.y() * c.z();
	}

	origins.clear();
	origins.reserve(num_rays);

	directions.clear();
	directions.reserve(num_rays);

	for(int i = 0; i < num_rays; ++i)
	{
		vertex org = {dis_x(gen), dis_y(gen), dis_z(gen)};
		vertex look_at = {dis_x(gen), dis_y(gen), dis_z(gen)};

		origins.push_back(org);
		directions.push_back(normalize(look_at - org));
	}

	return num_rays;
}

vec2 occam::raycast_random(	const rt::raytracing * rt
									, const vertex & min_vertex
									, const vertex & max_vertex
									, int num_rays
									, che ** out
									)
{
	std::vector<vertex> origins;
	std::vector<vertex> directions;

	num_rays = generate_random_rays(origins, directions, min_vertex, max_vertex, num_rays);

	float nohits = 0;
	float rays = 0;

	std::vector<float> p(num_rays);
	#pragma omp parallel for
	for(int i = 0; i < num_rays; ++i)
		p[i] = inside_ray(rt, origins[i]);

	std::vector<bool> ho(num_rays, 0);
//	gproshan_log_var(hist(p));

	#pragma omp parallel for
	for(int i = 0; i < num_rays; ++i)
	{
		const auto & hit = rt->intersect(origins[i], directions[i]);
		const float w = p[i];

		#pragma omp atomic
		rays += w;

		if((ho[i] = hit.primID == NIL))
		{
			#pragma omp atomic
			nohits += w;
		}
	}

	if(!out) return {nohits, rays};


	std::vector<int> inv;
	inv.reserve(num_rays);

	for(unsigned i = 0; i < size(ho); ++i)
		if(ho[i]) inv.push_back(i);

	gproshan_error_var(size(inv));

	delete out[0];
	out[0] = new che(size(inv));

	che & pc = *out[0];

	const int bin_res = 2;
	std::map<ivec3, float> bins;

	#pragma omp parallel for
	for(int v = 0; v < pc.n_vertices; ++v)
	{
		const int i = inv[v];
		pc.point(v) = origins[i];
		pc.heatmap(v) = p[i]; 

		auto o = origins[i] * bin_res;
		ivec3 b = {int(o.x()), int(o.y()), int(o.z())};

		#pragma omp critical
		bins[b] += pc.heatmap(v);
	}

	gproshan_error_var(min_vertex);
	gproshan_error_var(max_vertex);
	gproshan_log_var(size(bins));


	delete out[1];
	out[1] = new che(size(bins));

	float max_w = 0;
	for(const auto & p: bins)
		max_w = std::max(max_w, p.second);

	int i = 0;
	for(const auto & p: bins)
	{
		auto b = p.first;
		vertex v = {float(b.x()), float(b.y()), float(b.z())};
		out[1]->point(i) = (v + 0.5f) / bin_res;
		out[1]->heatmap(i) = p.second / max_w;
		++i;
	}

	return {nohits, rays};
}

float occam::occam_random(	const rt::raytracing * rt
									, const vertex & min_vertex
									, const vertex & max_vertex
									, const int num_rays
									, che ** out
									)
{
	auto occ = raycast_random(rt, min_vertex, max_vertex, num_rays, out);
	return occam_random(occ.x() / occ.y());
}

float occam::occam_random(const float ratio)
{
	return std::pow(ratio, 2.0 / 3.0);
}


int occam::main_test_bbr(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occam::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};
	const std::vector<float> vradius = {};//0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04;
//	bool scan_patterns = false;
	const size_t n_tests = 1;
//	const size_t n_rows = 1000;
//	const size_t n_cols = 1000;


	rt::embree::pc_opts pc_opts;

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

		che * pc = che::load_mesh(file);
		const occam::data scene(pc);

		for(auto & os: results)
			os << pc->name();

		for(unsigned i = 0; i < n_tracers; ++i)
		{
			if(!pc->n_trigs && !i) continue;

			rt::embree rt({pc}, {mat4::identity()}, rt_opts(scene, i));

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

				float sum_occam = 0;
				for(unsigned t = 0; t < n_tests; ++t)
					sum_occam += occam_random(&rt, scene.min_vertex, scene.max_vertex, num_rays);
				sum_occam /= n_tests;

				std::ofstream & os = results[i];

				for(float v: stats)
					os << " " << v;

				os << " " << sum_occam;
			}

/*
			if(!scan_patterns) continue;

			std::string filename = file + "_scan_" + std::to_string(i);
			for(int i = 0; i < 1; ++i)
			{
				che *out = scan(&rt, scene.patterns[i], n_rows, n_cols);
				out->heatmap_scale(pc->heatmap_scale());

				che_xyz::write_file(out, filename + "_" + std::to_string(i), true);

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

int occam::main_test_inside(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occam::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};


	rt::embree::pc_opts pc_opts;

	std::ifstream is(input);
	std::ofstream results[size(vnum_rays)][n_tracers];

	std::vector<vertex> orgs;
	std::vector<vertex> dirs;
	std::vector<vec<float,10>> inside;

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

		const che * pc = che::load_mesh(file);
		const occam::data scene(pc);

		for(unsigned t = 0; t < n_tracers; ++t)
		{
			if(!pc->n_trigs && !t) continue;

			rt::embree rt({pc}, {mat4::identity()}, rt_opts(scene, t));

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

int occam::main_test_intersection(const std::string & input)
{
	for(const auto & s: rt_opts_str)
		std::cout << s << "\n";

	if(input == "") return 0;


	// setup
	occam::volume_dynamic_rays = true;
	const std::vector<int> vnum_rays = {1000};


	rt::embree::pc_opts pc_opts;

	std::vector<std::ofstream> results(size(vnum_rays));
	for(unsigned i = 0; i < size(vnum_rays); ++i)
	{
		std::ofstream & os = results[i];
		os.open("hit_" + std::to_string(vnum_rays[i]) + "_" + input);
	}

	std::ifstream is(input);

	std::vector<vertex> orgs;
	std::vector<vertex> dirs;

	int miss_hit, extra_hit;


	std::string file;
	while(is >> file)
	{
		std::cerr << "processing: " << file << "\n";

		const che * pc = che::load_mesh(file);
		const rt::embree mrt({pc}, {mat4::identity()});
		const occam::data scene(pc);

		for(auto & os: results)
			os << pc->n_vertices << " " << pc->n_trigs;

		for(unsigned t = 1; t < n_tracers; ++t)
		{
			const rt::embree prt({pc}, {mat4::identity()}, rt_opts(scene, t));

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

int occam::main_test_reconstruction(const std::string & input)
{
	const std::string home = std::getenv("HOME");
	const std::string path = home + "/scannet++/";

	const int num_rays = 1000;
	occam::volume_dynamic_rays = true;

	std::string file;
	std::string filename;

	std::ifstream is(input);
	std::ofstream os("fscore_" + input);

	while(is >> file)
	{
		filename = path + file + "/scans/mesh_aligned_0.05.ply";
		gproshan_error_var(filename);

		che * orig = che::load_mesh(filename);
		che * scan = che::load_mesh(filename + "_scan_-2.000000_0.xyz");
		che * reco = che::load_mesh(filename + "_scan_-2.000000_0.xyz.off");	// reconstructed

		const occam::data scene(orig);
		const float d = 0.01 * scene.box_min;

		os << file << " " << f_score(&orig->point(0), orig->n_vertices, &reco->point(0), reco->n_vertices, d);

		for(unsigned i = 1; i < occam::n_tracers; ++i)
		{
			rt::embree rt({scan}, {mat4::identity()}, rt_opts(scene, i));

			os << " " << occam_random(&rt, scene.min_vertex, scene.max_vertex, num_rays);
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
	knn::k3tree nn(P, nP, Q, nQ, 1);

	int sum = 0;

	#pragma omp parallel for
	for(unsigned i = 0; i < nQ; ++i)
		if(length(Q[i] - P[nn(i, 0)]) < d)
		{
			#pragma omp atomic
			++sum;
		}

	return 100.f * sum / nQ;
}


} // namespace gproshan

