#ifndef OCCAM_H
#define OCCAM_H

#include <gproshan/app_viewer.h>
#include <gproshan/raytracing/embree.h>


namespace gp = gproshan;

using point = gp::vec3;

class occam
{
	public:
		struct data
		{
			const gp::che * mesh = nullptr;
			gp::vertex min_vertex;
			gp::vertex max_vertex;
			std::vector<gp::vertex> patterns[6];
			float box_min = 0;
			float area = 0;
			float volume = 0;
			float radius = 0;

			data(const gp::che * m = nullptr);
		};

		struct pn_sample
		{
			gp::vertex p;
			gp::vertex n;
		};

		static const std::vector<std::string> rt_opts_str;
		static const size_t n_tracers;

		static const std::vector<std::string> patterns_str;
		static const size_t n_patterns;

		static bool volume_dynamic_rays;
		static bool use_inside_ray;

	private:
		data scene;
		std::vector<gp::vertex> viewpoints;

	public:
		static gp::rt::embree::pc_opts rt_opts(const data & scene, const unsigned id);

		static void get_min_max_vertex(const gp::che * mesh, gp::vertex & min_vertex, gp::vertex & max_vertex);

		static float halton(const int index, const int base);

		static std::vector<pn_sample> halton_sample_trigs(	gp::partitions & trig_samples
															, const gp::che * mesh
															, const int total_samples
															);

		// return ints: 0 not visible, 1 is visible at least for 1 view point, 2 is occluded for all viewpoints
		static std::vector<int> raycast_vpoint( const gp::rt::raytracing * rt
												, const std::vector<pn_sample> & samples
												, const std::vector<gp::vertex> & pviews
												, const float min_dist
												);

		static gp::che * scan(	const gp::rt::raytracing * rt
								, const std::vector<gp::vertex> & vo
								, const size_t rows
								, const size_t cols
								);

		static float inside_ray(const gp::rt::raytracing * rt, const gp::vertex & org, const int n_inside_ray = 100);

		static int generate_random_rays(std::vector<gp::vertex> & origins
										, std::vector<gp::vertex> & destinations
										, const gp::vertex & min_vertex
										, const gp::vertex & max_vertex
										, int num_rays
										);

		static gp::vec2 raycast_random(	const gp::rt::raytracing * rt
										, const gp::vertex & min_vertex
										, const gp::vertex & max_vertex
										, const int num_rays
										, gp::che ** out = nullptr
										);

		static float occam_random(	const gp::rt::raytracing * rt
										, const gp::vertex & min_vertex
										, const gp::vertex & max_vertex
										, const int num_rays
										, gp::che ** out = nullptr
										);

		static float occam_random(const float ratio);

		static int main_test_bbr(const std::string & input);
		static int main_test_inside(const std::string & input);
		static int main_test_intersection(const std::string & input);
		static int main_test_reconstruction(const std::string & input);
};


float f_score(const point * G, const size_t nG, const point * R, const size_t nR, const float d);
float f_score_percent(const point * A, const size_t nA, const point * B, const size_t nB, const float d);


#endif // OCCLUSION_H

