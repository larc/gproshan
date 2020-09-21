#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include "viewer/viewer.h"
#include "include.h"
#include "che_off.h"
#include "che_obj.h"
#include "che_ply.h"
#include "che_ptx.h"
#include "che_img.h"
#include "che_sphere.h"
#include "laplacian.h"
#include "che_off.h"
#include "dijkstra.h"
#include "geodesics.h"
#include "fairing_taubin.h"
#include "fairing_spectral.h"
#include "sampling.h"
#include "che_fill_hole.h"
#include "che_poisson.h"
#include "decimation.h"
#include "mdict/denoising.h"
#include "mdict/super_resolution.h"
#include "mdict/inpainting.h"
#include "mdict/basis_dct.h"
#include "mdict/patch.h"
#include "key_points.h"
#include "key_components.h"


// geometry processing and shape analysis framework
namespace gproshan {


class app_viewer : public viewer
{
	private:
		double time;
		real_t * dist;
		size_t n_dist;

	public:
		app_viewer();
		~app_viewer();
		
		int main(int nargs, const char ** args);
	
	private:
		static bool process_fairing_taubin(viewer * p_view);
		static bool process_fairing_spectral(viewer * p_view);

		static bool process_fastmarching(viewer * p_view);
		static bool process_geodesics_fm(viewer * p_view);
		static bool process_geodesics_ptp_cpu(viewer * p_view);
		static bool process_geodesics_heat_flow(viewer * p_view);

	#ifdef GPROSHAN_CUDA
		static bool process_geodesics_ptp_gpu(viewer * p_view);
		static bool process_geodesics_heat_flow_gpu(viewer * p_view);
	#endif // GPROSHAN_CUDA

		static bool process_farthest_point_sampling(viewer * p_view);
		static bool process_farthest_point_sampling_radio(viewer * p_view);
		static bool compute_toplesets(viewer * p_view);
		static bool process_voronoi(viewer * p_view);

		static bool process_mdict_patch(viewer * p_view);
		static bool process_denoising(viewer * p_view);
		static bool process_super_resolution(viewer * p_view);
		static bool process_inpaiting(viewer * p_view);
		static bool process_iterative_inpaiting(viewer * p_view);

		static bool process_functional_maps(viewer * p_view);
		static bool process_gps(viewer * p_view);
		static bool process_hks(viewer * p_view);
		static bool process_wks(viewer * p_view);
		static bool process_key_points(viewer * p_view);
		static bool process_key_components(viewer * p_view);

		static bool process_poisson(viewer * p_view, const index_t & k);
		static bool process_poisson_laplacian_1(viewer * p_view);
		static bool process_poisson_laplacian_2(viewer * p_view);
		static bool process_poisson_laplacian_3(viewer * p_view);

		static bool process_threshold(viewer * p_view);
		static bool process_noise(viewer * p_view);
		static bool process_black_noise(viewer * p_view);
		static bool process_multiplicate_vertices(viewer * p_view);
		static bool process_fill_holes(viewer * p_view);
		static bool process_delete_vertices(viewer * p_view);
		static bool process_fill_holes_test(viewer * p_view);
		static bool process_delete_non_manifold_vertices(viewer * p_view);
		static bool process_fill_holes_biharmonic_splines(viewer * p_view);
		static bool process_gaussian_curvature(viewer * p_view);
		static bool process_edge_collapse(viewer * p_view);
		static bool select_multiple(viewer * p_view);
};


} // namespace gproshan


#endif //APP_VIEWER_H

