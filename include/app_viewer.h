#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include "include.h"
#include "viewer/viewer.h"

#include "mesh/che_off.h"
#include "mesh/che_obj.h"
#include "mesh/che_ply.h"
#include "mesh/che_ptx.h"
#include "mesh/che_img.h"
#include "mesh/che_sphere.h"
#include "mesh/che_fill_hole.h"
#include "mesh/che_poisson.h"
#include "mesh/simplification.h"

#include "laplacian/laplacian.h"
#include "laplacian/fairing_taubin.h"
#include "laplacian/fairing_spectral.h"

#include "geodesics/dijkstra.h"
#include "geodesics/geodesics.h"
#include "geodesics/sampling.h"

#include "mdict/denoising.h"
#include "mdict/super_resolution.h"
#include "mdict/inpainting.h"
#include "mdict/basis_dct.h"
#include "mdict/patch.h"

#include "features/key_points.h"
#include "features/key_components.h"


// geometry processing and shape analysis framework
namespace gproshan {


class app_viewer : public viewer
{
	private:
		double time;

	public:
		app_viewer();
		~app_viewer();
		
		int main(int nargs, const char ** args);
	
	private:
		static bool process_fairing_taubin(viewer * p_view);
		static bool process_fairing_spectral(viewer * p_view);

		static bool process_geodesics(viewer * p_view, const geodesics::algorithm & alg);
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

