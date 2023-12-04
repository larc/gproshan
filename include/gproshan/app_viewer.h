#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include <gproshan/include.h>
#include <gproshan/viewer/viewer.h>

#include <gproshan/mesh/che_obj.h>
#include <gproshan/mesh/che_off.h>
#include <gproshan/mesh/che_ply.h>
#include <gproshan/mesh/che_ptx.h>
#include <gproshan/mesh/che_xyz.h>
#include <gproshan/mesh/che_pts.h>
#include <gproshan/mesh/che_pcd.h>
#include <gproshan/mesh/che_img.h>
#include <gproshan/mesh/che_sphere.h>
#include <gproshan/mesh/che_fill_hole.h>
#include <gproshan/mesh/che_poisson.h>
#include <gproshan/mesh/simplification.h>

#include <gproshan/laplacian/laplacian.h>
#include <gproshan/laplacian/fairing_taubin.h>
#include <gproshan/laplacian/fairing_spectral.h>

#include <gproshan/scenes/scanner.h>

#include <gproshan/geometry/points.h>
#include <gproshan/geometry/convex_hull.h>

#include <gproshan/geodesics/dijkstra.h>
#include <gproshan/geodesics/geodesics.h>
#include <gproshan/geodesics/sampling.h>

#include <gproshan/mdict/msparse_coding.h>
#include <gproshan/mdict/basis_dct.h>
#include <gproshan/mdict/patch.h>

#include <gproshan/features/descriptor.h>
#include <gproshan/features/key_points.h>
#include <gproshan/features/key_components.h>


// geometry processing and shape analysis framework
namespace gproshan {


class app_viewer : public viewer
{
	protected:
		double time;

	public:
		app_viewer(const char * title = "gproshan", const int & width = 1600, const int & height = 900);
		virtual ~app_viewer();

		int main(int nargs, const char ** args);

	protected:
		virtual void init();

		// Point Cloud
		static bool process_knn(viewer * p_view);
		static bool process_compute_normals(viewer * p_view);

		// Scenes
		static bool process_simulate_scanner(viewer * p_view);
		static bool process_scatter(viewer * p_view);

		// Geometry
		static bool process_sampling_4points(viewer * p_view);
		static bool process_convex_hull(viewer * p_view);
		static bool process_connected_components(viewer * p_view);
		static bool process_gaussian_curvature(viewer * p_view);
		static bool process_edge_collapse(viewer * p_view);
		static bool process_multiplicate_vertices(viewer * p_view);
		static bool process_delete_vertices(viewer * p_view);
		static bool process_delete_non_manifold_vertices(viewer * p_view);

		// Fairing
		static bool process_fairing_taubin(viewer * p_view);
		static bool process_fairing_spectral(viewer * p_view);

		// Geodesics
		static bool process_geodesics(viewer * p_view);
		static bool process_farthest_point_sampling(viewer * p_view);
		static bool process_voronoi(viewer * p_view);
		static bool process_compute_toplesets(viewer * p_view);

		// Mesh Sparse Coding
		static bool process_msparse_coding(viewer * p_view);
		static bool process_mdict_patch(viewer * p_view);
		static bool process_mask(viewer * p_view);
		static bool process_pc_reconstruction(viewer * p_view);

		// Features
		static bool process_eigenfuntions(viewer * p_view);
		static bool process_descriptor_heatmap(viewer * p_view);
		static bool process_key_points(viewer * p_view);
		static bool process_key_components(viewer * p_view);

		// Hole Filling
		static bool process_poisson(viewer * p_view, const index_t & k);
		static bool process_poisson_laplacian_1(viewer * p_view);
		static bool process_poisson_laplacian_2(viewer * p_view);
		static bool process_poisson_laplacian_3(viewer * p_view);
		static bool process_fill_holes(viewer * p_view);
		static bool process_fill_holes_biharmonic_splines(viewer * p_view);

		// Others
		static bool process_select_multiple(viewer * p_view);
		static bool process_threshold(viewer * p_view);
		static bool process_noise(viewer * p_view);
		static bool process_black_noise(viewer * p_view);
};


} // namespace gproshan


#endif //APP_VIEWER_H

