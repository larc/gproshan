#ifndef APP_VIEWER_H
#define APP_VIEWER_H

#include "include.h"
#include "viewer/viewer.h"

#include "mesh/che_off.h"
#include "mesh/che_obj.h"
#include "mesh/che_ply.h"
#include "mesh/che_ptx.h"
#include "mesh/che_xyz.h"
#include "mesh/che_img.h"
#include "mesh/che_sphere.h"
#include "mesh/che_fill_hole.h"
#include "mesh/che_poisson.h"
#include "mesh/simplification.h"

#include "laplacian/laplacian.h"
#include "laplacian/fairing_taubin.h"
#include "laplacian/fairing_spectral.h"

#include "geometry/convex_hull.h"

#include "geodesics/dijkstra.h"
#include "geodesics/geodesics.h"
#include "geodesics/sampling.h"

#include "mdict/msparse_coding.h"
#include "mdict/basis_dct.h"
#include "mdict/patch.h"

#include "features/descriptor.h"
#include "features/key_points.h"
#include "features/key_components.h"


// geometry processing and shape analysis framework
namespace gproshan {


class app_viewer : public viewer
{
	protected:
		double time;

	public:
		app_viewer() = default;
		~app_viewer();

		int main(int nargs, const char ** args);

	protected:
		virtual void init();

		che * load_mesh(const string & file_path);

		// Geometry
		static bool process_convex_hull(viewer * p_view);
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

