#ifndef PATCH_H
#define PATCH_H

#include "include.h"
#include "che.h"
#include "include_arma.h"
#include "geodesics.h"

#include <vector>
#include <map>
#include <CImg.h>
#include <algorithm> 

#ifdef Success
	#undef Success
#endif

using namespace cimg_library;

using namespace std;


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class dictionary;

typedef function<bool(const index_t &)> fmask_t;
typedef function<bool(const index_t &, size_t tam)> fmask_local_t;
typedef std::map<index_t, index_t> vpatches_t;

/// 
class patch
{
	public:
		std::vector<index_t> vertices;		///< Vertices of the patch.
		a_mat T;							///< Transformation matrix.
		a_vec x;							///< Center point.
		a_mat xyz;						///< Matrix of points.
		a_mat phi;
		double avg_dist; // Average distance betweenn points in a patch
		real_t radio; // radio of a patch
		size_t min_nv;
	
	public:
		static size_t expected_nv;		///< Expected number of patch vertices.

	public:
		patch() = default;
		~patch() = default;
		
		void init(	che * mesh,						///< input mesh.
					const index_t & v,				///< center vertex of the patch.
					const size_t & n_toplevels,		///< number of toplevels to jet fitting.
					const real_t & radio_,		///< euclidean radio in XY of the patch.
					index_t * _toplevel = nullptr		///< aux memory to gather toplevel vertices.
					);
		void init_disjoint(che * mesh,
					const index_t & v,
					const size_t & n_toplevels,
					vector<index_t> & _vertices, 
					index_t * _toplevel = nullptr);
		void init_radial_disjoint(vector<index_t> & idxs_he,
					che * mesh,
					const real_t & radio_,
					const index_t & v,
					real_t & euc_radio,
					real_t & geo_radio,
					double delta,		
					double sum_thres,
					double area_thres);
		void init_random(vertex c, arma::mat T, real_t radio, real_t max_radio, real_t per, real_t fr);
		void recover_radial_disjoint(che * mesh,
					const real_t & radio_,
					const index_t & v);

		void transform();
		
		void itransform();
		
		void reset_xyz(	che * mesh,
						std::vector<vpatches_t> & vpatches,
						const index_t & p,
						const fmask_t & mask = nullptr
						);
		void reset_xyz_disjoint(che * mesh,
						real_t * dist,
						size_t M,
						std::vector<vpatches_t> & vpatches,
						const index_t & p,
						const fmask_t & mask = nullptr
						);
		void remove_extra_xyz_disjoint(size_t & max_points);
		void add_extra_xyz_disjoint(che * mesh, vector<vpatches_t> & vpatches, const index_t & p);
		const a_vec normal();
		bool is_covered( bool * covered);

		void save(const real_t & radio, const size_t & imsize, CImgList<real_t> & imlist);
		void update_heights(real_t & min, real_t & max, bool flag);
		void compute_avg_distance(che * mesh, vector<vpatches_t> & vpatches, const index_t & p);
		void scale_xyz(const real_t & radio_f);
		void iscale_xyz(const real_t & radio_f);
		bool add_vertex_by_faces(const vertex & c,
								vertex & n,
								vector<vertex> & N,
								index_t * indexes,
								size_t nc,
								double thr_angle,
								const geodesics & geo,
								che * mesh,
								const index_t & v,
								double & area,
								double & proj_area,
								double deviation);
		

	private:
		/// Gather the vertices needed to compute the jet_fit_directions of the patch.
		void gather_vertices(	che * mesh,
								const index_t & v,
								const size_t & n_toplevels,
								index_t * toplevel
								);
		
		/// Gather the vertices filter by radio in the local coordinates require initialize T and x.
		void gather_vertices(	che * mesh,
								const index_t & v,
								const real_t & radio,
								index_t * toplevel
								);
		bool exists(index_t idx);
		
		/// Initialize transformation matrix T and translation vector x, using CGAL jet_fitting.
		void jet_fit_directions(che * mesh,
								const index_t & v
								);
		void normal_fit_directions(che * mesh,
								const index_t & v
								);
		real_t get_min_z();
		real_t get_max_z();

		void save_z(ostream & os);
		index_t find(index_t * indexes, size_t nc, index_t idx_global);
		

	friend class dictionary;
};


} // namespace gproshan::mdict

#endif // PATCH_H

