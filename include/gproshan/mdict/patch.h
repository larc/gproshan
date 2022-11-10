#ifndef PATCH_H
#define PATCH_H

#include <gproshan/include.h>
#include <gproshan/mesh/che.h>

#include <vector>
#include <map>
#include <algorithm>


#include <gproshan/include_arma.h>


#ifdef Success
	#undef Success
#endif


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class msparse_coding;

typedef function<bool(const index_t &)> fmask_t;
typedef std::map<index_t, index_t> vpatches_t;

///
class patch
{
	public:
		std::vector<index_t> vertices;		///< Vertices of the patch.
		a_mat T;							///< Transformation matrix.
		a_vec x;							///< Center point.
		a_mat xyz;							///< Matrix of points.
		a_mat phi;							///< Projected basis.
		double avg_dist;					///< Average distance between points.
		real_t radio;						///< Radio.
		size_t min_nv;						///<

	public:
		static size_t expected_nv;			///< Expected number of patch vertices.
		static real_t nyquist_factor;		///< nyquist factor

	public:
		patch() = default;
		~patch() = default;

		void init(	che * mesh,						///< input mesh.
					const index_t & v,				///< center vertex of the patch.
					const size_t & n_toplevels,		///< number of toplevels to jet fitting.
					const real_t & radio_,			///< euclidean radio in XY of the patch.
					index_t * _toplevel = nullptr	///< aux memory to gather toplevel vertices.
					);

		void init_disjoint(che * mesh,
					const index_t & v,
					const size_t & n_toplevels,
					vector<index_t> & _vertices,
					index_t * _toplevel = nullptr);

		void init_radial_disjoint(	real_t & euc_radio,
									real_t & geo_radio,
									che * mesh,
									const index_t & v,
									const real_t & delta,
									const real_t & sum_thres,
									const real_t & area_thres,
									const real_t & area_mesh
									);

		void init_random(const vertex & c, const a_mat & T, const real_t & radio, const real_t & max_radio, const real_t & percent, const real_t & fr);

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

//		void save(const real_t & radio, const size_t & imsize, CImgList<real_t> & imlist);
		void update_heights(real_t & min, real_t & max, bool flag);
		void compute_avg_distance(che * mesh, vector<vpatches_t> & vpatches, const index_t & p);
		void scale_xyz(const real_t & radio_f);
		void iscale_xyz(const real_t & radio_f);
		bool add_vertex_by_faces(	vertex & n,
									vector<vertex> & N,
									double thr_angle,
									const real_t * geo,
									che * mesh,
									const index_t & v,
									real_t & area,
									real_t & proj_area,
									real_t deviation
									);


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
		index_t find(const index_t * indexes, size_t nc, index_t idx_global);


	friend class msparse_coding;
};


} // namespace gproshan::mdict

#endif // PATCH_H

