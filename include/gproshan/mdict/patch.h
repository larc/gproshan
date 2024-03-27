#ifndef PATCH_H
#define PATCH_H

#include <gproshan/include.h>
#include <gproshan/mesh/che.h>

#include <vector>
#include <map>
#include <algorithm>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class msparse_coding;

typedef std::function<bool(const index_t)> fmask_t;
typedef std::map<index_t, index_t> vpatches_t;

///
class patch
{
	public:
		std::vector<index_t> vertices;		///< Vertices of the patch.
		arma::fmat T;							///< Transformation matrix.
		arma::fvec x;							///< Center point.
		arma::fmat xyz;							///< Matrix of points.
		arma::fmat phi;							///< Projected basis.
		double avg_dist;					///< Average distance between points.
		float radio;						///< Radio.
		size_t min_nv;						///<

	public:
		static size_t expected_nv;			///< Expected number of patch vertices.
		static float nyquist_factor;		///< nyquist factor

	public:
		patch() = default;
		~patch() = default;

		void init(	che * mesh,						///< input mesh.
					const index_t v,				///< center vertex of the patch.
					const size_t n_toplevels,		///< number of toplevels to jet fitting.
					const float radio_,			///< euclidean radio in XY of the patch.
					index_t * _toplevel = nullptr	///< aux memory to gather toplevel vertices.
					);

		void init_disjoint(che * mesh,
					const index_t v,
					const size_t n_toplevels,
					std::vector<index_t> & _vertices,
					index_t * _toplevel = nullptr);

		void init_radial_disjoint(	float & euc_radio,
									float & geo_radio,
									che * mesh,
									const index_t v,
									const float delta,
									const float sum_thres,
									const float area_thres,
									const float area_mesh
									);

		void init_random(const vertex & c, const arma::fmat & T, const float radio, const float max_radio, const float percent, const float fr);

		void recover_radial_disjoint(che * mesh,
					const float radio_,
					const index_t v);

		void transform();

		void itransform();

		void reset_xyz(	che * mesh,
						std::vector<vpatches_t> & vpatches,
						const index_t p,
						const fmask_t & mask = nullptr
						);
		void reset_xyz_disjoint(che * mesh,
						float * dist,
						size_t M,
						std::vector<vpatches_t> & vpatches,
						const index_t p,
						const fmask_t & mask = nullptr
						);
		void remove_extra_xyz_disjoint(size_t & max_points);
		void add_extra_xyz_disjoint(che * mesh, std::vector<vpatches_t> & vpatches, const index_t p);
		const arma::fvec normal();
		bool is_covered( bool * covered);

//		void save(const float radio, const size_t imsize, CImgList<float> & imlist);
		void update_heights(float & min, float & max, bool flag);
		void compute_avg_distance(che * mesh, std::vector<vpatches_t> & vpatches, const index_t p);
		void scale_xyz(const float radio_f);
		void iscale_xyz(const float radio_f);
		bool add_vertex_by_trigs(	vertex & n,
									std::vector<vertex> & N,
									double thr_angle,
									const float * geo,
									che * mesh,
									const index_t v,
									float & area,
									float & proj_area,
									float deviation
									);


	private:
		/// Gather the vertices needed to compute the jet_fit_directions of the patch.
		void gather_vertices(	che * mesh,
								const index_t v,
								const size_t n_toplevels,
								index_t * toplevel
								);

		/// Gather the vertices filter by radio in the local coordinates require initialize T and x.
		void gather_vertices(	che * mesh,
								const index_t v,
								const float radio,
								index_t * toplevel
								);
		bool exists(index_t idx);

		/// Initialize transformation matrix T and translation vector x, using CGAL jet_fitting.
		void normal_fit_directions(che * mesh, const index_t v);
		float get_min_z();
		float get_max_z();

		void save_z(std::ostream & os);
		index_t find(const index_t * indexes, size_t nc, index_t idx_global);


	friend class msparse_coding;
};


} // namespace gproshan::mdict

#endif // PATCH_H

