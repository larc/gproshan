#ifndef PATCH_H
#define PATCH_H

#include "include.h"
#include "che.h"

#include <vector>
#include "include_arma.h"

#include <CImg.h>

#ifdef Success
  #undef Success
#endif

using namespace cimg_library;

using namespace std;

/// Mesh dictionary learning and sparse coding namespace
namespace mdict {

class dictionary;

typedef function<bool(const index_t &)> fmask_t;
typedef vector<pair<index_t, index_t> > vpatches_t;

/// 
class patch
{
	public:
		vector<index_t> vertices;		///< Vertices of the patch.
		a_mat T;							///< Transformation matrix.
		a_vec x;							///< Center point.
		a_mat xyz;						///< Matrix of points.
		a_mat phi;
	
	public:
		static size_t expected_nv;		///< Expected number of patch vertices.

	public:
		patch() = default;
		~patch() = default;
		
		void init(	che * mesh,						///< input mesh.
					const index_t & v,				///< center vertex of the patch.
					const size_t & n_toplevels,		///< number of toplevels to jet fitting.
					const distance_t & radio,		///< euclidean radio in XY of the patch.
					index_t * _toplevel = NULL		///< aux memory to gather toplevel vertices.
					);

		void transform();
		
		void itransform();
		
		void reset_xyz(	che * mesh,
						vector<vpatches_t> & vpatches,
						const index_t & p,
						const fmask_t & mask = nullptr
						);
		
		const a_vec normal();

		void save(const real_t & radio, const size_t & imsize, CImgList<real_t> & imlist);

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
								const distance_t & radio,
								index_t * toplevel
								);
		
		/// Initialize transformation matrix T and translation vector x, using CGAL jet_fitting.
		void jet_fit_directions(che * mesh,
								const index_t & v
								);
		

	friend class dictionary;
};

} // mdict

#endif // PATCH_H

