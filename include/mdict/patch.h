#ifndef PATCH_H
#define PATCH_H

#include "include.h"
#include "che.h"

#include <vector>
#include <armadillo>

#ifndef CGAL_PATCH_DEFS
	#define CGAL_PATCH_DEFS
	#define CGAL_EIGEN3_ENABLED
	#define CGAL_USE_BOOST_PROGRAM_OPTIONS
	#define CGAL_USE_GMP
	#define DCGAL_USE_MPFR
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>

using namespace std;
using namespace arma;

/// Mesh dictionary learning and sparse coding namespace
namespace mdict {

class dictionary;

/// 
class patch
{
	private:
		vector<index_t> vertices;		///< Vertices of the patch.
		mat T;							///< Transformation matrix.
		vec x;							///< Center point.
		mat xyz;						///< Matrix of points.
	
	public:
		static size_t expected_nv;		///< Expected number of patch vertices.

	public:
		patch() = default;
		~patch() = default;
		void init(che * mesh, const index_t & v, index_t * _level = NULL);
		operator const vector<index_t> & () const;

	private:
		void clear();
	friend class dictionary;
};

} // mdict

#endif // PATCH_H

