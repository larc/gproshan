#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "che.h"
#include "patch.h"
#include "d_dict_learning.h"
#include "d_basis.h"
#include "d_mesh.h"

#include "include_arma.h"

/// Mesh dictionary learning and sparse coding namespace
namespace mdict {

class dictionary
{
	protected:
		che * mesh;								///< input mesh.
		size_t n_vertices;						///< number of vertices.

		basis * phi_basis;						///< continuous basis.

		size_t m;								///< number of dictionary atoms.
		size_t M;								///< number of patches.
		a_mat A;									///< dictionary continuous matrix.
		a_mat alpha;								///< sparse coding matrix.
		
		distance_t f;
		distance_t s_radio;						///< sampling geodesic radio.
		vector<index_t> sampling;				///< samples, center of patches if sampling.
		vector<patch> patches;				///< vector of patches.
		vector<vpatches_t> patches_map;		///< invert index vertex to patches.

		double d_time;							///< time of operations.
		bool d_plot;
									///< plot atoms and basis with gnuplot.
		real_t minz;
		real_t maxz;
	
	public:
		static size_t L;					///< sparsity, norm L_0, default 10.
		static size_t T;					///< factor of patches' size, default 5 toplesets.

	protected:
		dictionary(	che *const & _mesh, 		///< pointer to input mesh.
					basis *const &_phi_basis,	///< pointer to continuous basis.
					const size_t & _m,			///< number of dictionary atoms.
					const size_t & _M,			///< number of patches.
					const distance_t & _f,		///< deprecated
					const bool & _plot			///< flag to plot basis and atoms with gnuplot.
					);

		virtual ~dictionary();

		virtual void execute() = 0;

		void learning();
		void sparse_coding();
		void init_sampling();
		void init_patches(	const bool & reset = 1,
							const fmask_t & mask = nullptr
							);

		void mesh_reconstruction();
		void update_alphas(a_mat & alpha, size_t threshold);

		index_t sample(const index_t & s);
};

} // mdict

#endif // DICTIONARY_H

