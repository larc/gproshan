#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "che.h"
#include "patch.h"
#include "mdict.h"
#include "basis.h"
#include "d_mesh.h"

#include "include_arma.h"


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class dictionary
{
	public:
		struct params
		{
			size_t m = 144;				///< number of dictionary atoms
			size_t M = 0;				///< number of patches
			size_t avg_p = 36;			///< avg number of vertices per patch
			size_t percentage = 0;		///< mask percentage
			real_t f = 1;				///<
			real_t delta = M_PI / 6;	///<
			real_t sum_thres = 1.01;	///<
			real_t area_thres = 0.005;	///<
			bool learn = false;			///<
			bool plot = false;
		};

	protected:
		che * mesh;								///< input mesh.
		size_t n_vertices;						///< number of vertices.

		basis * phi_basis;						///< continuous basis.

		size_t m;								///< number of dictionary atoms.
		size_t M;								///< number of patches.
		a_mat A;								///< dictionary continuous matrix.
		a_mat alpha;							///< sparse coding matrix.
		
		real_t f;
		real_t s_radio;							///< sampling geodesic radio.
		std::vector<index_t> sampling;			///< samples, center of patches if sampling.
		std::vector<patch> patches;				///< vector of patches.
		std::vector<vpatches_t> patches_map;	///< invert index vertex to patches.

		double d_time;							///< time of operations.
		bool learn;
		bool d_plot;							///< plot atoms and basis with gnuplot.
		real_t * dist;
	
	public:
		static size_t K;						///< number of iterations KSVD.
		static size_t L;						///< sparsity, norm L_0, default 10.
		static size_t T;						///< factor of patches' size, default 5 toplesets.
	
	public:
		const real_t & operator[](const index_t & i) const;
		void draw_patches(index_t i);

	protected:
		dictionary(	che *const & _mesh, 		///< pointer to input mesh.
					basis *const &_phi_basis,	///< pointer to continuous basis.
					const size_t & _m,			///< number of dictionary atoms.
					const size_t & _M,			///< number of patches.
					const real_t & _f,		///< deprecated
					const bool & _learn,		
					const bool & _plot			///< flag to plot basis and atoms with gnuplot.
					);

		virtual ~dictionary();

		virtual real_t execute() = 0;

		void learning();
		void sparse_coding();
		void init_sampling();
		void load_curvatures(a_vec & curvatures);
		void load_features(vector<index_t> & v_feat, size_t & featsize);
		void init_patches(	const bool & reset = 1,
							const fmask_t & mask = nullptr
							);

		real_t mesh_reconstruction(const fmask_t & mask = nullptr);
		void update_alphas(a_mat & alpha, size_t threshold);
		void save_alpha(string file);

		index_t sample(const index_t & s);
		
};


} // namespace gproshan::mdict

#endif // DICTIONARY_H

