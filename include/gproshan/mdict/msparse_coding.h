#ifndef MSPARSE_CODING_H
#define MSPARSE_CODING_H

#include <gproshan/mesh/che.h>
#include <gproshan/mdict/patch.h>
#include <gproshan/mdict/mdict.h>
#include <gproshan/mdict/basis.h>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class msparse_coding
{
	public:
		struct params
		{
			size_t n_atoms		= 144;			///< number of dictionary atoms
			size_t n_patches	= 0;			///< number of patches
			size_t avg_p		= 36;			///< avg number of vertices per patch
			size_t percent		= 0;			///< mask percentage
			float f			= 1;			///<
			float delta		= M_PI / 6;		///<
			float sum_thres	= 1.01;			///<
			float area_thres	= 0.005;		///<
			bool learn			= false;		///<
			bool plot			= false;		///<
		};

	private:
		che * mesh;								///< input mesh.
		size_t n_vertices;						///< number of vertices.

		basis * phi_basis;						///< continuous basis.
		params m_params;						///<

		arma::fmat A;								///< dictionary continuous matrix.
		arma::fmat alpha;							///< sparse coding matrix.

		float s_radio;							///< sampling geodesic radio.
		std::vector<index_t> sampling;			///< samples, center of patches if sampling.
		std::vector<patch> patches;				///< vector of patches.
		std::vector<vpatches_t> patches_map;	///< invert index vertex to patches.
		std::vector<std::pair<float, index_t> > patches_error;

		double d_time;							///< time of operations.
		float * dist;

		bool * mask = nullptr;
		std::string key_name;

	public:
		static size_t K;						///< number of iterations KSVD.
		static size_t L;						///< sparsity, norm L_0, default 10.
		static size_t T;						///< factor of patches' size, default 5 toplesets.

	public:
		msparse_coding(	che * _mesh, 			///< pointer to input mesh.
						basis * _phi_basis,		///< pointer to continuous basis.
						const params & p		///<
						);

		virtual ~msparse_coding();

		float operator[](const index_t i) const;
		index_t draw_patches(const index_t p) const;

		operator const std::string & () const;
		operator const float * () const;

		float execute();
		void load_mask(const std::vector<index_t> * vertices, const index_t * clusters);
		void load_mask();
		void init_voronoi_patches();
		void init_radial_feature_patches();
		void load_sampling();
		che * point_cloud_reconstruction(float per, float fr);
		std::vector<index_t> sort_indexes(const std::vector<float> &v);

		float execute_tmp();

	private:
		void learning();
		void sparse_coding();
		void init_sampling();
		void load_features(std::vector<index_t> & v_feat, size_t & featsize);
		void init_patches(	const bool reset = 1,
							const fmask_t & mask = nullptr
							);

		float mesh_reconstruction(const fmask_t & mask = nullptr);
		void update_alphas(arma::fmat & alpha, size_t threshold);
		void save_alpha(std::string file);

		index_t sample(const index_t s);
};


} // namespace gproshan::mdict

#endif // MSPARSE_CODING_H

