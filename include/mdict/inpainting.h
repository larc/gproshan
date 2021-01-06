#ifndef INPAINTING_H
#define INPAINTING_H

#include "mdict/msparse_coding.h"
#include "mesh/che_poisson.h"
#include "mesh/che_fill_hole.h"

#include <random>

// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class inpainting : public msparse_coding
{
	public:
		inpainting(che *const & _mesh, basis *const & _phi_basis, const params & p);
		virtual ~inpainting();

		operator const std::string & () const;

		real_t execute();
		void load_mask(const std::vector<index_t> * vertices, const index_t * clusters);
		void load_mask();
		void init_voronoi_patches();
		void init_radial_feature_patches();
		void load_sampling(bool save_all);
		che * point_cloud_reconstruction(real_t per, real_t fr);
		vector<index_t> sort_indexes(const vector<real_t> &v);

		real_t execute_tmp();

	private:
		bool * mask;
		std::string key_name;
};


} // namespace gproshan::mdict

#endif // INPAINTING_H

