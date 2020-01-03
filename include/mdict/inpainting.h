#ifndef INPAINTING_H
#define INPAINTING_H

#include "dictionary.h"
#include "../che_poisson.h"
#include "../che_fill_hole.h"
#include "sampling.h"
#include "geodesics.h"
#include "geodesics_ptp.h"
#include <random>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class inpainting : public dictionary
{
	public:
		inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _learn, size_t _avg_p = 36, size_t _perc = 50, const bool & _plot = false);
		virtual ~inpainting() = default;

		distance_t execute();
		void load_mask(const std::vector<index_t> * vertices, const index_t * clusters);
		void init_patches_disjoint();
		void init_voronoi_sampling(const std::vector<index_t> * vertices, const index_t * clusters);
		distance_t execute_tmp();
	private:
		size_t avg_p;
		size_t percent;
		bool * mask;
	
};


} // namespace gproshan::mdict

#endif // INPAINTING_H
