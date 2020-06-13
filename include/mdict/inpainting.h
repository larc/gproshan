#ifndef INPAINTING_H
#define INPAINTING_H

#include "dictionary.h"
#include "../che_poisson.h"
#include "../che_fill_hole.h"
#include "sampling.h"
#include "geodesics.h"
#include "geodesics_ptp.h"
#include <random>
#define PI 3.14159265

// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class inpainting : public dictionary
{
	public:
		inpainting(che *const & _mesh, basis *const & _phi_basis, const size_t & _m, const size_t & _M,
		 const real_t & _f, const bool & _learn, size_t _avg_p = 36, size_t _perc = 0, double _delta=PI/6, double _sum_thres = 0.001 , double _area_thres = 0.001, const bool & _plot = false);
		virtual ~inpainting() = default;

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
		size_t avg_p;
		size_t percent;
		double delta;
		double sum_thres;
		bool * mask;
		double area_thres;
	
};


} // namespace gproshan::mdict

#endif // INPAINTING_H
