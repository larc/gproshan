#ifndef KEY_COMPONENTS_H
#define KEY_COMPONENTS_H

#include "mesh/che.h"
#include "features/key_points.h"

#include <map>


// geometry processing and shape analysis framework
namespace gproshan {


class key_components
{
	private:
		real_t radio;
		size_t n_vertices;
		size_t n_comp;
		index_t * comp;
		size_t * comp_size;
		std::map<index_t, index_t> comp_idx;

	public:
		key_components(che * mesh, const key_points & kps, const real_t & r);
		~key_components();
		index_t operator()(const index_t & i);
		operator const size_t & () const;
	
	private:
		void compute_kcs(che * mesh, const key_points & kps);
		index_t find(const index_t & x);
		bool join(index_t x, index_t y);
};


} // namespace gproshan

#endif // KEY_COMPONENTS_H

