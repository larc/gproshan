#ifndef KEY_COMPONENTS_H
#define KEY_COMPONENTS_H

#include <gproshan/mesh/che.h>

#include <map>


// geometry processing and shape analysis framework
namespace gproshan {


class key_components
{
	private:
		real_t radio		= 0;
		size_t n_vertices	= 0;
		size_t n_comp		= 0;
		index_t * comp		= nullptr;
		size_t * comp_size	= nullptr;
		std::map<index_t, index_t> comp_idx;

	public:
		key_components(che * mesh, const std::vector<index_t> & kps, const real_t r);
		~key_components();
		index_t operator()(const index_t i);
		operator size_t () const;

	private:
		void compute_kcs(che * mesh, const std::vector<index_t> & kps);
		index_t find(const index_t x);
		bool join(index_t x, index_t y);
};


} // namespace gproshan

#endif // KEY_COMPONENTS_H

