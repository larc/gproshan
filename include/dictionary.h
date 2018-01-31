#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "che.h"
#include "d_basis.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

class dictionary
{
	protected:
		che * mesh;
		size_t n_vertices;
		
		basis * phi_basis;
		size_t m;		// number of atoms
		size_t M; 		// number of patches
		mat A;			// dictionary continuous matrix
		mat alpha;		// sparse coding matrix

		distance_t s_radio;
		vector<index_t> sampling;
		vector<patch_t> patches;
		vector<patches_map_t> patches_map;
		
		double d_time;
		bool d_plot;

		static const size_t min_nvp;
		static const size_t L;

	protected:
		dictionary(che *const & _mesh, basis *const &_phi_basis, const size_t & _m, const size_t & _M, const distance_t & f, const bool & _plot = false);
		virtual ~dictionary();

	public:
		void learning();
		virtual void execute() = 0;
		
	protected:
		void init_patches();
		index_t sample(const index_t & s);
};

#endif // DICTIONARY_H

