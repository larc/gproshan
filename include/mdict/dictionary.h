#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "che.h"
#include "d_basis.h"
#include "d_mesh.h"

#include <armadillo>

using namespace arma;

// mesh dictionary learning and sparse coding namespace
namespace mdict {

class dictionary
{
	protected:
		che * mesh;
		size_t n_vertices;
		
		basis * phi_basis;

		distance_t f;	// overlapping factor
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

		static const size_t L;

	protected:
		dictionary(che *const & _mesh, basis *const &_phi_basis, const size_t & _m, const size_t & _M, const distance_t & _f, const bool & _plot);
		virtual ~dictionary();

		virtual void execute() = 0;
		
		void learning();
		void sparse_coding();
		void init_sampling();
		void init_patches(const bool & reset = 1, const size_t & threshold = NIL);
		void mesh_reconstruction();

		index_t sample(const index_t & s);
};

} // mdict

#endif // DICTIONARY_H

