#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "che.h"
#include "d_basis.h"

class dictionary
{
	private:
		che * mesh;
		basis * phi_basis;
		size_t m;		// number of atoms
		size_t M; 		// number of patches

	public:
		dictionary(che *const & _mesh, basis *const &_phi_basis, const size_t & _m, const size_t & _M);
		virtual ~dictionary();
};

#endif // DICTIONARY_H

