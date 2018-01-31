#ifndef DENOISING_H
#define DENOISING_H

#include "dictionary.h"

#include <armadillo>

using namespace arma;

class denoising : public dictionary
{
	public:
		denoising(che *const & _mesh, basis *const &_phi_basis, const size_t & _m, const size_t & _M,
         const distance_t & f, const bool & _plot = false): dictionary(_mesh, phi_basis, _m, _M, f, _plot){}
		virtual ~denoising();
        void execute();

};

#endif // DICTIONARY_H
