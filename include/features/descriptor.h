#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include "mesh/che.h"
#include "include_arma.h"


// geometry processing and shape analysis framework
namespace gproshan {


class descriptor
{
	public:
		enum signature { GPS, HKS, WKS };
		
	private:
		a_sp_mat L, A;
		a_vec eigval;
		a_mat eigvec;
	
	public:
		descriptor(const signature & sig, const che * mesh, const size_t & n_eigs);

	private:
};


} // namespace gproshan

#endif // DESCRIPTOR_H

