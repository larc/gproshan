#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include <gproshan/mesh/che.h>

#include <armadillo>


// geometry processing and shape analysis framework
namespace gproshan {


class descriptor
{
	public:
		enum signature { GPS, HKS, WKS };

	private:
		arma::sp_mat L, A;
		arma::vec eigval;
		arma::mat eigvec;
		arma::mat features;

	public:
		descriptor(const signature & sig, const che * mesh, const size_t n_eigs);
		size_t n_eigs();

		///< return true if the features were computed
		operator bool () const;

		///< return norm of the descriptor for the vertex v
		float operator () (const index_t v) const;

	private:
		void compute_gps();
		void compute_hks();
		void compute_wks();
};


} // namespace gproshan

#endif // DESCRIPTOR_H

