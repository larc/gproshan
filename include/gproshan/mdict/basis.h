#ifndef BASIS_H
#define BASIS_H

#include <gproshan/include.h>

#include <fstream>

#include <armadillo>


// geometry processing and shape analysis framework
// mesh dictionary learning and sparse coding namespace
namespace gproshan::mdict {


class basis
{
	protected:
		float _radio;
		size_t _dim;

	public:
		basis(const float r, const size_t d);
		virtual ~basis() = default;
		virtual void discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y) = 0;
		virtual void d_discrete(arma::fmat & phi, const arma::fvec & x, const arma::fvec & y, const bool b) = 0;
		virtual float freq(const index_t idx) = 0;
		float & radio();
		size_t dim() const;
		void plot_basis();
		void plot_atoms(const arma::fmat & A);
		void plot_patch(const arma::fmat & A, const arma::fmat & xyz, const index_t p);

	private:
		virtual void plot_basis(std::ostream & os) = 0;
		virtual void plot_atoms(std::ostream & os, const arma::fvec & A) = 0;
};


} // namespace gproshan::mdict

#endif // BASIS_H

