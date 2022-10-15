#include <gproshan/geometry/mat.h>

#include <gproshan/include_arma.h>


// geometry processing and shape analysis framework
namespace gproshan {


mat4 inverse(const mat4 & m)
{
	mat4 inv;
	mat4 mt = mat4::transpose(m);

	a_mat a_inv((real_t *) &inv, 4, 4, false, true);
	a_mat a_m((real_t *) &mt, 4, 4, false, true);

	arma::inv(a_inv, a_m);

	return mat4::transpose(inv);
}


} // namespace gproshan

