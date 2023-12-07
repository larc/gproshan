#ifndef CHE_SPHERE_H
#define CHE_SPHERE_H

#include <gproshan/mesh/che.h>


// geometry processing and shape analysis framework
namespace gproshan {


class che_sphere : public che
{
	public:
		che_sphere(const real_t r = 1, const size_t & n = 6);
};


} // namespace gproshan

#endif // CHE_SPHERE_H

