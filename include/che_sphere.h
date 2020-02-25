#ifndef CHE_SPHERE_H
#define CHE_SPHERE_H

#include "che.h"


// geometry processing and shape analysis framework
namespace gproshan {


class che_sphere : public che
{
	private:
		real_t radio;

	public:
		che_sphere(const real_t & r = 1, const size_t & n = 10);
		che_sphere(const che_sphere & mesh);
		virtual ~che_sphere();
};


} // namespace gproshan

#endif // CHE_SPHERE_H

