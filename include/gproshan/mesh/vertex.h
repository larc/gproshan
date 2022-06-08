#ifndef VERTEX_H
#define VERTEX_H

#include <gproshan/include.h>
#include <gproshan/geometry/vec.h>

#include <iostream>


#define glm_vec3(v) glm::vec3((v)[0], (v)[1], (v)[2])


// geometry processing and shape analysis framework
namespace gproshan {


/*!
	The vertex class represents a 3D point and implements 3D vector operations.
*/
using vertex = vec<real_t, 3>;


} // namespace gproshan

#endif // VERTEX_H

