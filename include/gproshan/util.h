#ifndef UTIL_H
#define UTIL_H

#include <gproshan/include.h>


// geometry processing and shape analysis framework
namespace gproshan {


void copy_real_t_array(float * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(float * destination, const double * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const float * source, const size_t & n_elem);

void copy_real_t_array(double * destination, const double * source, const size_t & n_elem);


} // namespace gproshan

#endif // UTIL_H

