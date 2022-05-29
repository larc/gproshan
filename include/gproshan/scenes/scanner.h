#ifndef SCANNER_H
#define SCANNER_H

#include <gproshan/raytracing/raytracing.h>

#include "mesh/che.h"


// geometry processing and shape analysis framework
// raytracing approach
namespace gproshan::rt {


che * scanner_ptx(const raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam);
che * scanner_ptx(const che * mesh, raytracing * rt, const size_t & n_rows, const size_t & n_cols, const vertex & cam, const std::string & file_jpg = "");


} // namespace gproshan

#endif // SCANNER_H

