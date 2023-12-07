#ifndef SCANNER_H
#define SCANNER_H

#include <gproshan/mesh/che.h>
#include <gproshan/raytracing/raytracing.h>


// geometry processing and shape analysis framework
namespace gproshan {


che * scanner_ptx(const rt::raytracing * rt, const size_t n_rows, const size_t n_cols, const vertex & cam_pos);

che * scanner_ptx_jpg(const rt::raytracing * rt, const size_t n_rows, const size_t n_cols, const vertex & cam_pos, const std::string & file_jpg = "");


} // namespace gproshan

#endif // SCANNER_H

