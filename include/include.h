#ifndef INCLUDE_H
#define INCLUDE_H

#include "config.h"

#include <map>
#include <iomanip>
#include <string>

#include <omp.h>

#define NIL (0u - 1)


// geometry processing and shape analysis framework
namespace gproshan {


typedef unsigned int index_t;

#ifdef SINGLE_P
	typedef float real_t;
#else
	typedef double real_t;
#endif

typedef real_t distance_t;
typedef real_t area_t;
typedef real_t angle_t;


#define TMP_DIR	"../tmp/"
#define tmp_file_path(file) (std::string(TMP_DIR) + file)


#ifndef NDEBUG
	#define debug(vari) std::cerr << "\033[0;33m" << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << std::endl;
	#define debug_me(message) fprintf(stderr, "\033[1;31m%s: %s (%s:%d)\n\033[0m", #message, __FUNCTION__, __FILE__, __LINE__);
	#define d_message(message) fprintf(stderr, "\033[0;33m%s\n\033[0m", #message);
#else
	#define debug(vari) ;
	#define debug_me(message) ;
	#define d_message(message) ;
#endif


#define TIC(t) (t) = omp_get_wtime();
#define TOC(t) (t) = omp_get_wtime() - (t);


} // namespace gproshan

#endif // INCLUDE_H

