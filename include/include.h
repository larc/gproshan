#ifndef INCLUDE_H
#define INCLUDE_H

#include "config.h"

#include <map>
#include <iomanip>
#include <string>

#define NIL (0u - 1)

typedef unsigned int index_t;
typedef unsigned int length_t;
typedef double matrix_t;
typedef double percent_t;

#ifdef SINGLE_P
	typedef float vertex_t;
#else
	typedef double vertex_t;
#endif

typedef vertex_t distance_t;
typedef vertex_t area_t;
typedef vertex_t angle_t;
typedef std::pair<area_t,size_t> areapair_t;

#define PATH_MDATA	std::string("../DATA/")
#define PATH_DATA	std::string("../DATA/shapes_grande/")
#define PATH_CORR	std::string("../DATA/corr_shapes/")
#define PATH_KPS	std::string("../TEST/keypoints/")
#define PATH_KCS	std::string("../TEST/keycomponents/")
#define PATH_HOLES	std::string("../TEST/holes/")
#define PATH_TEST	std::string("../TEST/")


#ifndef NDEBUG
	#define debug(vari) cerr << "\033[0;33m" << setprecision(3) << scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << endl;
	#define debug_me(message) fprintf(stderr, "\033[1;31m%s: %s (%s:%d)\n\033[0m", #message, __FUNCTION__, __FILE__, __LINE__);
	#define d_message(message) fprintf(stderr, "\033[0;33m%s\n\033[0m", #message);
#else
	#define debug(vari)
	#define debug_me(message)
	#define d_message(message)
#endif


#define TIC(t) (t) = omp_get_wtime();
#define TOC(t) (t) = omp_get_wtime() - (t);

#endif // INCLUDE_H

