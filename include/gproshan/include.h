#ifndef INCLUDE_H
#define INCLUDE_H

#include <iostream>
#include <iomanip>
#include <string>

#include <omp.h>

#define NIL (0u - 1)


// geometry processing and shape analysis framework
namespace gproshan {


typedef unsigned int index_t;

#ifdef GPROSHAN_FLOAT
	typedef float real_t;
#else
	typedef double real_t;
#endif


#define tmp_file_path(file) (std::string(GPROSHAN_DIR) + "/tmp/" + file)
#define shaders_path(file) (std::string(GPROSHAN_DIR) + "/shaders/" + file)

#ifdef GPROSHAN_LOG
	#define gproshan_log_var(vari) std::cerr << "\033[0;33m[LOG] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << std::endl
	#define gproshan_log(message) fprintf(stderr, "\033[1;31m[LOG] %s: %s\n\033[0m", #message, __FUNCTION__)
#else
	#define gproshan_log_var(vari)
	#define gproshan_log(message)
#endif

#define gproshan_input(input) fprintf(stderr, "\033[38;5;210m[INPUT] %s: \033[0m", #input)

#define gproshan_error_var(vari) std::cerr << "\033[0;33m[ERROR] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << std::endl
#define gproshan_error(message) fprintf(stderr, "\033[1;31m[ERROR] %s: %s (%s:%d)\n\033[0m", #message, __FUNCTION__, __FILE__, __LINE__)

#ifndef NDEBUG
	#define gproshan_debug_var(vari) std::cerr << "\033[0;33m[DEBUG] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << std::endl
	#define gproshan_debug(message) fprintf(stderr, "\033[1;31m[DEBUG] %s: %s (%s:%d)\n\033[0m", #message, __FUNCTION__, __FILE__, __LINE__)
#else
	#define gproshan_debug_var(vari)
	#define gproshan_debug(message)
#endif


#define TIC(t) (t) = omp_get_wtime();
#define TOC(t) (t) = omp_get_wtime() - (t);


} // namespace gproshan

#endif // INCLUDE_H

