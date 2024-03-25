#ifndef INCLUDE_H
#define INCLUDE_H

#include <gproshan/config.h>

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <string>

#include <omp.h>

#define NIL (0u - 1)

#ifdef __CUDACC__
	#define __host_device__ __host__ __device__
#else
	#define __host_device__
#endif // __CUDACC__


// geometry processing and shape analysis framework
namespace gproshan {


using index_t = unsigned int;


inline std::string tmp_file_path(const std::string & file)
{
	const std::string & gproshan_home = std::string(getenv("HOME")) + "/.gproshan";
	std::filesystem::create_directory(gproshan_home);

	return gproshan_home + "/" + file;
}


#define shaders_path(file) (std::string(GPROSHAN_DIR) + "/shaders/" + file)


#ifdef GPROSHAN_LOG
	#define gproshan_log_var(vari) std::cerr << "\033[0;33m[LOG] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << std::endl
	#define gproshan_log(message) fprintf(stderr, "\033[1;31m[LOG] %s: %s\n\033[0m", #message, __FUNCTION__)
#else
	#define gproshan_log_var(vari)
	#define gproshan_log(message)
#endif

#define gproshan_input(input) fprintf(stderr, "\033[38;5;210m[INPUT] %s: \033[0m", #input)

#define gproshan_error_var(vari) std::cerr << "\033[1;31m[ERROR] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << std::endl
#define gproshan_error(message) fprintf(stderr, "\033[1;31m[ERROR]\033[0m %s\033[38;5;210m\t%s:%d > %s\n\033[0m", #message, __FILE__, __LINE__, __FUNCTION__)

#ifndef NDEBUG
	#define gproshan_debug_var(vari) std::cerr << "\033[1;31m[DEBUG] " << std::setprecision(3) << std::scientific << #vari << ":\033[0m " << (vari) << "\t\033[38;5;210m" << __FILE__ ":" << __LINE__ << " > " << __FUNCTION__ << "\033[0m" << std::endl
	#define gproshan_debug(message) fprintf(stderr, "\033[1;31m[DEBUG]\033[0m %s\033[38;5;210m\t%s:%d > %s\n\033[0m", #message, __FILE__, __LINE__, __FUNCTION__)
#else
	#define gproshan_debug_var(vari)
	#define gproshan_debug(message)
#endif


#define TIC(t) (t) = omp_get_wtime();
#define TOC(t) (t) = omp_get_wtime() - (t);


} // namespace gproshan

#endif // INCLUDE_H

