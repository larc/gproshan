cmake_minimum_required(VERSION 3.28)

project(gproshan VERSION 4.0)

list(APPEND CMAKE_MODULE_PATH "${gproshan_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wno-unused-result)


find_package(CUDAToolkit 12)
if(CUDAToolkit_FOUND)
	enable_language(CUDA)
	include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})

	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	if(NOT DEFINED ENV{GITHUB_ACTIONS})
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()

	find_package(OptiX 8)
	if(OptiX_INCLUDE)
		include_directories(SYSTEM ${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDAToolkit_FOUND)


set(GPROSHAN_DIR "${gproshan_SOURCE_DIR}")

option(GPROSHAN_LOG "Enables logging gproshan messages." ON)

set(GPROSHAN_CUDA ${CUDAToolkit_FOUND})
set(GPROSHAN_OPTIX ${OptiX_INCLUDE})


configure_file(	${gproshan_SOURCE_DIR}/include/gproshan/config.h.in
				${gproshan_SOURCE_DIR}/include/gproshan/config.h
				)


find_package(embree 4 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(flann REQUIRED)


include_directories(SYSTEM ${embree_INCLUDE_DIRS})
include_directories(SYSTEM ${GLEW_INCLUDE_DIRS})
include_directories(SYSTEM ${GLFW3_INCLUDE_DIRS})
include_directories(SYSTEM ${AMADILLO_INCLUDE_DIR})
include_directories(SYSTEM ${SuiteSparse_INCLUDE_DIRS})
include_directories(SYSTEM ${flann_INCLUDE_DIRS})


add_subdirectory(src)			# gproshan library
add_subdirectory(apps)			# gproshan executables

