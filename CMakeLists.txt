cmake_minimum_required(VERSION 3.22)

project(gproshan VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${gproshan_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wno-unused-result)


find_package(CUDAToolkit 11.7)
if(CUDAToolkit_FOUND)
	enable_language(CUDA)
	include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})

	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	set(CMAKE_CUDA_ARCHITECTURES 52 60 61 70 75 80 86)

	find_package(OptiX 7.4)
	if(OptiX_INCLUDE)
		include_directories(SYSTEM ${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDAToolkit_FOUND)


set(GPROSHAN_DIR "${gproshan_SOURCE_DIR}")

option(GPROSHAN_LOG "Enables logging gproshan messages." ON)
option(GPROSHAN_FLOAT "Compile gproshan with float precision." ON)

set(GPROSHAN_CUDA ${CUDAToolkit_FOUND})
set(GPROSHAN_OPTIX ${OptiX_INCLUDE})


configure_file(	${gproshan_SOURCE_DIR}/include/gproshan/config.h.in
				${gproshan_SOURCE_DIR}/include/gproshan/config.h
				)


find_package(embree 3.13 REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(Boost COMPONENTS thread system)


include_directories(SYSTEM ${embree_INCLUDE_DIRS})
include_directories(SYSTEM ${GLEW_INCLUDE_DIRS})
include_directories(SYSTEM ${GLFW3_INCLUDE_DIRS})
include_directories(SYSTEM ${AMADILLO_INCLUDE_DIR})
include_directories(SYSTEM ${SuiteSparse_INCLUDE_DIRS})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})


add_subdirectory(src)			# gproshan library
add_subdirectory(apps)			# gproshan executables

