cmake_minimum_required(VERSION 3.20)

project(gproshan VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${gproshan_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra -Wno-unused-result)
add_definitions(-DGPROSHAN_FLOAT)
add_definitions(-DGPROSHAN_LOG)


find_package(CUDAToolkit 11.4)
if(CUDAToolkit_FOUND)
	enable_language(CUDA)
	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	set(CMAKE_CUDA_ARCHITECTURES 52 60 61 70 75 80 86)

	add_definitions(-DGPROSHAN_CUDA)
	include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})

	find_package(OptiX 7.3)
	if(OptiX_INCLUDE)
		add_definitions(-DGPROSHAN_OPTIX)
		include_directories(SYSTEM ${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDAToolkit_FOUND)


find_package(embree 3.13 REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(Boost COMPONENTS thread system)


include_directories(SYSTEM ${embree_INCLUDE_DIRS})
include_directories(SYSTEM ${GLEW_INCLUDE_DIRS})
include_directories(SYSTEM ${GLFW3_INCLUDE_DIRS})
include_directories(SYSTEM ${AMADILLO_INCLUDE_DIR})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(SYSTEM ${SuiteSparse_INCLUDE_DIRS})
include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})


add_subdirectory(imgui)
include_directories(SYSTEM ${gproshan_SOURCE_DIR}/imgui)


include_directories(${gproshan_SOURCE_DIR}/include/gproshan)
add_compile_definitions(GPROSHAN_DIR="${gproshan_SOURCE_DIR}")


FILE(GLOB_RECURSE cpp_sources ${gproshan_SOURCE_DIR}/src/*.cpp)

add_library(gproshan SHARED ${cpp_sources})

target_link_libraries(gproshan embree)
target_link_libraries(gproshan OpenMP::OpenMP_CXX)
target_link_libraries(gproshan GLEW::GLEW)
target_link_libraries(gproshan glfw)
target_link_libraries(gproshan ${X11_X11_LIB})
target_link_libraries(gproshan ${ARMADILLO_LIBRARIES})
target_link_libraries(gproshan ${SuiteSparse_LIBRARIES})
target_link_libraries(gproshan CGAL::CGAL)
target_link_libraries(gproshan ${OptiX_LIBRARY})
target_link_libraries(gproshan imgui)

if(CUDAToolkit_FOUND)
	FILE(GLOB_RECURSE cu_sources ${gproshan_SOURCE_DIR}/src/*.cu)
	list(REMOVE_ITEM cu_sources ${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu)

	add_library(gproshan_cuda SHARED ${cu_sources})

	target_link_libraries(gproshan_cuda CUDA::cudart)
	target_link_libraries(gproshan_cuda CUDA::cuda_driver)
	target_link_libraries(gproshan_cuda CUDA::cublas)
	target_link_libraries(gproshan_cuda CUDA::cusolver)
	target_link_libraries(gproshan_cuda CUDA::cusparse)

	target_link_libraries(gproshan gproshan_cuda)
endif(CUDAToolkit_FOUND)


add_subdirectory(apps)

