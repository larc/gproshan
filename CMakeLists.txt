cmake_minimum_required(VERSION 3.18)

project(gproshan VERSION 3.0)

list(APPEND CMAKE_MODULE_PATH "${gproshan_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CUDAToolkit 11)
if(CUDAToolkit_FOUND)
	enable_language(CUDA)
	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	set(CMAKE_CUDA_ARCHITECTURES 52 53 60 61 62 70 72 75 80)

	add_definitions(-DGPROSHAN_CUDA)
	include_directories(${CUDAToolkit_INCLUDE_DIRS})
	
	find_package(OptiX 7.2)
	if(OptiX_INCLUDE)
#		set(CUDA_HOST_COMPILER ${CMAKE_CUDA_HOST_COMPILER})
		add_definitions(-DGPROSHAN_OPTIX)
		include_directories(${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDAToolkit_FOUND)


find_package(embree 3.10)
if(embree_FOUND)
	add_definitions(-DGPROSHAN_EMBREE)
	include_directories(${embree_INCLUDE_DIRS})
endif(embree_FOUND)


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW 2 REQUIRED)
find_package(glfw3 3 REQUIRED)
find_package(glm REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(Boost COMPONENTS thread system)


set(THREADS_PREFER_PTHREAD_FLAG ON)


include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${GLFW3_INCLUDE_DIRS})
include_directories(${X11_INCLUDE_DIR})
include_directories(${AMADILLO_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${SuiteSparse_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

include_directories(${gproshan_SOURCE_DIR}/include)
include_directories(${gproshan_SOURCE_DIR}/imgui)


FILE(GLOB_RECURSE imgui_sources ${gproshan_SOURCE_DIR}/imgui/*.cpp)
add_library(imgui SHARED ${imgui_sources})


set(ptx_code "")
if(OptiX_INCLUDE)
	cuda_compile_ptx(ptx_files "${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu")
	list(GET ptx_files 0 ptx_file)
	set(ptx_code ptx_code.c)
	add_custom_command(	OUTPUT ${ptx_code} 
						COMMAND bin2c --const --padd 0 --type char --name ptx_code ${ptx_file} > ${ptx_code}
						DEPENDS ${ptx_file} 
						)
endif(OptiX_INCLUDE)


FILE(GLOB_RECURSE cpp_sources ${gproshan_SOURCE_DIR}/src/*.cpp)
if(CUDAToolkit_FOUND)
	FILE(GLOB_RECURSE cu_sources ${gproshan_SOURCE_DIR}/src/*.cu)
endif(CUDAToolkit_FOUND)

add_library(gproshan SHARED ${cpp_sources} ${cu_sources} ${ptx_code})

target_compile_options(gproshan PRIVATE -Wall -Wno-unused-result)

target_link_libraries(gproshan OpenMP::OpenMP_CXX)
target_link_libraries(gproshan OpenGL::GL)
target_link_libraries(gproshan GLEW::GLEW)
target_link_libraries(gproshan glfw)
target_link_libraries(gproshan ${X11_X11_LIB})
target_link_libraries(gproshan ${ARMADILLO_LIBRARIES})
target_link_libraries(gproshan ${SuiteSparse_LIBRARIES})
target_link_libraries(gproshan CGAL::CGAL)
target_link_libraries(gproshan ${OptiX_LIBRARY})
target_link_libraries(gproshan imgui)

if(embree_FOUND)
	target_link_libraries(gproshan embree)
endif(embree_FOUND)

if(CUDAToolkit_FOUND)
	target_link_libraries(gproshan CUDA::cudart)
	target_link_libraries(gproshan CUDA::cublas)
	target_link_libraries(gproshan CUDA::cusolver)
	target_link_libraries(gproshan CUDA::cusparse)
endif(CUDAToolkit_FOUND)

add_executable(gproshan_exe gproshan.cpp)
set_target_properties(gproshan_exe PROPERTIES OUTPUT_NAME gproshan)

add_executable(test_geodesics test_geodesics.cpp)
add_executable(test_image_denoising test_image_denoising.cpp)

target_link_libraries(gproshan_exe gproshan)
target_link_libraries(test_geodesics gproshan)
target_link_libraries(test_image_denoising gproshan)

file(MAKE_DIRECTORY tmp)

