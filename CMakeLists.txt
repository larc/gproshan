cmake_minimum_required(VERSION 3.12)

project(gproshan VERSION 2.0)

list(APPEND CMAKE_MODULE_PATH "${gproshan_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CUDA 10.1)
if(CUDA_FOUND)
	enable_language(CUDA)
	add_definitions(-DGPROSHAN_CUDA)
	include_directories(${CUDA_INCLUDE_DIRS})
	
	find_package(OptiX 7.0)
	if(OptiX_INCLUDE)
		add_definitions(-DGPROSHAN_OPTIX)
		include_directories(${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDA_FOUND)


find_package(Embree 3.7)
if(EMBREE_LIBRARY)
	add_definitions(-DGPROSHAN_EMBREE)
	include_directories(${EMBREE_INCLUDE_DIRS})
endif(EMBREE_LIBRARY)


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW 2 REQUIRED)
find_package(glfw3 3 REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SuiteSparse REQUIRED)


set(THREADS_PREFER_PTHREAD_FLAG ON)

add_definitions(${EIGEN3_DEFINITIONS})

# ARMA warning ARMA_ALLOW_FAKE_GCC with g++-9
add_definitions(-DARMA_ALLOW_FAKE_GCC)

include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${GLFW3_INCLUDE_DIRS})
include_directories(${X11_INCLUDE_DIR})
include_directories(${AMADILLO_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${SuiteSparse_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

include_directories(${gproshan_SOURCE_DIR}/include)
include_directories(${gproshan_SOURCE_DIR}/include/viewer)
include_directories(${gproshan_SOURCE_DIR}/include/mdict)
include_directories(${gproshan_SOURCE_DIR}/include/cuda)
include_directories(${gproshan_SOURCE_DIR}/include/raytracing)
include_directories(${gproshan_SOURCE_DIR}/imgui)

FILE(GLOB_RECURSE imgui_sources ${gproshan_SOURCE_DIR}/imgui/*.cpp)
add_library(imgui STATIC ${imgui_sources})

set(ptx_code "")
if(CUDA_FOUND)
	cuda_compile_ptx(ptx_files "${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu")
	list(GET ptx_files 0 ptx_file)
	set(ptx_code ptx_code.c)
	add_custom_command(	OUTPUT ${ptx_code} 
						COMMAND bin2c --const --padd 0 --type char --name ptx_code ${ptx_file} > ${ptx_code}
						DEPENDS ${ptx_file} 
						)
endif(CUDA_FOUND)

FILE(GLOB_RECURSE cpp_sources ${gproshan_SOURCE_DIR}/src/*.cpp)
add_library(gproshan_cpp STATIC ${cpp_sources} ${ptx_code})

target_compile_options(gproshan_cpp PRIVATE -Wall -Wno-unused-result)

target_link_libraries(gproshan_cpp OpenMP::OpenMP_CXX)
target_link_libraries(gproshan_cpp OpenGL::GL)
target_link_libraries(gproshan_cpp GLEW::GLEW)
target_link_libraries(gproshan_cpp glfw)
target_link_libraries(gproshan_cpp ${X11_X11_LIB})
target_link_libraries(gproshan_cpp ${ARMADILLO_LIBRARIES})
target_link_libraries(gproshan_cpp ${SuiteSparse_LIBRARIES})
target_link_libraries(gproshan_cpp CGAL::CGAL)
target_link_libraries(gproshan_cpp ${EMBREE_LIBRARY})
target_link_libraries(gproshan_cpp ${OptiX_LIBRARY})
target_link_libraries(gproshan_cpp imgui)

if(CUDA_FOUND)
	FILE(GLOB_RECURSE cu_sources ${gproshan_SOURCE_DIR}/src/cuda/*.cu)
	add_library(gproshan_cu STATIC ${cu_sources})
	target_compile_options(gproshan_cu PRIVATE -Xcompiler -fopenmp)
	set_target_properties(gproshan_cu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	target_link_libraries(gproshan_cu gproshan_cpp)
	target_link_libraries(gproshan_cu ${CUDA_LIBRARIES})
	target_link_libraries(gproshan_cu ${CUDA_CUDA_LIBRARY})
	target_link_libraries(gproshan_cu ${CUDA_cublas_LIBRARY})
	target_link_libraries(gproshan_cu ${CUDA_cusolver_LIBRARY})
	target_link_libraries(gproshan_cu ${CUDA_cusparse_LIBRARY})
endif(CUDA_FOUND)

add_executable(gproshan gproshan.cpp)
add_executable(test_geodesics test_geodesics.cpp)
add_executable(test_image_denoising test_image_denoising.cpp)

target_link_libraries(gproshan gproshan_cpp)
target_link_libraries(test_geodesics gproshan_cpp)
target_link_libraries(test_image_denoising gproshan_cpp)

if(CUDA_FOUND)
	target_link_libraries(gproshan gproshan_cu)
	target_link_libraries(test_geodesics gproshan_cu)
	target_link_libraries(test_image_denoising gproshan_cu)
endif(CUDA_FOUND)

file(MAKE_DIRECTORY tmp)

