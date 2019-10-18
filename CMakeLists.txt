cmake_minimum_required(VERSION 3.13)

project(gproshan VERSION 2.0 LANGUAGES CXX)

find_package(CUDA 10.1)
if(CUDA_FOUND)
	enable_language(CUDA)
endif(CUDA_FOUND)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SuiteSparse REQUIRED)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-fopenmp -Wall")
set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")

set(THREADS_PREFER_PTHREAD_FLAG ON)

add_definitions(${EIGEN3_DEFINITIONS})

# ARMA warning ARMA_ALLOW_FAKE_GCC with g++-9
add_definitions(-DARMA_ALLOW_FAKE_GCC)

if(CUDA_FOUND)
	include_directories(${CUDA_INCLUDE_DIRS})
endif(CUDA_FOUND)

include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(${X11_INCLUDE_DIR})
include_directories(${AMADILLO_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${SuiteSparse_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

include_directories(${gproshan_SOURCE_DIR}/include)
include_directories(${gproshan_SOURCE_DIR}/include/viewer)
include_directories(${gproshan_SOURCE_DIR}/include/mdict)
include_directories(${gproshan_SOURCE_DIR}/include/cuda)

FILE(GLOB_RECURSE cpp_sources ${gproshan_SOURCE_DIR}/src/*.cpp)
add_library(gproshan_cpp STATIC ${cpp_sources})

target_link_libraries(gproshan_cpp ${OpenMP_CXX_LIBRARIES})
target_link_libraries(gproshan_cpp ${OPENGL_LIBRARIES})
target_link_libraries(gproshan_cpp ${GLUT_LIBRARIES})
target_link_libraries(gproshan_cpp ${X11_X11_LIB})
target_link_libraries(gproshan_cpp ${ARMADILLO_LIBRARIES})
target_link_libraries(gproshan_cpp CGAL::CGAL)
target_link_libraries(gproshan_cpp ${SuiteSparse_LIBRARIES})

if(CUDA_FOUND)
	FILE(GLOB_RECURSE cu_sources ${gproshan_SOURCE_DIR}/src/*.cu)
	add_library(gproshan_cu STATIC ${cu_sources})
	set_target_properties(gproshan_cu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	target_link_libraries(gproshan_cu gproshan_cpp)
	target_link_libraries(gproshan_cu ${CUDA_LIBRARIES} ${CUDA_cublas_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY})
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

