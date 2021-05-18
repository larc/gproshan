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

find_package(CUDAToolkit 11.3)
if(CUDAToolkit_FOUND)
	enable_language(CUDA)
	set(CMAKE_CUDA_STANDARD 17)
	set(CMAKE_CUDA_STANDARD_REQUIRED ON)

	set(CMAKE_CUDA_FLAGS "-Xcompiler -fopenmp")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	set(CMAKE_CUDA_ARCHITECTURES 52 53 60 61 62 70 72 75 80)

	add_definitions(-DGPROSHAN_CUDA)
	include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})

	find_package(OptiX 7.3)
	if(OptiX_INCLUDE)
		add_definitions(-DGPROSHAN_OPTIX)
		include_directories(SYSTEM ${OptiX_INCLUDE})
	endif(OptiX_INCLUDE)
endif(CUDAToolkit_FOUND)


find_package(embree 3.12)
if(embree_FOUND)
	add_definitions(-DGPROSHAN_EMBREE)
	include_directories(SYSTEM ${embree_INCLUDE_DIRS})
endif(embree_FOUND)


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
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


include_directories(SYSTEM ${GLEW_INCLUDE_DIRS})
include_directories(SYSTEM ${GLFW3_INCLUDE_DIRS})
include_directories(SYSTEM ${X11_INCLUDE_DIR})
include_directories(SYSTEM ${AMADILLO_INCLUDE_DIR})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(SYSTEM ${SuiteSparse_INCLUDE_DIRS})
include_directories(SYSTEM ${CGAL_INCLUDE_DIRS})


include_directories(SYSTEM ${gproshan_SOURCE_DIR}/imgui)


FILE(GLOB_RECURSE imgui_sources ${gproshan_SOURCE_DIR}/imgui/*.cpp)
add_library(imgui SHARED ${imgui_sources})
target_link_libraries(imgui GLEW::GLEW)
target_link_libraries(imgui glfw)


include_directories(${gproshan_SOURCE_DIR}/include/gproshan)
add_compile_definitions(GPROSHAN_DIR="${gproshan_SOURCE_DIR}")


FILE(GLOB_RECURSE cpp_sources ${gproshan_SOURCE_DIR}/src/*.cpp)
if(CUDAToolkit_FOUND)
	FILE(GLOB_RECURSE cu_sources ${gproshan_SOURCE_DIR}/src/*.cu)
	list(REMOVE_ITEM cu_sources ${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu)
endif(CUDAToolkit_FOUND)


add_library(gproshan SHARED ${cpp_sources} ${cu_sources})

target_link_libraries(gproshan OpenMP::OpenMP_CXX)
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
	target_link_libraries(gproshan CUDA::cuda_driver)
	target_link_libraries(gproshan CUDA::cublas)
	target_link_libraries(gproshan CUDA::cusolver)
	target_link_libraries(gproshan CUDA::cusparse)
endif(CUDAToolkit_FOUND)


add_executable(test_gproshan gproshan.cpp)
set_target_properties(test_gproshan PROPERTIES OUTPUT_NAME gproshan)
target_link_libraries(test_gproshan gproshan)

if(OptiX_INCLUDE)
	add_library(ptx_sources OBJECT ${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu)
	set_property(TARGET ptx_sources PROPERTY CUDA_PTX_COMPILATION ON)
	add_custom_target(ptx_cp_sources DEPENDS ptx_sources COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_OBJECTS:ptx_sources> $<TARGET_FILE_DIR:gproshan>)
	add_dependencies(test_gproshan ptx_cp_sources)
endif(OptiX_INCLUDE)


add_executable(test_geodesics test_geodesics.cpp)
target_link_libraries(test_geodesics gproshan)


add_executable(test_image_denoising test_image_denoising.cpp)
target_link_libraries(test_image_denoising gproshan)


install(TARGETS gproshan imgui test_gproshan test_geodesics test_image_denoising
		ARCHIVE DESTINATION ${gproshan_SOURCE_DIR}/lib
		LIBRARY DESTINATION ${gproshan_SOURCE_DIR}/lib
		RUNTIME DESTINATION ${gproshan_SOURCE_DIR}/bin)


file(MAKE_DIRECTORY tmp)

