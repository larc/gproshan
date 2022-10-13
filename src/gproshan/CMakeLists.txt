include_directories(${gproshan_SOURCE_DIR}/include/)

FILE(GLOB_RECURSE cpp_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(gproshan SHARED ${cpp_sources})

target_link_libraries(gproshan embree)
target_link_libraries(gproshan OpenMP::OpenMP_CXX)
target_link_libraries(gproshan GLEW::GLEW)
target_link_libraries(gproshan OpenGL::GL)
target_link_libraries(gproshan glfw)
target_link_libraries(gproshan ${X11_X11_LIB})
target_link_libraries(gproshan ${ARMADILLO_LIBRARIES})
target_link_libraries(gproshan ${SuiteSparse_LIBRARIES})
target_link_libraries(gproshan CGAL::CGAL)
target_link_libraries(gproshan ${OptiX_LIBRARY})
target_link_libraries(gproshan imgui)

if(CUDAToolkit_FOUND)
	FILE(GLOB_RECURSE cu_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)
	list(REMOVE_ITEM cu_sources ${CMAKE_CURRENT_SOURCE_DIR}/raytracing/rt_optix.cu)

	add_library(gproshan_cuda SHARED ${cu_sources})

	target_link_libraries(gproshan_cuda CUDA::cudart)
	target_link_libraries(gproshan_cuda CUDA::cuda_driver)
	target_link_libraries(gproshan_cuda CUDA::cublas)
	target_link_libraries(gproshan_cuda CUDA::cusolver)
	target_link_libraries(gproshan_cuda CUDA::cusparse)

	target_link_libraries(gproshan gproshan_cuda)
endif(CUDAToolkit_FOUND)


if(OptiX_INCLUDE)
	add_custom_target(	optix_ptx ALL
						COMMAND ${CMAKE_CUDA_COMPILER} --optix-ir --use_fast_math
								-isystem ${OptiX_INCLUDE}
								-isystem ${gproshan_SOURCE_DIR}/include
								-o ${gproshan_SOURCE_DIR}/src/rt_optix.ptx
								${gproshan_SOURCE_DIR}/src/gproshan/raytracing/rt_optix.cu
						DEPENDS ${gproshan_SOURCE_DIR}/src/gproshan/raytracing/rt_optix.cu
						)
endif(OptiX_INCLUDE)



set(installable_libs gproshan imgui)
if(TARGET gproshan_cuda)
	list(APPEND installable_libs gproshan_cuda)
endif()

install(TARGETS ${installable_libs}
		EXPORT gproshanTargets
		DESTINATION lib
		)

install(EXPORT gproshanTargets
		FILE gproshanTargets.cmake
		DESTINATION lib/cmake/gproshan
		)

set(INCLUDE_INSTALL_DIR "${gproshan_SOURCE_DIR}/include/")

include(CMakePackageConfigHelpers)

configure_package_config_file(
		${gproshan_SOURCE_DIR}/gproshanConfig.cmake.in
		"${CMAKE_CURRENT_BINARY_DIR}/gproshanConfig.cmake"
		INSTALL_DESTINATION "lib/cmake/gproshan"
		PATH_VARS INCLUDE_INSTALL_DIR
		)

write_basic_package_version_file(
		"${CMAKE_CURRENT_BINARY_DIR}/gproshanConfigVersion.cmake"
		VERSION "${gproshan_VERSION_MAJOR}.${gproshan_VERSION_MINOR}"
		COMPATIBILITY AnyNewerVersion
		)

install(FILES
		${CMAKE_CURRENT_BINARY_DIR}/gproshanConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/gproshanConfigVersion.cmake
		DESTINATION lib/cmake/gproshan
		)

export(	EXPORT gproshanTargets
		FILE "${CMAKE_CURRENT_BINARY_DIR}/gproshanTargets.cmake"
		)

