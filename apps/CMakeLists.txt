include_directories(${gproshan_SOURCE_DIR}/include/)


add_executable(test_gproshan gproshan.cpp)
set_target_properties(test_gproshan PROPERTIES OUTPUT_NAME gproshan)
target_link_libraries(test_gproshan gproshan)

if(OptiX_INCLUDE)
	find_package(CUDA)
	cuda_compile_ptx(ptx_sources ${gproshan_SOURCE_DIR}/src/gproshan/raytracing/rt_optix.cu)
	list(GET ptx_sources 0 ptx_file)
	add_custom_target(	ptx_tmp DEPENDS ${ptx_file}
						COMMAND ${CMAKE_COMMAND} -E copy
						${ptx_file}
						${gproshan_SOURCE_DIR}/tmp/rt_optix.ptx
						)
	add_dependencies(test_gproshan ptx_tmp)
endif(OptiX_INCLUDE)


add_executable(test_geodesics test_geodesics.cpp)
target_link_libraries(test_geodesics gproshan)


add_executable(test_image_denoising test_image_denoising.cpp)
target_link_libraries(test_image_denoising gproshan)


file(MAKE_DIRECTORY ${gproshan_SOURCE_DIR}/tmp)

