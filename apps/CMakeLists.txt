add_executable(test_gproshan gproshan.cpp)
set_target_properties(test_gproshan PROPERTIES OUTPUT_NAME gproshan)
target_link_libraries(test_gproshan gproshan)

if(OptiX_INCLUDE)
	add_library(ptx_sources OBJECT ${gproshan_SOURCE_DIR}/src/raytracing/rt_optix.cu)
	set_property(TARGET ptx_sources PROPERTY CUDA_PTX_COMPILATION ON)
	set_property(TARGET ptx_sources PROPERTY CUDA_ARCHITECTURES OFF)
	target_compile_options(ptx_sources PRIVATE "--keep")
	add_custom_target(	ptx_tmp DEPENDS ptx_sources
						COMMAND ${CMAKE_COMMAND} -E copy
						${gproshan_BINARY_DIR}/apps/rt_optix.ptx
						${gproshan_SOURCE_DIR}/tmp/rt_optix.ptx
						)
	add_dependencies(test_gproshan ptx_tmp)
endif(OptiX_INCLUDE)


add_executable(test_geodesics test_geodesics.cpp)
target_link_libraries(test_geodesics gproshan)


add_executable(test_image_denoising test_image_denoising.cpp)
target_link_libraries(test_image_denoising gproshan)


file(MAKE_DIRECTORY ${gproshan_SOURCE_DIR}/tmp)

