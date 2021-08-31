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


file(MAKE_DIRECTORY ${gproshan_SOURCE_DIR}/tmp)

