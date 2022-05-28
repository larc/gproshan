include_directories(${gproshan_SOURCE_DIR}/include/)


add_executable(test_gproshan gproshan.cpp)
set_target_properties(test_gproshan PROPERTIES OUTPUT_NAME gproshan)
target_link_libraries(test_gproshan gproshan)


add_executable(test_geodesics test_geodesics.cpp)
target_link_libraries(test_geodesics gproshan)


add_executable(test_image_denoising test_image_denoising.cpp)
target_link_libraries(test_image_denoising gproshan)


file(MAKE_DIRECTORY ${gproshan_SOURCE_DIR}/tmp)

