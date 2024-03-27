FILE(GLOB_RECURSE imgui_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(${gproshan_SOURCE_DIR}/include/imgui/)

add_library(imgui SHARED ${imgui_sources})
target_link_libraries(imgui GLEW::GLEW)
target_link_libraries(imgui glfw)

