FILE(GLOB_RECURSE imgui_sources ${gproshan_SOURCE_DIR}/imgui/*.cpp)

add_library(imgui SHARED ${imgui_sources})
target_link_libraries(imgui GLEW::GLEW)
target_link_libraries(imgui glfw)

