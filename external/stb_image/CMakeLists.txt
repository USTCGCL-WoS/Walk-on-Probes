add_library(stb_image_headers INTERFACE)
target_include_directories(stb_image_headers INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/stb
)

add_library(stb_image_impl OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/stb_image_impl.cpp)
target_link_libraries(stb_image_impl PUBLIC stb_image_headers)
