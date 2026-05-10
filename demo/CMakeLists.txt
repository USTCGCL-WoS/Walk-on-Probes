add_executable(demo main.cpp)
target_link_libraries(demo core scene solver stb_image_impl CLI11::CLI11)

add_executable(debug_wost debug_wost.cpp)
target_link_libraries(debug_wost core scene solver)
