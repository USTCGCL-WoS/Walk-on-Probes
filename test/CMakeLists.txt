add_executable(test_sampling test_sampling.cpp)
target_link_libraries(test_sampling PRIVATE core gtest_main)
add_test(NAME sampling COMMAND test_sampling)

add_executable(test_scalar test_scalar.cpp)
target_link_libraries(test_scalar PRIVATE core gtest_main)
add_test(NAME scalar COMMAND test_scalar)

add_executable(test_probe test_probe.cpp)
target_link_libraries(test_probe PRIVATE solver gtest_main)
add_test(NAME sh_probe COMMAND test_probe)