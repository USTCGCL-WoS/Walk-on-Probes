# CLI11
add_subdirectory(cli11)

# Eigen
add_subdirectory(eigen)

# FCPW
add_subdirectory(fcpw)
# Workaround for FCPW CMake bug: FCPW_EIGEN_INCLUDES and FCPW_ENOKI_INCLUDES
# are set via PARENT_SCOPE but not locally, so FCPW_DEPS_INCLUDES ends up empty.
# Manually add them to the fcpw INTERFACE target so they propagate to consumers.
target_include_directories(fcpw INTERFACE ${FCPW_EIGEN_INCLUDES} ${FCPW_ENOKI_INCLUDES})

# nlohmann/json
set(JSON_BuildTests OFF CACHE INTERNAL "" FORCE)
set(JSON_Install OFF CACHE INTERNAL "" FORCE)
add_subdirectory(nlohmann_json)

# spdlog
add_subdirectory(spdlog)

# tinyxml2
set(tinyxml2_BUILD_TESTING OFF CACHE INTERNAL "" FORCE)
add_subdirectory(tinyxml2)

# pcg_cpp
add_subdirectory(pcg_cpp)

# stb_image
add_subdirectory(stb_image)

# Google Test (built static to avoid DLL export issues with Clang on Windows)
set(BUILD_SHARED_LIBS_SAVE "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(gtest)
set(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS_SAVE}")

# nanoflann
set(NANOFLANN_BUILD_EXAMPLES OFF CACHE INTERNAL "" FORCE)
set(NANOFLANN_BUILD_TESTS OFF CACHE INTERNAL "" FORCE)
add_subdirectory(nanoflann)