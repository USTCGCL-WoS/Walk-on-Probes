cmake_minimum_required(VERSION 3.15)
project(WoP VERSION 1.0.0 LANGUAGES CXX)

# C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# Define NOMINMAX to prevent Windows min/max macros
if(WIN32)
    add_compile_definitions(NOMINMAX)
endif()

# Export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Make external libraries SHARED by default
set(BUILD_SHARED_LIBS ON)

# MSVC specific settings
if(MSVC)
    # utf-8
    add_compile_options("/utf-8")
    
    # Add /PROFILE linker flag for RelWithDebInfo to enable profiling
    set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} /PROFILE")
    set(CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO} /PROFILE")
    set(CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO} /PROFILE")
endif()

# Options
option(WOS_BUILD_DEMO "Build the demo application" ON)
option(WOS_BUILD_TESTS "Build the test suite" OFF)
option(WOS_USE_FAST_MATH "Enable fast math optimizations" OFF)
option(WOS_USE_NATIVE_ARCH "Enable native CPU architecture optimizations (-march=native). Improves performance but reduces portability." ON)

# Set default build type for single-configuration generators
if(NOT CMAKE_CONFIGURATION_TYPES)
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
    endif()
endif()

# Set output directories to ensure DLLs and EXEs are in the same location
# Different build types go to different directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/$<CONFIG>)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/$<CONFIG>)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/$<CONFIG>)

# CTest
include(CTest)

# CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompDB)
include(AddLibrary)

# Find OpenMP (parallel processing support)
include(FindOpenMPCustom)

# External Libraries
add_subdirectory(external)

# Main Library
add_subdirectory(core)
add_subdirectory(scene)
add_subdirectory(solver)

# Demo Application
if(WOS_BUILD_DEMO)
    add_subdirectory(demo)
endif()

if(WOS_BUILD_TESTS)
    add_subdirectory(test)
endif()