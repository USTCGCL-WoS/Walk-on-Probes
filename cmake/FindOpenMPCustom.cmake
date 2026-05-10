# FindOpenMPCustom.cmake
# OpenMP finder with special handling for Clang on Windows
#
# This module attempts to find OpenMP.
# For Clang on Windows, it automatically searches LLVM installation for libomp.
#
# Output Variables:
#   OpenMP::OpenMP_CXX - Target will be available if OpenMP is found

# Special handling for Clang on Windows
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND WIN32)
    message(STATUS "Detecting Clang on Windows - searching for LLVM OpenMP...")
    
    # Try to locate libomp in LLVM installation directory
    get_filename_component(LLVM_BIN_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    get_filename_component(LLVM_ROOT_DIR ${LLVM_BIN_DIR} DIRECTORY)
    
    set(LIBOMP_INCLUDE_DIR "${LLVM_ROOT_DIR}/include")
    set(LIBOMP_LIBRARY "${LLVM_ROOT_DIR}/lib/libomp.lib")
    
    if(EXISTS ${LIBOMP_LIBRARY})
        # Manually configure OpenMP variables for CMake's FindOpenMP
        # For clang-cl, use /openmp. For clang++, use -fopenmp.
        if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
             set(OpenMP_CXX_FLAGS "/openmp")
        else()
             set(OpenMP_CXX_FLAGS "-fopenmp")
        endif()
        
        set(OpenMP_CXX_LIB_NAMES "libomp")
        set(OpenMP_libomp_LIBRARY ${LIBOMP_LIBRARY})
        
        message(STATUS "Found LLVM OpenMP library: ${LIBOMP_LIBRARY}")
    else()
        message(STATUS "LLVM OpenMP library not found at: ${LIBOMP_LIBRARY}")
    endif()
endif()

# Attempt to find OpenMP using standard CMake module
find_package(OpenMP QUIET)

# Error if OpenMP is not found
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP or ensure your compiler supports it.")
endif()