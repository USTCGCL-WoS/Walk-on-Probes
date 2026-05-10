function(ADD_LIB LIB_NAME LIB_TYPE)
    set(multiValueArgs PUBLIC_LIBS PRIVATE_LIBS INTERFACE_LIBS PUBLIC_DEFS PRIVATE_DEFS INTERFACE_DEFS)
    cmake_parse_arguments(ADD_LIB "" "" "${multiValueArgs}" ${ARGN})

    if(NOT LIB_TYPE MATCHES "^(SHARED|STATIC|INTERFACE)$")
        message(FATAL_ERROR "Invalid library type: ${LIB_TYPE}. Must be SHARED, STATIC, or INTERFACE")
    endif()

    set(current_dir ${CMAKE_CURRENT_SOURCE_DIR})

    if(LIB_TYPE STREQUAL "INTERFACE")
        add_library(${LIB_NAME} INTERFACE)
        
        if(EXISTS ${current_dir}/include)
            target_include_directories(${LIB_NAME}
                INTERFACE
                $<BUILD_INTERFACE:${current_dir}/include>
                $<INSTALL_INTERFACE:include>
            )
        endif()
    else()
        file(GLOB_RECURSE sources 
            ${current_dir}/src/*.cpp
            ${current_dir}/src/*.c
            ${current_dir}/src/*.cc
            ${current_dir}/src/*.cxx
        )
        
        file(GLOB_RECURSE headers
            ${current_dir}/include/*.h
            ${current_dir}/include/*.hpp
            ${current_dir}/src/*.h
            ${current_dir}/src/*.hpp
        )

        add_library(${LIB_NAME} ${LIB_TYPE} ${sources} ${headers})
        
        target_compile_features(${LIB_NAME} PUBLIC cxx_std_20)
        
        # Define BUILD_WOS_${LIB_NAME_UPPER}_MODULE for source files in this library, so they can use #ifdef to control symbol visibility
        string(TOUPPER ${LIB_NAME} LIB_NAME_UPPER)
        target_compile_definitions(${LIB_NAME} PRIVATE BUILD_WOS_${LIB_NAME_UPPER}_MODULE)
        
        if(ADD_LIB_PUBLIC_DEFS)
            target_compile_definitions(${LIB_NAME} PUBLIC ${ADD_LIB_PUBLIC_DEFS})
        endif()
        
        if(ADD_LIB_PRIVATE_DEFS)
            target_compile_definitions(${LIB_NAME} PRIVATE ${ADD_LIB_PRIVATE_DEFS})
        endif()
        
        if(ADD_LIB_INTERFACE_DEFS)
            target_compile_definitions(${LIB_NAME} INTERFACE ${ADD_LIB_INTERFACE_DEFS})
        endif()
        
        if(EXISTS ${current_dir}/include)
            target_include_directories(${LIB_NAME}
                PUBLIC
                $<BUILD_INTERFACE:${current_dir}/include>
                $<INSTALL_INTERFACE:include>
            )
        endif()
        
        if(EXISTS ${current_dir}/src)
            target_include_directories(${LIB_NAME}
                PRIVATE
                ${current_dir}/src
            )
        endif()
    endif()

    if(ADD_LIB_PUBLIC_LIBS)
        target_link_libraries(${LIB_NAME} PUBLIC ${ADD_LIB_PUBLIC_LIBS})
    endif()
    
    if(ADD_LIB_PRIVATE_LIBS)
        target_link_libraries(${LIB_NAME} PRIVATE ${ADD_LIB_PRIVATE_LIBS})
    endif()
    
    if(ADD_LIB_INTERFACE_LIBS)
        target_link_libraries(${LIB_NAME} INTERFACE ${ADD_LIB_INTERFACE_LIBS})
    endif()
endfunction()
