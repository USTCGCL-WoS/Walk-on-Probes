# Enrich compile_commands.json with header entries via compdb (if available)
find_program(COMPDB_PROGRAM compdb QUIET)
if(COMPDB_PROGRAM)
    add_custom_target(compdb ALL
        COMMAND "${COMPDB_PROGRAM}" -p "${CMAKE_BINARY_DIR}" list
                > "${CMAKE_BINARY_DIR}/compile_commands.json.tmp"
        COMMAND ${CMAKE_COMMAND} -E rename
                "${CMAKE_BINARY_DIR}/compile_commands.json.tmp"
                "${CMAKE_BINARY_DIR}/compile_commands.json"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Running compdb to add header entries to compile_commands.json"
    )
endif()