include_directories( ${PROJECT_SOURCE_DIR} )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_subdirectory (component_types)
add_subdirectory (clhash)

SET(SOURCES alt_component_analyzer.cpp
	component_cache.cpp
	component_management.cpp
	instance.cpp
	solver.cpp
	statistics.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp)

SET(LIB
	clhash
	component_types)

ADD_LIBRARY(sym_ganak ${SOURCES} ${LIB})
target_link_libraries (sym_ganak ${LIB})

if (NOT WIN32)
    set_target_properties(sym_ganak PROPERTIES
        PUBLIC_HEADER "${sym_ganak_public_headers}"
        VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
else()
    set_target_properties(sym_ganak PROPERTIES
        OUTPUT_NAME sym_ganak1win
        PUBLIC_HEADER "${sym_ganak_public_headers}"
        VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
        SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
endif()

add_executable (sym_ganak-bin
    main.cpp
)

set(sym_ganak_exec_link_libs
    ${Boost_LIBRARIES}
    ${GMP_LIB}
    ${GMPXX_LIB}
    ${NAUTY_LIB}
    ${IGRAPH_LIB}
)

target_link_libraries (sym_ganak-bin
sym_ganak ${sym_ganak_exec_link_libs}
)

set_target_properties(sym_ganak-bin PROPERTIES
    OUTPUT_NAME sym_ganak
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)
