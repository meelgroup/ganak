include_directories( ${PROJECT_SOURCE_DIR} )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DENABLE_TBDD)
# add_definitions(-DENABLE_BTRACE)
add_definitions(-DCACHE_STATS)

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()


set(SOURCES
src/bddio.c
src/bddop.c
src/bvec.c
src/cache.c
src/fdd.c
src/imatrix.c
src/kernel.c
src/pairs.c
src/prime.c
src/reorder.c
src/tree.c
src/cppext.cxx
)

# set(buddy_public_headers bdd.h)

add_library(buddy ${SOURCES})

# set_target_properties(buddy PROPERTIES
#   #     PUBLIC_HEADER "${buddy_public_headers}"
#     VERSION 1.0
#     SOVERSION 1.0
# )

install(TARGETS buddy
    EXPORT ${GANAK_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ganak"
)
