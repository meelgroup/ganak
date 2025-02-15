# Copyright (c) 2017, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/buddy-2.4/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/buddy-2.4/src/)
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})
include_directories(${ARJUN_INCLUDE_DIRS})
include_directories( ${GMP_INCLUDE_DIR} )
include_directories( ${MPFR_INCLUDE_DIR} )
include_directories( ${BREAKID_INCLUDE_DIRS} )
include_directories( ${APPROXMC_INCLUDE_DIRS} )

if (NOT WIN32)
    #add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    add_cxx_flag_if_supported("-Wduplicated-cond")
    add_cxx_flag_if_supported("-Wduplicated-branches")
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

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

SET(SOURCES
    comp_analyzer.cpp
    comp_cache.cpp
    comp_types/comp_archetype.cpp
    comp_management.cpp
    counter.cpp
    statistics.cpp
    clauseallocator.cpp
    IFlowCutter.cpp
    graph.cpp
    flow-cutter-pace17/src/tree_decomposition.cpp
    flow-cutter-pace17/src/cell.cpp
    flow-cutter-pace17/src/greedy_order.cpp
    TreeDecomposition.cpp
    clhash/clhash.cpp
    buddy-2.4/src/bddio.c
    buddy-2.4/src/bddop.c
    buddy-2.4/src/bvec.c
    buddy-2.4/src/cache.c
    buddy-2.4/src/fdd.c
    buddy-2.4/src/imatrix.c
    buddy-2.4/src/kernel.c
    buddy-2.4/src/pairs.c
    buddy-2.4/src/prime.c
    buddy-2.4/src/reorder.c
    buddy-2.4/src/tree.c
    buddy-2.4/src/cppext.cxx
    ganak.cpp
	${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

ADD_LIBRARY(ganak ${SOURCES} ${LIB})
target_link_libraries (ganak
  LINK_PUBLIC ${CRYPTOMINISAT5_LIBRARIES}
  LINK_PUBLIC ${ARJUN_LIBRARIES}
  LINK_PUBLIC ${MPFR_LIBRARIES}
  LINK_PUBLIC ${APPROXMC_LIBRARIES}
  LINK_PUBLIC ${BREAKID_LIBRARIES}
  LINK_PUBLIC ${GMP_LIBRARIES}
  )


set_target_properties(ganak PROPERTIES
    PUBLIC_HEADER "${ganak_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

add_executable (ganak-bin main.cpp)
add_executable (example example.cpp)

set(ganak_exec_link_libs
  ganak
)

IF (ZLIB_FOUND)
    SET(ganak_exec_link_libs ${ganak_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(ganak-bin ${ganak_exec_link_libs})
target_link_libraries(example ganak)

set_target_properties(ganak-bin PROPERTIES
  OUTPUT_NAME ganak
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  INSTALL_RPATH_USE_LINK_PATH TRUE)

ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/ganak.hpp )
ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/counter_config.hpp )
ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/mpreal.h )
ganak_add_public_header(ganak ${CMAKE_CURRENT_SOURCE_DIR}/lit.hpp )

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The ganakConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/ganak")
add_custom_target(CopyPublicHeaders ALL)
get_target_property(ganak_public_headers ganak PUBLIC_HEADER)
foreach(public_header ${ganak_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD
                       COMMAND ${CMAKE_COMMAND} -E make_directory
                               "${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E echo
                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E
                           copy_if_different
                           ${public_header}
                           "${HEADER_DEST}"
                      )
endforeach()


install(TARGETS ganak
    EXPORT ${GANAK_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ganak"
)

install(TARGETS ganak-bin
    EXPORT ${GANAK_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
