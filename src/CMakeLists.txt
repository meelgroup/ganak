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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/buddy-2.4/src/)
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})
include_directories(${ARJUN_INCLUDE_DIRS})
include_directories(${LOUVAIN_COMMUNITIES_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories( ${GMP_INCLUDE_DIR} )

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

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_subdirectory (comp_types)
add_subdirectory (clhash)

SET(SOURCES
    alt_comp_analyzer.cpp
	comp_cache.cpp
	comp_management.cpp
	instance.cpp
	counter.cpp
	statistics.cpp
    clauseallocator.cpp
    IFlowCutter.cpp
    graph.cpp
    flow-cutter-pace17/src/tree_decomposition.cpp
    flow-cutter-pace17/src/cell.cpp
    flow-cutter-pace17/src/greedy_order.cpp
    TreeDecomposition.cpp
	${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp)

ADD_LIBRARY(ganak ${SOURCES} ${LIB})
target_link_libraries (ganak
  LINK_PUBLIC buddy
  LINK_PUBLIC clhash
  LINK_PUBLIC component_types
  LINK_PUBLIC ${CRYPTOMINISAT5_LIBRARIES}
  LINK_PUBLIC ${ARJUN_LIBRARIES}
  LINK_PUBLIC ${LOUVAIN_COMMUNITIES_LIBRARIES}
  LINK_PUBLIC ${Boost_LIBRARIES} # due to serialization
  LINK_PUBLIC ${GMP_LIBRARIES})
)

set_target_properties(ganak PROPERTIES
    PUBLIC_HEADER "${ganak_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

add_executable (ganak-bin main.cpp)

set(ganak_exec_link_libs
  ganak
  ${Boost_LIBRARIES} # due to program_options
  ${GMP_LIBRARIES}
)

IF (ZLIB_FOUND)
    SET(ganak_exec_link_libs ${ganak_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(ganak-bin ${ganak_exec_link_libs})

if (NOT WIN32)
    set_target_properties(ganak-bin PROPERTIES
    OUTPUT_NAME ganak
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
else()
    set_target_properties(ganak-bin PROPERTIES
    OUTPUT_NAME ganak-win
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    )
endif()
