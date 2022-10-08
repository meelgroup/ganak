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

include_directories( ${PROJECT_SOURCE_DIR} )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_subdirectory (component_types)
add_subdirectory (clhash)

SET(SOURCES
    alt_component_analyzer.cpp
	component_cache.cpp
	component_management.cpp
	instance.cpp
	new_component_analyzer.cpp
	solver.cpp
	statistics.cpp
	${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp)

ADD_LIBRARY(ganak ${SOURCES} ${LIB})
target_link_libraries (ganak
  LINK_PUBLIC clhash
  LINK_PUBLIC component_types
  LINK_PUBLIC ${CRYPTOMINISAT5_LIBRARIES})

set_target_properties(ganak PROPERTIES
    PUBLIC_HEADER "${ganak_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

add_executable (ganak-bin
    main.cpp
)

target_link_libraries (ganak-bin
  ganak
  ${GMP_LIB} ${GMPXX_LIB}
)

set_target_properties(ganak-bin PROPERTIES
    OUTPUT_NAME ganak
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)

